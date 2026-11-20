#!/usr/bin/env python3
"""Per-column comparison of VEP output files in every output format.

Three modes, one JSON report with a ``verdict`` each, exit 0 on ``green`` and 1 on ``red``:

``--format {default,tab,vcf,json,parquet,recoder} --expected F --actual F``
    One Ensembl VEP output (expected) against one candidate output (actual) in one format.
    Every row of the default and tab formats, every CSQ entry and record line of the VCF
    format, and every consequence object and record object of the JSON format is one
    comparable entry; every field of every entry is compared verbatim. A Parquet actual (a
    file or a directory of part files) is read back through
    ``scripts/adapters/parquet_to_vep_tab.py``, which needs ``pyarrow``, and compared with a
    tab expected as a multiset, both sides in the adapter's canonical form (a frequency or
    score cell as a number, so ``0.0010`` and ``0.001`` agree; every other cell verbatim).
    ``--input`` (the input VCF) keys the VCF and JSON entries by
    input record (a line's CHROM, POS, ID, REF and ALT, which both engines write back
    unchanged; a JSON object's ``input`` names its record the same way, since Ensembl VEP
    writes a structural variant's empty or ``0`` QUAL, FILTER or INFO back as ``.``), so a
    record one side skipped keeps every later record aligned; beside ``--manifest`` its file
    name selects the flag set documented over that input, in the default and tab formats
    too, which read nothing else of it. ``recoder`` compares two
    ``variant_recoder`` outputs (``vep recode``'s is the candidate's): the JSON array's
    objects are matched by their ``input`` text in order, each allele of an object is one
    entry, and every key of the allele is compared as a sorted set of its values (Ensembl VEP
    writes the array's values and the object's keys in hash order), a ``var_synonyms`` entry
    with its identifiers sorted, a ``mane_select`` triple and a ``ga4gh_vrs`` object as
    canonical JSON; a ``null`` object and an object of warnings alone are entries of their
    own. The built-in classes ``refseq_transcript_entries`` and ``lrg_transcript_entries``
    drop, on both sides, the ``hgvsc`` and ``hgvsp`` values and the ``mane_select`` triples
    naming a RefSeq (``NM_``, ``NR_``, ``XM_``, ``XR_``, ``NP_``, ``XP_``) or LRG model, and the
    ``hgvsg``, ``spdi`` and ``vcf_string`` values on an LRG region, which the recoder reads from
    its merged database and an Ensembl cache does not carry;
    ``archive_dbsnp_synonyms`` drops the ``var_synonyms`` entries of the ``Archive dbSNP``
    source, the merged identifiers the cache's synonyms column never holds; and
    ``identifier_indel_transcript_notations`` drops the ``hgvsc``, ``hgvsp`` and
    ``mane_select`` keys of an allele of an identifier input that is an insertion or a
    deletion, whose transcript notations the recoder leaves unshifted (the stored class of the
    database row, ``indel`` for a non-minimal allele string, survives its minimisation, and
    ``_return_3prime`` shifts only an insertion or a deletion) where every other input kind
    and ``vep recode`` shift them 3' along the transcript.

``--self-parity --default F --tab F --vcf F --json F [--parquet P]``
    One engine's formats against each other, the default format as the pivot: every record,
    every (allele, feature) entry and every column the two formats both declare must carry
    the same value once the format's own encoding is undone.

``--warnings --expected F --actual F [--expected-skipped F --actual-skipped F]``
    Two engines' ``_warnings.txt`` files, and optionally their ``--skipped_variants_file``
    outputs, as sets of (input line number, description) with the quoted input line and any
    path masked.

Masks, floors and field pairs (comparison mode):

* ``--manifest`` is a golden-corpus ``manifest.json`` whose ``divergences``,
  ``vep_rs_only_tuples`` and ``field_divergences`` name keys documented to differ, each
  under an ``expected_divergence`` class; ``--masks CLASS,...`` names the classes applied.
  A documented consequence set or field value of an applied class is checked against the
  actual instead of the expected (an actual value other than the documented one is a stale
  divergence, red); a key documented as candidate-only is set aside, and a record the
  reference skipped keeps its record-level line and object uncompared; a key of a class not
  named is compared like any other. An entry that names a second class under
  ``second_expected_divergence`` documents a value two causes change and is applied only
  where both classes are named: with either alone the value is compared like any other,
  whichever of the runs read documents it. The three lists of the manifest's own run and those of
  each of its ``flag_sets`` entries are read alike, so a run under a flag set's flags applies
  the keys that set documents once their class is named. ``--run NAME`` names the one run the
  pair belongs to, ``default`` for the manifest's own run or a flag set's name: the lists of
  that run alone are read (a flag set without ``divergences`` of its own reads the
  manifest's), so an item of another run over the same input neither is expected on the
  pair's rows nor sets one aside, and a record-level row the run documents (the one row of a
  ``--summary`` run) has nothing to observe on a VCF line or a JSON object, which print no
  consequence set of their own. ``--reference-sequence present|absent`` states whether the
  candidate's run had a reference sequence. A ``divergences`` item may state the candidate's
  sets for such a run under ``vep_rs_consequence_sets_with_reference_sequence``, beside
  ``vep_rs_consequence_sets``, its sets in a run without one: under ``present`` the item
  documents the former, and where they are Ensembl VEP's own set it documents nothing, is
  not counted as applied, and its row is compared like any other with no class accepting a
  difference of its Consequence; under ``absent`` and without the option the item documents
  ``vep_rs_consequence_sets``, and an item without the key reads the same under all three.
  An item is counted as applied as it is read, except in a comparison none of whose entries
  names a Location or an input record the lists number (the VCF and JSON formats of an
  input that is not a VCF, or that no run of the manifest names, whose record ordinals the
  manifest's own lists do not count): no item meets a row there, and a row counts once,
  where a class's own predicate takes it. The class
  ``motif_score_reverse_strand_mnv`` is built in as well as documented: under it a differing
  ``MOTIF_SCORE_CHANGE`` on a reverse-strand ``MotifFeature`` entry whose allele is two or
  more bases is accepted whether or not a manifest documents the key. The class
  ``regulatory_feature_rows`` is built in rather
  than read from a manifest: it drops the expected ``RegulatoryFeature`` and ``MotifFeature``
  rows in every format for a candidate cache without regulatory data, recomputes the
  ``most_severe_consequence`` of each JSON record that lost a row from the rows that remain,
  and in ``--warnings`` mode requires the candidate's notice about missing regulatory data
  (or the notice of releases 0.1.0 to 0.5.0) and leaves it out of the comparison. The class
  ``order_only`` is built in too: the set-valued columns
  ``DOMAINS``, ``FREQS``, ``ZYG`` and ``NearestExonJB`` (and the ``Extra`` text carrying them)
  agree when their items agree as a multiset in any order, since Ensembl VEP prints the
  protein features in the order of the cache's array and the frequencies, the per-sample
  zygosities and the plugin's exon items in hash order, and the ``5UTR_annotation`` groups
  agree when each group's ``key=value`` pairs agree as a multiset; every such agreement is
  counted per column as ``order_only`` in the report, never passed as equal. The class
  ``individual_allele_order`` is built in as well: Ensembl VEP
  lists the two alternate alleles of an individual's copy of a record in hash order
  (``Parser/VCF.pm:639``), so an entry carrying ``IND`` (``ind`` on a JSON record object)
  agrees when the copy's alternate rows come in any order within each feature (a reference
  row its number or its reference allele identifies keeping its place), when its
  ``allele_string`` and ``UPLOADED_ALLELE`` name
  the same reference allele and the same alternate alleles in any order, and when its
  ``ALLELE_NUM`` values differ by one permutation of the alternates' numbers applied to every
  row of the copy, the reference allele keeping 0; every such agreement is counted per column
  as ``individual_allele_order``, the reordered entries under ``counts``, and an entry
  without ``IND`` is compared like any other. A consistent renumbering whose allele strings
  do not move is accepted too, since ``UPLOADED_ALLELE`` prints the record's own string
  wherever an original allele string is kept. Two built-in classes rewrite the reference's
  text before it is parsed (``reference_text_masked``): ``conservation_short_name_spelling``
  spells the key of a ``--phyloP`` source ``phyloP<N>way``, as the candidate and the file it
  reads do, where Ensembl VEP's option set spells it ``phlyoP<N>way`` (``Config.pm:555``),
  and ``pc_description_text`` gives the header line of a custom source's ``PC`` field the
  candidate's description of the share the field reports, in place of the two texts of
  ``AnnotationSourceAdaptor.pm:395-397``, each of which names the share the other setting of
  ``reciprocal`` computes (``File.pm:556-564``). The two registry pair classes of the
  SNP/indel comparator, ``start_cooccurrence_swap`` and ``splice_family_swap``, are built in
  so that a run without a manifest sets aside what the concordance figures set aside: a key
  both engines write whose consequence sets differ by the shape the registry's excluding rule
  admits (``compare_vep_outputs.py`` ``EXCLUDING_RULES``) is skipped whole, counted under the
  class and under ``counts.registry_pairs_set_aside``, and the JSON record's
  ``most_severe_consequence`` is left uncompared as for a documented key. The field classes
  of ``FIELD_CLASSES`` (``gxf_reverse_strand_utr_sequence``, ``conservation_max_text_order``,
  ``gwas_haplotype_cross_product``, ``alphamissense_two_digit_transcript_version``,
  ``cadd_structural_allele_match``, ``splicevault_sibling_allele_rows``,
  ``splicevault_json_later_rows``, ``utrannotator_gained_start_window``,
  ``utrannotator_consequence_precedence``, ``utrannotator_structural_allele_text``,
  ``utrannotator_lost_start_after_insertion``,
  ``utrannotator_gained_start_remade_in_place``,
  ``utrannotator_frameshift_at_stop_boundary``,
  ``utrannotator_lost_stop_written_before_stop_codon``,
  ``loeuf_zero_minimum``,
  ``maxentscan_overlapping_windows``, ``alphamissense_sibling_allele_rows``,
  ``revel_sibling_allele_rows``, ``revel_other_reference_amino_acid_rows``,
  ``gnomadmt_info_value_second_equals``,
  ``loftee_anc_allele_multiallelic_record``, ``loftee_first_alt_donor_filter``,
  ``nearest_exon_jb_stretch_priming``, ``nearest_exon_jb_shifted_span``,
  ``nearest_exon_jb_id_less_exon_fold``, ``downstream_reverse_strand_rotation``,
  ``spliceai_chopped_multiallelic``, ``spliceai_multi_entry_cell``,
  ``plugin_contig_mitochondrial_spelling``,
  ``loftee_gerp_zoom_level_summary``, ``splice_family_swap`` on the Downstream plugin's fields,
  ``gwas_interaction_first_snp_allele_whitespace``,
  ``refseq_edited_model_hgvsc_offset``, ``refseq_edited_model_shifted_insertion_bases``,
  ``refseq_edited_model_seam_insertion_hgvsc``,
  ``refseq_edited_model_unminimised_allele_hgvsc``,
  ``refseq_edited_model_hgvs_offset``,
  ``hgvs_shifted_insertion_clamped_end``, ``hgvs_shifted_indel_written_span_gate``,
  ``hgvsp_shifted_indel_written_span_gate``,
  ``hgvs_shifted_deletion_clamped_span``, ``hgvs_clipped_reference_against_whole_alternate``,
  ``hgvs_added_bases_left_at_the_span_end``,
  ``hgvs_insertion_left_where_the_shift_stops``,
  ``hgvs_added_bases_written_outside_the_transcript``,
  ``hgvs_intron_offset_abutting_insertion``,
  ``hgvsp_stop_codon_protein_comparison``,
  ``hgvsp_stop_named_on_incomplete_coding_sequence``,
  ``hgvsp_intron_end_insertion``, ``hgvsp_intron_start_insertion``,
  ``shift_intron_term_unshifted_span``,
  ``downstream_utr_frame_past_coding_end``, ``downstream_start_past_own_codon``,
  ``refseq_edited_model_deleted_span_codons``,
  ``utrannotator_annotation_pair_order``, ``nearest_exon_jb_item_order``,
  ``frequency_filter_zero_value``, ``plugin_hgvs_shifted_cds_coordinates``,
  ``plugin_value_follows_splice_term_past_the_record_span``,
  ``plugin_value_follows_stop_codon_term_pair``) each accept
  one shape
  of difference in the fields of one feature or plugin, the constants' comments state the
  shapes (``hgvsp_intron_end_insertion`` and ``hgvsp_intron_start_insertion`` read the exons and
  coding bounds of the row's transcript from the candidate JSON cache named by ``--json-cache``
  and are refused without it;
  ``frequency_filter_zero_value`` also sets aside the rows and record entry of a
  record the candidate alone keeps whose every ``FREQS`` item is 0;
  ``utrannotator_consequence_precedence`` sets aside a ``5UTR_consequence`` whose two values
  are both among the plugin's five terms and differ, with the row's ``5UTR_annotation``, where
  the runs of ``--utrannotator-sweep`` write both terms on the row;
  ``utrannotator_structural_allele_text`` sets aside the ``5UTR_consequence`` and
  ``5UTR_annotation`` the reference alone fills on a row whose Allele is a structural class
  word, the reference's term one of the plugin's five, and those the candidate fills on the
  default or tab row of a symbolic deletion whose span lies inside one exon of the 5' UTR
  and of a plain symbolic insertion that carries the frameshift its stated length decides,
  where the reference's differ from them or are empty;
  ``utrannotator_lost_start_after_insertion`` sets aside a ``5UTR_consequence`` that is the
  lost-start term on the reference's row alone, with the row's ``5UTR_annotation``, where the
  annotation's one ``CapDistanceToStart`` is the lower cDNA position of an insertion's row, or
  lies one or two below it on a row whose allele, read on the transcript's strand, ends with
  the start codon's first base or its first two;
  ``utrannotator_gained_start_remade_in_place`` sets aside the two fields of an insertion's
  row on which the reference reports a gained start one or two below the lower cDNA position,
  itself or under ``--utrannotator-sweep``, where the allele so read opens with the start
  codon's last two bases or its last and the candidate reports no gained start there;
  ``utrannotator_frameshift_at_stop_boundary`` sets aside the ``5UTR_annotation`` of a row on
  which both report the frameshift and the candidate lists, beyond the reference's groups,
  only uORFs whose length changes and whose stop codon begins directly after the insertion
  point the row's HGVSc and ``HGVS_OFFSET`` give, so it applies to runs under ``--hgvs``;
  ``utrannotator_lost_stop_written_before_stop_codon`` sets aside the two fields of an
  insertion's row on which the reference reports the frameshift and the candidate the lost
  stop, each of the candidate's groups carrying the Kozak pairs of one uORF the reference
  lists and the alternative stop codon that uORF's printed ``alt_type_length`` gives; the
  two ``--shift_3prime 1`` classes read a row's ``SHIFT_LENGTH``, so they apply to runs under
  ``--shift_length``, and the two shifted-indel notation classes read it where the run prints
  it and a non-zero ``HGVS_OFFSET`` in any case, so they apply under ``--hgvs`` alone too;
  ``refseq_edited_model_shifted_insertion_bases`` reads the row's non-zero ``HGVS_OFFSET``, its
  ``Allele`` and its ``STRAND``, so it applies wherever both engines print the allele as written:
  every row under ``--hgvs`` alone and the non-coding rows under ``--shift_3prime 1``);
  an accepted entry is counted once under the class whatever the format, so the count agrees
  across formats (the classes that read a record's alleles, the three sibling-allele classes,
  the two ``LoF_filter`` classes, ``spliceai_chopped_multiallelic`` and the two classes of a
  later allele's row, ``hgvs_shift_shared_across_alleles`` and
  ``plugin_hgvs_shifted_cds_coordinates``, read them from the VCF
  and JSON entries and from a default or tab row whose ``Uploaded_variation`` spells them, and
  leave a row under an identifier a difference), a field two classes accept is the first's in
  ``FIELD_CLASSES`` order, the
  report lists it per column under ``masked`` and the first entries under
  ``mask_examples``. The row class ``merged_same_stable_id_collapse`` accepts a transcript row
  the candidate alone writes under ``--refseq`` or ``--merged``, of a model the reference wrote
  no row for on a record it wrote a transcript row for: ``AnnotationType/Transcript.pm:313-343``
  keeps one transcript per stable identifier among the models a buffer loads, so a model the
  cache carries at two loci loses its rows at one of them; the row's ``SOURCE``, where the run
  writes one, is ``RefSeq``, a repeat of an accepted row is a duplicate, and the JSON record's
  ``most_severe_consequence`` is left uncompared as for a documented key. The row class
  ``refseq_edited_model_deleted_span_term`` accepts the Consequence of a transcript row of a
  sequence allele on a model edited to the genome (``BAM_EDIT=OK``) that the reference
  writes as ``intergenic_variant`` and the candidate as exactly one of
  ``coding_sequence_variant``, ``5_prime_UTR_variant`` and ``3_prime_UTR_variant``, neither
  row with a cDNA position and both with one IMPACT: the allele lies wholly inside bases the
  model's alignment deleted, which map to no transcript coordinate, so no predicate fires in
  Ensembl VEP and ``BaseVariationFeatureOverlapAllele.pm:285`` assigns the default term;
  every other field of the row compares exactly, the row is counted once whatever the
  format, and the JSON record's ``most_severe_consequence`` is left uncompared as for a
  documented key. The row class ``refseq_edited_model_tail_term`` accepts, in the same way,
  the Consequence of such a row that the reference writes as ``intergenic_variant`` and the
  candidate as ``3_prime_UTR_variant``, both with a cDNA position and one IMPACT: the allele
  lies on one of the model's final bases, whose edited cDNA coordinate exceeds the unedited
  exon total that ``within_cdna`` (``Utils/VariationEffect.pm:680``) reads from
  ``Transcript::length``, so the UTR predicate is false. The row classes
  ``mature_mirna_unshifted_span`` and ``stop_codon_unshifted_span`` accept, in the same way,
  the Consequence of a transcript row with a positive ``SHIFT_LENGTH`` whose two sets are
  equal but for ``mature_miRNA_variant`` against ``non_coding_transcript_exon_variant`` in
  either direction, or for the reference's ``coding_sequence_variant`` against the candidate's
  ``stop_retained_variant``: under ``--shift_3prime 1`` the predicates read the shifted span
  from the variation feature's offset coordinates, but ``within_mature_miRNA``
  (``Utils/VariationEffect.pm:545-549``) and ``_overlaps_stop_codon_cil`` (``:1403``,
  ``:1424-1427``) overlap their feature with the sequence-region coordinates
  ``BaseVariationFeature.pm:279-282`` and ``:304-307`` cached before the shift, so those two
  tests alone read the unshifted span. The row classes of the stop codon
  tests behind ``stop_retained`` (``Utils/VariationEffect.pm:1290-1324``) accept, in the same
  way, the Consequence of an insertion's or deletion's transcript row, each by the shape the
  comment above its predicate states. On a row one of them takes, the cells
  ``term_pair_following_cell`` reads from the two sets of terms are counted under it as well:
  an IMPACT that is each side's most severe term's, a field of the Downstream, NMD or LOFTEE
  plugin written exactly on the side whose terms let the plugin run, and an HGVSp that is a
  frameshift's notation exactly on the side without ``stop_retained_variant``. They are
  ``stop_retained_without_stop_codon_test`` (the reference's ``stop_retained_variant``
  against nothing, ``coding_sequence_variant`` or ``frameshift_variant`` on a row that shows
  the stop codon test is not made: an allele letter that names no base, a cDNA position
  missing at one end, or no CDS start);
  ``stop_codon_overlap_genomic_span`` (the reference's ``coding_sequence_variant`` against
  the candidate's ``stop_retained_variant``, both beside ``splice_region_variant``, on a row
  with both cDNA positions and a CDS start);
  ``stop_retained_three_base_stop_codon_test`` (on a row whose test is made, the
  reference's ``stop_retained_variant`` against ``coding_sequence_variant``, or it, alone or
  beside an in-frame or protein-altering term, against ``frameshift_variant`` with or without
  ``stop_gained``);
  ``stop_codon_removed_without_stop_codon_test`` (the reference's
  ``stop_retained_variant`` against the candidate's ``stop_lost`` on a row without a cDNA
  position at one end or a CDS start);
  ``stop_codon_overlap_insertion_range`` (an insertion's row with both cDNA positions and
  a CDS start: the reference's ``coding_sequence_variant`` away from ``splice_region_variant``,
  or its ``frameshift_variant`` with or without ``stop_gained``, against the candidate's
  ``stop_retained_variant``, the second with or without ``inframe_insertion``);
  a row that ``stop_codon_unshifted_span`` takes as well is that class's, which stands before
  them in ``TERM_PAIR_CLASSES``.
  The row class ``splice_term_past_the_record_span`` accepts, in the same way and with the
  same cells, the Consequence of a transcript row. On a record whose alternate allele is
  longer than its reference allele, Ensembl VEP writes a splice donor, splice acceptor, donor
  fifth base, donor region, polypyrimidine tract or intron term from positions past the end
  of the record's span; where no base of the span lies in a window of that term, vep-rs does
  not write it. The six terms are ``splice_donor_variant``, ``splice_acceptor_variant``,
  ``splice_donor_5th_base_variant``, ``splice_donor_region_variant``,
  ``splice_polypyrimidine_tract_variant`` and ``intron_variant``, whose window is the intron,
  its end bases included. The class reads the record where an entry carries the reference
  bases its allele stands in for (``Entry.reference``): a default or tab row whose
  ``Uploaded_variation`` spells the alleles, a VCF entry by its line's POS, REF and ALT, a
  JSON entry by its object's ``start``, ``end`` and ``allele_string``, and a default or tab
  row under an identifier by the records of ``--input``, read by position
  (``read_entry_references``). With the record a row has the shape where the record is no
  pure insertion, the reference bases are two or more and as many as the Location's
  positions, the allele is longer, the cDNA and CDS positions are the same on both sides,
  and one intron can stand beside the span, at a place the row's cDNA positions allow, such
  that the record's differing regions give the reference's terms of the six and its
  ``splice_region_variant`` and the span gives the candidate's: the span may lie inside an
  exon, across its end or in an intron and may be two positions, and the donor region and
  ``splice_region_variant`` come where the terms that withhold them went. Without the record
  a row has the shape where its Location holds three or more positions, its allele more
  nucleotide letters than the span has bases and its two cDNA positions, the same on both
  sides as the CDS positions are, numbers as far apart as the span's two ends, so that the
  span holds no base of an intron: the candidate's row holds none of the six, with or
  without ``splice_region_variant`` where a donor, acceptor, donor-region or fifth-base term
  went, and a Location of two neighbouring positions, which an insertion's row shares, and a
  span that holds an intron's base stay differences. A pure insertion is taken in neither
  reading. A row that differs by the terms of a stop codon class as well is taken where that
  class's rule takes the rest of the difference; the stop codon classes stand before this
  one in ``TERM_PAIR_CLASSES``. A VCF or JSON entry is read at the Location its record gives
  its allele (``Entry.allele_location``), so the class takes a row in every format.
  The field class
  ``refseq_edited_model_deleted_span_codons`` accepts the Codons of a ``BAM_EDIT=OK``
  transcript row whose USED_REF has more bases than its CDS span, the alternate codons equal
  and the reference codons differing in an upper-case base alone, Ensembl VEP's upper-case
  bases spelling USED_REF's leading bases on the transcript's strand:
  ``TranscriptVariationAllele.pm:872-874`` splice the genomic reference allele over the
  shorter span before reading the reference codons.
  The field class ``hgnc_id_copy`` accepts an ``HGNC_ID`` one side leaves empty and the
  other fills with an id on a transcript row whose ``SYMBOL_SOURCE`` is not HGNC, in either
  direction: ``AnnotationType/Transcript.pm:290`` and ``:305`` lend a same-symbol gene's id
  to the transcripts of one buffer, so the field follows the buffer's other records; a run
  writes it under ``--symbol`` or in the VCF format, whose option set implies the flag. The
  field class ``organelle_peptide_paths_translation_table`` is built in as well as documented:
  it accepts a differing ``HGVSp`` on a row whose transcript the candidate JSON cache named by
  ``--json-cache`` files under a codon table that reads some codon unlike the standard one,
  when both values name one protein and the row's consequence set and the two descriptions
  place the difference in a frameshift's description, a lost stop's extension count or the
  choice between a duplication and an insertion: ``TranscriptVariationAllele.pm:2333``,
  ``:2450`` and ``:2492`` translate those sequences with the standard code against a
  reference peptide of the transcript's table. The class changes the HGVSp alone, so the two
  rows state one coding change (one consequence set and the same cDNA, CDS and protein
  positions, amino acids, codons, ``HGVSc`` and ``HGVS_OFFSET``), and an inserted
  allele of two or more bases on a reverse-strand transcript, shifted past its own length by
  a length that is no multiple of it, stays a difference: ``:1337-1339`` leave such an allele
  unrotated, so its HGVSp is read from bases the record does not produce. The rows do
  not state the table, so without the cache the class applies from a manifest alone, and a
  row a manifest documents keeps its recorded value. The
  field classes ``downstream_utr_frame_past_coding_end`` and
  ``downstream_start_past_own_codon`` accept the ``DownstreamProtein`` of the Downstream
  plugin on a frameshift row of a deletion, with or without inserted bases, whose first
  coding base is a codon's third base: ``Downstream.pm:112-115`` start the downstream
  translation at the codon after the first one the change touches, the candidate at that
  codon. The rows do not state the sequences, so the classes read the candidate JSON cache
  named by ``--json-cache`` and recompute both values from the transcript's coding sequence
  and 3' UTR under the transcript's codon table, with the ``CDS_position`` as printed, the
  record's own coordinates. The first takes
  a deletion without inserted bases
  that runs to the last coding base, where Ensembl VEP reads past the changed sequence's end
  and translates the 3' UTR from its first base (``:115-119``), the second every other such
  row, so the two share no row. Each changes the ``DownstreamProtein`` alone: both entries
  carry the value, one consequence set or two that differ as ``start_cooccurrence_swap``
  admits, one ``CDS_position`` (or the two of a later allele's row that
  ``plugin_hgvs_shifted_cds_coordinates`` reads, below) and one
  ``ProteinLengthChange``, a transcript absent from the cache or under a codon table
  ``crates/vep-core/src/codon.rs`` lacks stays a difference, and without the cache the
  classes accept nothing. The field class ``plugin_hgvs_shifted_cds_coordinates`` accepts a
  ``DownstreamProtein``, a ``ProteinLengthChange`` or an ``NMD`` on a row whose reference
  entry carries a non-zero ``HGVS_OFFSET``: under ``--hgvs`` the protein notation leaves the
  row's CDS coordinates and allele at the 3'-shifted position
  (``TranscriptVariationAllele.pm:844-858``, ``:2406-2410``, ``:1330-1345``) and Ensembl
  VEP's plugins read them there (``Downstream.pm:108-110``, ``NMD.pm:88-103``), where the
  candidate reads the record's own with and without the flag. The class reads the same
  cache with its mapper pairs, maps the record's span moved by the offset to the shifted
  coordinates, and accepts a Downstream value where the candidate's is its rule on the own
  coordinates, the reference's is Ensembl VEP's rule on the shifted ones and that rule gives
  another value on the own ones, and an ``NMD`` one side alone carries where Ensembl VEP's
  rule holds on the shifted coordinates exactly when the reference carries it and on the own
  ones exactly when the candidate does. On the row of a second or later alternate allele of
  its record the reference may print the shifted pair itself as its ``CDS_position``, the
  pair an earlier allele's notation left on the object the record's alleles share
  (``OutputFactory.pm:1737``): the own coordinates are then the candidate's ``CDS_position``,
  the record's span, and the reference's that span moved by the row's ``HGVS_OFFSET`` or by
  that of the last earlier allele's row that carries one (``later_allele_prints_shifted_pair``,
  ``earlier_allele_shift``); two ``CDS_position`` values on any other row stay a
  difference. The plugin reads the pair at the row's own offset whatever the row prints, and
  on such a row without an offset the Downstream plugin reads the printed pair, the one the
  earlier allele's notation left, with the allele as written: a Downstream value there is
  the class's where it is Ensembl VEP's rule on that pair. A row whose shifted and own
  readings give Ensembl
  VEP one value is left to the two Downstream classes, so the three share no row. A row
  whose candidate value is not Ensembl VEP's on the own coordinates either is changed by the
  start rule as well as the flag: it is this class's only where the Downstream class of its
  own coordinates is named beside it, and stays a difference with either named alone. A row
  without the offset that is no later allele's of a record whose earlier allele's row
  carries one, an ``NMD`` on any row without it, a shift that leaves the coding sequence and
  a value planted on either side stay differences. A text entry's Location is the record's span; a VCF or JSON entry's
  is the one its CDS coordinates map back to, and a row whose coordinates name no one span
  has no class. The field class ``plugin_value_follows_stop_codon_term_pair`` accepts a
  field of Downstream, NMD or LOFTEE that one side alone writes on a row whose two
  Consequence sets differ as a stop codon class admits
  (``STOP_CODON_PLUGIN_VALUE_PAIRS``): the plugins run on a row's terms, so the field stands
  on the side whose terms let the plugin run (``term_pair_following_cell``). The field class
  ``plugin_value_follows_splice_term_past_the_record_span`` accepts such a field, by the same
  rule, on a row ``splice_term_past_the_record_span`` takes, with the record where an entry
  carries it and without it where none does: an NMD or LOFTEE field on the
  reference's side where the donor or acceptor term it alone holds is the row's only
  truncating term, a Downstream field on the candidate's side where the reference holds
  ``frameshift_variant`` beside that term and the candidate's row holds no splice term; the
  intron term decides none of those fields. The field class
  ``blosum62_hgvs_shift_shared_across_alleles`` accepts a ``BLOSUM62`` the reference alone
  writes on a later allele's row that ``hgvs_shift_shared_across_alleles`` takes by its CDS
  and protein positions, where the reference's ``Amino_acids`` is one residue over one, read
  at the moved span, the candidate's own pair is none the plugin scores, and the value is
  the matrix cell of the reference's pair. The field class
  ``plugin_value_on_reference_skipped_record`` names a DosageSensitivity or UTRAnnotator
  field the candidate alone writes on a row of a record the reference dropped: such a row
  has no reference entry, so this comparison sets it aside whole under
  ``reference_skipped_record`` and takes no pair for the class, which the plugin tuple
  scorer decides from the rows of the value's Location. The
  three classes of the structural-variant comparator (``scripts/validation/
  compare_sv_concordance.py``) are built in through its own predicates over the entries read
  as its tuples, so a structural-variant case sets aside what the adjusted structural-variant
  figures set aside: ``reference_skipped_record`` (every candidate-only row of a record
  Ensembl VEP drops before annotation, ``perl_skip_route``, and wrote nothing for; it reads
  the input VCF, so ``--input`` is needed in every format), ``breakend_mate_local_read`` (a
  bracket-mate pair whose consequence sets are the local-coordinate and the mate-coordinate
  derivations over the transcript model, skipped whole as a registry pair is) and
  ``bnd_own_chromosome_rows_lost`` (a candidate-only own-end transcript row at a Location
  where the reference kept another chromosome's mate rows and no own-end row, the transcript
  within 5,000 bases on the Location's chromosome); the two breakend classes read the
  candidate JSON cache named by ``--json-cache``. A VCF or JSON entry's Location, which those
  formats do not print, is the one its input record gives (``vep_rs_record_rows``).
  In ``--warnings`` mode ``cadd_structural_allele_match`` and
  ``dosage_sensitivity_regulatory_row_warnings`` set aside the reference's ``Plugin '<name>'
  went wrong`` lines, and ``fork_warning_pid_prefix`` removes the ``<pid> : `` Ensembl VEP's
  parent wrote before a forked child's line. The classes of ``MANIFEST_CLASSES`` have no
  predicate over the two outputs and apply from a manifest; two of them are also decided
  without a manifest. ``hgvs_shift_shared_across_alleles`` is on three columns of a later
  allele's row when ``--json-cache`` names the candidate's cache: a differing
  ``CDS_position``, ``Protein_position`` or ``Amino_acids`` on the row of a second or later
  alternate allele of its record whose first allele's reference row on the feature carries a
  non-zero ``HGVS_OFFSET``, the candidate's coordinates the record's own and the reference's
  those moved by that offset, the two ``Amino_acids`` sharing their alternate part
  (``later_allele_shared_shift``). ``hgvs_shift_sequence_edge_window`` is when
  ``--reference-fasta`` names the reference the run read and ``--json-cache`` the candidate's
  cache, since the sequence and the transcript's span state what the outputs do not: a
  differing ``HGVSc`` or ``HGVS_OFFSET`` on an insertion's or deletion's row within 1,000
  bases of its region's edge on the walked side, each side's two values those its own walk
  leaves on the row's transcript, the candidate's walk over the allele's own flank and the
  reference's over the clipped window (``SequenceEdges``). A manifest's candidate-only tuples
  under a class
  of ``RECORD_DROPPING_CLASSES`` mark their record as one the reference wrote nothing for, so
  its VCF line and JSON object are counted, not compared, and a key a manifest documents
  under a built-in class is the manifest's. A requested class that accepts
  nothing is listed under ``masks.unused``, unless a value set aside under another class
  needs it named beside that class (the second class of a value two causes change): such a
  class is used, and the value stays counted under its own class alone.
* ``--floors COL=MIN,...``: the populated (present and not ``-``) count of column ``COL``
  must reach ``MIN`` on both sides. ``COL`` is the default-format name in every format; in
  the JSON format the count is taken under the column's JSON key (lower-cased and renamed
  as the JSON factory renames it, SIFT and PolyPhen under either of their two keys, a
  position under either of its start and end keys), or, for a column outside Ensembl VEP's
  own (``CORE_COLUMNS``) when no entry of either side carries that key, under the key of the
  same name inside a plugin's object (``alphamissense`` holds ``am_class``; a plugin's own copy
  of a core column, SpliceAI's ``SYMBOL``, never stands in for the column), and a column with
  no row count there has no floor there: a record-level column, which the JSON format
  carries once per record or not at all, ``Feature_type``, which names the consequence list
  and is no key of its entries, ``Extra``, the text formats' container of the keyed fields,
  and a custom source's columns (``<source>``, ``<source>_<field>``, ``<source>_FILTER``),
  which the JSON format carries once per record under ``custom_annotations``
  (``JSON.pm:249-258``); the report lists such a column under ``floors_not_applicable`` and
  the text formats of the same case check its floor. In the tab and VCF formats a floor on a
  column neither header declares (the tab header line and the CSQ ``Format:`` list carry the
  fields the run's options name, where the default format's ``Extra`` carries every key a
  row has, a VCF source's ``FILTER`` included) is not applicable there in the same way; a
  column one header declares and the other does not stays a difference the field comparison
  reports. Under ``conservation_short_name_spelling`` a floor spelt with the misspelt key is
  read under the renamed one.
* ``--field-pairs 'COL=EXPECTED=>ACTUAL;...'`` (repeatable): a value pair accepted as
  agreeing wherever it occurs in ``COL``.
* ``--gwas-catalog <file>``: the GWAS catalog associations file the run's GWAS plugin read; the
  two GWAS classes read its rows of two or more SNPs and accept nothing without it.
* ``--utrannotator-sweep <file>`` (repeatable): a default-format Ensembl VEP output of the
  same input under a fixed Perl hash seed (``PERL_HASH_SEED`` with ``PERL_PERTURB_KEYS=0``),
  one per seed; ``utrannotator_consequence_precedence`` reads the terms those runs write on a
  row and accepts nothing without them, and ``utrannotator_gained_start_remade_in_place``
  reads them on a row whose reference writes another term than the gained start.
* ``--nearest-sites <json_cache>``: ``NEAREST`` values (and the JSON ``nearest``) agree
  through the candidate cache's protein-coding transcription start sites when every expected
  value is among the actual values and every actual value starts at a site of an expected
  value. Ensembl VEP names one transcript per winning site, whichever its interval tree held
  first, and the candidate names every transcript at the site, so the ids differ where the
  sites agree.
* ``--keep-only regulatory_feature_rows`` in comparison mode inverts the built-in row mask:
  both inputs retain only ``RegulatoryFeature`` and ``MotifFeature`` consequence rows,
  plus the VCF/JSON record entries of records retaining such a row. A JSON record that
  loses non-regulatory rows has its ``most_severe_consequence`` recomputed from its retained
  rows. Other record fields compare verbatim. This cannot be combined with
  ``--masks regulatory_feature_rows`` or used in warnings or self-parity mode.

A comparison is green when no expected entry is missing, no actual entry is unexpected
(the two one-sided key counts), every compared field agrees after masks and field pairs,
the shared entries come in the same order (not for Parquet, whose rows the adapter reads in
file order, sorted by chromosome and position within each row group, and
not under ``--unordered``, for a run whose reference rows Ensembl VEP emits in hash order:
the allele buckets of ``--pick_allele`` and ``--pick_allele_gene`` and the gene buckets of
``--per_gene``, ``OutputFactory.pm:593``, ``:605`` and ``:855``; a golden corpus marks such
a flag set ``unordered``), no documented divergence is stale and every floor is met.

With the environment variable ``VEP_FORMATS_REMAINDER`` naming the executable of
``tools/vep-formats-remainder`` (a path, or a command on ``PATH``), a ``default`` or ``tab``
comparison is handed over to that front end: it reads both outputs on every core, keeps
back the record groups (every row of a side that shares a Location) that are the same text
in the same order on both sides, and this comparator compares the rest, the remainder, and
adds the kept-back rows to the counts, the populated columns and the place where the order
of the shared keys differs. The comparator hands it the Locations the manifest documents and
what marks a reference row that ``regulatory_feature_rows`` drops or
``conservation_short_name_spelling`` rewrites; a group holding one is in the remainder
whole. So is a group with a row that holds a byte outside ASCII, a run of more than 640
digits or more than 64 opening brackets: this comparator converts some values before it
compares them (a whole number, an object's text), and a value the interpreter does not
convert ends the run with this comparator's error, which such a row therefore meets on
the whole outputs and on the remainder alike. The
report and the standard output are the same bytes with and without the
variable. A comparison naming a structural-variant class or ``individual_allele_order``, one
under ``--keep-only``, every other format and the other two modes compare the whole
outputs; one line on the standard error, written once the report is, says which was done.
Unset or empty, nothing is handed over and that line is not written. A value naming no
executable, a front end that fails and a hand-over file that cannot be written end the run
with one error line and exit 2, before the report file is opened.

The parsers follow the comparators of ``crates/vep-cli/tests/common/mod.rs``. The output
grammar they read is Ensembl VEP's at release/116.2: ``OutputFactory/VEP_output.pm:99-136``
(the thirteen fixed columns and the ``Extra`` column of ``KEY=value`` pairs joined by ``;``,
``;`` in a value as ``%3B``), ``OutputFactory/Tab.pm:125-129`` and ``:187-198`` (every field
a column, ``-`` for an undefined value), ``OutputFactory/VCF.pm:100-119`` and ``:379-420``
(the CSQ fields in the ``##INFO`` header's ``Format:`` list, lists joined with ``&``, ``,``
and ``|`` rewritten to ``&``, ``;`` to ``%3B``, whitespace to ``_``, an undefined value
empty), ``OutputFactory/JSON.pm:93-136`` and ``:229-389`` (keys lower-cased, ``-`` values
dropped, ``YES`` as 1, ``*_position`` split into ``*_start`` and ``*_end``, SIFT and
PolyPhen split into prediction and score, ``DOMAINS`` into ``db``/``name`` objects, the
renames of ``%RENAME_KEYS``, the consequence lists named by feature type, numbers
numberified), ``OutputFactory.pm:883-931`` with ``Utils.pm:141-166`` (``Uploaded_variation``
and ``Location``), ``BaseVEP.pm:802-863`` and ``BaseRunner.pm:364-380`` (the warnings and
skipped-variants files), and ``%OVERLAP_CONSEQUENCES`` of ensembl-variation
``Utils/Constants.pm`` at release/116 (the term ranks behind ``most_severe_consequence``).
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import zlib
from collections import Counter, defaultdict, deque
from collections.abc import Callable, Collection, Iterator
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import NamedTuple, TypeVar

# What this comparator shares with the tuple comparator is defined in that script, once
# for both.
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_vep_outputs import (
    ABSENT,
    CODING_SEQUENCE_TERM,
    FRAMESHIFT_TERM,
    STOP_GAINED_TERM,
    STOP_RETAINED_TERM,
    ComparatorError,
    allele_is_indel,
    blank,
    consequence_set,
    dash,
    feature_kind,
    indel_row,
    insertion_row,
    integer_of,
    location_span,
    perl_split,
    read_record_references,
    reference_allele_is_empty,
    replaced_reference,
    sequence_location,
    uploaded_alleles,
    uploaded_position,
    vep_rs_pos,
)
from compare_vep_outputs import (
    splice_term_past_the_record_span_pair as _splice_term_past_the_record_span_pair,
)
from compare_vep_outputs import (
    stop_codon_overlap_genomic_span_pair as _stop_codon_overlap_genomic_span_pair,
)
from compare_vep_outputs import (
    stop_codon_overlap_insertion_range_pair as _stop_codon_overlap_insertion_range_pair,
)
from compare_vep_outputs import (
    stop_codon_removed_without_stop_codon_test_pair as _stop_codon_removed_without_stop_codon_test_pair,
)
from compare_vep_outputs import (
    stop_codon_test_coordinates as _stop_codon_test_coordinates,
)
from compare_vep_outputs import (
    stop_codon_test_not_made as _stop_codon_test_not_made,
)
from compare_vep_outputs import (
    stop_retained_three_base_stop_codon_test_pair as _stop_retained_three_base_stop_codon_test_pair,
)
from compare_vep_outputs import (
    stop_retained_without_stop_codon_test_pair as _stop_retained_without_stop_codon_test_pair,
)

COMPONENT = "compare_vep_formats"
ADAPTER = Path(__file__).resolve().parent.parent / "adapters" / "parquet_to_vep_tab.py"
# The annotation formats, and with the recoder's JSON every format the comparison mode takes.
VEP_FORMATS = ("default", "tab", "vcf", "json", "parquet")
FORMATS = (*VEP_FORMATS, "recoder")

# The class a vep_rs_only_tuples entry carries when every record it names is one the
# reference dropped before annotation.
REFERENCE_SKIPPED_RECORD = "reference_skipped_record"
# The built-in mask class dropping the reference's regulatory rows (``without_regulatory_rows``).
REGULATORY_FEATURE_ROWS = "regulatory_feature_rows"
REGULATORY_KINDS = {"regulatory_feature", "motif_feature"}
# The built-in field class ``motif_score_reverse_strand_mnv``: MOTIF_SCORE_CHANGE on a
# reverse-strand MotifFeature entry whose allele is two or more bases, which Ensembl VEP
# writes from the motif position of the variant's genomic start, the allele's last position
# along a reverse-strand motif (MotifFeatureVariationAllele.pm:245), so its score, or its
# absence, differs from the score of the allele written over the positions it covers. An
# insertion or a length-changing allele carries no score on either side, so a difference on
# such an allele is this class alone.
MOTIF_SCORE_REVERSE_STRAND_MNV = "motif_score_reverse_strand_mnv"
MOTIF_SCORE_FIELDS = ("MOTIF_SCORE_CHANGE", "motif_score_change")
STRAND_FIELDS = ("STRAND", "strand")
# A motif score change is the difference of two relative similarity scores, each in [0, 1],
# printed to three decimals; the JSON number drops trailing zeros.
MOTIF_SCORE_TEXT = re.compile(r"-?[01](?:\.[0-9]{1,3})?")
# The input a corpus manifest's own run annotates; a flag set names any other.
CORPUS_INPUT = "variants.vcf"
# The key under which a manifest entry names the class of a second cause that changes the value
# it documents: the entry is applied only where both of its classes are named.
SECOND_EXPECTED_DIVERGENCE = "second_expected_divergence"
# The name ``--run`` gives the manifest's own run, as the golden harness names it; a flag set
# is named by its ``name``, which the harness holds to another word.
OWN_RUN = "default"
# The optional key of a ``divergences`` item: the candidate's consequence sets on the key in
# a run that has a reference sequence, ``vep_rs_consequence_sets`` being those of a run
# without one.
WITH_REFERENCE_SETS = "vep_rs_consequence_sets_with_reference_sequence"
# The values of ``--reference-sequence``: whether the candidate's run had a reference sequence.
REFERENCE_PRESENT, REFERENCE_ABSENT = "present", "absent"
# The 0.1.0-to-0.5.0 and the cache-aware warnings-file phrases for --regulatory on a cache
# without regulatory data.
REGULATORY_NOTICE = "regulatory-feature annotation is not implemented"
MISSING_REGULATORY_NOTICE = "this cache holds no regulatory data"
EXAMPLES_PER_GROUP = 4
# The built-in classes of the recoder format: the prefixes each drops from the values of both
# sides (the hgvsc and hgvsp values and the mane_select triples of a RefSeq model; those and
# the hgvsg, spdi and vcf_string values of an LRG region), since variant_recoder reads the
# merged transcript set and LRG mappings from its database (VariantRecoder.pm sets `merged`
# and `lrg`) and an Ensembl cache carries neither.
RECODER_MASKS = {
    "refseq_transcript_entries": ("NM_", "NR_", "XM_", "XR_", "NP_", "XP_"),
    "lrg_transcript_entries": ("LRG_",),
}
RECODER_TRANSCRIPT_KEYS = ("hgvsc", "hgvsp")
# The keys whose values name a transcript, protein or region by prefix, the ones the classes read.
RECODER_NOTATION_KEYS = (*RECODER_TRANSCRIPT_KEYS, "hgvsg", "spdi", "vcf_string")
# The recoder format's two other built-in classes: the var_synonyms source the cache's column
# never carries (the merged dbSNP identifiers), and the keys the recoder writes unshifted for an
# identifier input of an insertion or deletion (TranscriptVariationAllele.pm `_return_3prime`
# shifts only a feature whose class is insertion or deletion, and the database row's stored
# class, `indel` for a non-minimal allele string, survives Parser.pm `minimise_alleles`).
ARCHIVE_DBSNP_SYNONYMS = "archive_dbsnp_synonyms"
ARCHIVE_DBSNP_SOURCE = "Archive dbSNP: "
IDENTIFIER_INDEL_NOTATIONS = "identifier_indel_transcript_notations"
IDENTIFIER_INDEL_KEYS = ("hgvsc", "hgvsp", "mane_select")
INDEL_NOTATION = re.compile(r"del|dup|ins")

# OutputFactory/JSON.pm:98-112 at release/116.2.
JSON_RENAME = {
    "consequence": "consequence_terms",
    "gene": "gene_id",
    "allele": "variant_allele",
    "symbol": "gene_symbol",
    "symbol_source": "gene_symbol_source",
    "overlapbp": "bp_overlap",
    "overlappc": "percentage_overlap",
    "refseq": "refseq_transcript_ids",
    "ensp": "protein_id",
    "chr": "seq_region_name",
    "variation_name": "id",
    "sv": "colocated_structural_variants",
    "clinical_impact": "clinvar_somatic_classification",
}
JSON_CONSEQUENCE_LISTS = (
    ("transcript_consequences", "transcript"),
    ("intergenic_consequences", "intergenic"),
    ("regulatory_feature_consequences", "regulatory_feature"),
    ("motif_feature_consequences", "motif_feature"),
)
# The engine identity vep-rs writes in every JSON object and Ensembl VEP never writes, set
# aside like ``input`` rather than compared; its header forms (``## vep-rs <version>``,
# ``##vep-rs=<version>``) fall with the other ``##`` lines the text parsers skip.
JSON_IDENTITY_KEY = "vep_rs_version"
# Columns the default format prints on every row of a record while the JSON format
# carries them once per record, in colocated_variants or at the top level, or not at all:
# every key VariationFeature_to_output_hash sets (OutputFactory.pm:883-1023: SV, VARIANT_CLASS,
# IND, ZYG, MINIMISED, AMBIGUITY, NEAREST and CHECK_REF) is copied once to the record object
# (JSON.pm:262), never into a consequence entry; the JSON factory's add_colocated_variant_info
# returns its hash unchanged (JSON.pm:493-495), so FREQS (OutputFactory.pm:1149-1163) has no
# JSON key, and JSON.pm:419-426 copies the AF keys alone out of the per-allele hash that
# carries MAX_AF and MAX_AF_POPS.
RECORD_LEVEL_COLUMNS = {
    "Uploaded_variation",
    "Location",
    "Existing_variation",
    "SV",
    "VARIANT_CLASS",
    "IND",
    "ZYG",
    "MINIMISED",
    "AMBIGUITY",
    "NEAREST",
    "CHECK_REF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "VAR_SYNONYMS",
    "FREQS",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLINVAR_SOMATIC_CLASSIFICATION",
}
RECORD_LEVEL_COLUMN_RE = re.compile(r"(^|_)AF$")
# The built-in class under which a set-valued column agrees when its items agree as a
# multiset: Ensembl VEP prints DOMAINS in the order of the cache's protein_features array
# (ensembl-variation BaseTranscriptVariation.pm:634-643 at release/116 copies it from the
# core API's Translation.pm:772-819, which buckets a translation's features by analysis and
# flattens the buckets in Perl hash order, the database's rows inside each), FREQS in hash
# order (OutputFactory.pm:1152-1153) and the ZYG list of --individual_zyg in the hash order
# of its samples (OutputFactory.pm:983-987), so none of the orders carries meaning.
ORDER_ONLY = "order_only"
# The NearestExonJB plugin lists its exon items in the order of a hash's keys as well
# (NearestExonJB.pm:267-271).
ORDER_ONLY_COLUMNS = ("DOMAINS", "FREQS", "ZYG", "NearestExonJB")
# The UTRAnnotator plugin writes its ``5UTR_annotation`` as groups of ``key=value`` pairs
# joined by ``:``, the pairs in hash order within each group (UTRAnnotator.pm:282-288); under
# ``order_only`` the groups agree as a multiset of groups whose pairs agree as a multiset.
ORDER_ONLY_PAIR_GROUP_COLUMN = "5UTR_annotation"
# The marker ``order_only_fields`` gives a JSON array field in place of a list separator, and
# the prefix it gives the pair-group column's separator.
JSON_ARRAY = "json"
PAIR_GROUPS = "pairs:"
# The JSON field of each order-only column with its separator: the ``domains`` array of a
# consequence entry, the record entry's ``zyg`` array and the consequence's ``nearestexonjb``
# string, whose items ``json_scalar`` joins with ``,``; FREQS has no JSON key.
ORDER_ONLY_JSON_FIELDS = {"domains": JSON_ARRAY, "zyg": ",", "nearestexonjb": ","}
# The built-in class under which an individual's alleles agree as a set. Under --individual,
# Parser/VCF.pm:635-639 at release/116.2 collects the alternate alleles an individual's
# genotype names in a hash and writes the copy's allele string as the reference allele
# followed by the hash's keys, so a copy carrying two alternate alleles lists them in hash
# order, and that order runs through the copy's rows within each feature, its ALLELE_NUM
# values (the allele's index in the string, BaseVariationFeatureOverlapAllele.pm:217-218 at
# release/116) and the JSON record's allele_string (JSON.pm:241). UPLOADED_ALLELE follows it
# only on a record with no original allele string (OutputFactory.pm:1360 prefers that
# string, which minimise_alleles sets on an indel record and under --minimal, Parser.pm:852,
# 881, 978-982).
INDIVIDUAL_ALLELE_ORDER = "individual_allele_order"
# The record-level token of an individual's copy: IND on a default, tab or CSQ entry
# (OutputFactory.pm:971), ind on a JSON record object.
INDIVIDUAL_FIELDS = ("IND", "ind")
# The fields carrying the copy's allele string, and the allele's index in it.
ALLELE_STRING_FIELDS = ("allele_string", "UPLOADED_ALLELE", "uploaded_allele")
ALLELE_NUMBER_FIELDS = ("ALLELE_NUM", "allele_num")
# The fields carrying the reference allele of a --show_ref_allele row.
REFERENCE_ALLELE_FIELDS = ("REF_ALLELE", "ref_allele")
# The built-in class spelling the reference's conservation key as the candidate does. Each
# --phyloP <N> expands into a --custom source whose short name is ``phlyoP<N>way`` while the
# file it names is spelt ``phyloP<N>way`` (Config.pm:555 at release/116.2), so the misspelt key
# heads the reference's header line, tab column, Extra item, CSQ field and JSON
# ``custom_annotations`` member where the candidate writes ``phyloP<N>way``; the class renames
# the key throughout the reference text before it is parsed.
CONSERVATION_SHORT_NAME_SPELLING = "conservation_short_name_spelling"
MISSPELT_CONSERVATION_KEY = re.compile(r"phlyoP(\d+)way")
CONSERVATION_KEY_SPELLING = r"phyloP\1way"
# The built-in class mapping the reference's two descriptions of a custom source's ``PC`` field
# onto the shares the field reports. AnnotationSourceAdaptor.pm:395-397 at release/116.2 names
# the field under reciprocal=0 as the share of the input variant covered and under reciprocal=1
# as a minimum reciprocal overlap, while File.pm:556-564 computes the share of the annotation
# record the variant covers without reciprocal and the share of the variant the record covers
# with it. The candidate's header names the share computed; the class puts that text into the
# reference's ``## <key>_PC : <text> from <file>`` and ``##INFO=<ID=<key>_PC,...>`` lines.
PC_DESCRIPTION_TEXT = "pc_description_text"
PC_DESCRIPTIONS = {
    "Percentage of input variant covered by reference variant": (
        "Percentage of the annotation record covered by the input variant"
    ),
    "Percentage of minimum reciprocal overlap between input variant and reference variant": (
        "Percentage of the input variant covered by the annotation record"
    ),
}
PC_FIELD_SUFFIX = "_PC"
# The two pair classes of the SNP/indel comparator's registry (compare_vep_outputs.py
# EXCLUSION_REGISTRY, whose predicates decide them here as well), built in so that a
# population case, which has no manifest beside its input, sets aside the pairs the
# concordance figures set aside: a key both engines write whose consequence sets differ by
# the one shape each rule admits. ``start_cooccurrence_swap`` is Ensembl VEP's emission of
# ``start_lost`` and ``start_retained_variant`` together on one allele (Utils/VariationEffect.pm
# :851-881, :947-963 and :1028-1075), ``splice_family_swap`` its loss of ``splice_region_variant``,
# with ``intron_variant`` or alone, to the assignment at BaseTranscriptVariationAllele.pm:215.
# A matched pair is skipped whole, as a manifest's documented key is: its terms, IMPACT and
# HGVSp follow the term the engines disagree on.
START_COOCCURRENCE_SWAP = "start_cooccurrence_swap"
SPLICE_FAMILY_SWAP = "splice_family_swap"
REGISTRY_CLASSES = (START_COOCCURRENCE_SWAP, SPLICE_FAMILY_SWAP)
# The built-in field classes of the corrected behaviours, each one shape of disagreement in
# the fields of one feature or plugin, applied without a manifest; every other field of the
# entry compares exactly. A value is empty when it is ``-``, the empty string or absent.
#
# ``gxf_reverse_strand_utr_sequence``: on a transcript Ensembl VEP reads from a GTF or GFF
# file through --gtf or --gff, BaseTranscriptVariation.pm:1101-1113 (``_generic_utr``)
# transfers the transcript onto its feature slice and reads the UTR from the copy; Exon.pm:985
# drops the exons' cached bases on the transfer, and FastaSequence.pm:369-371 maps a
# reverse-strand slice's relative coordinates from its genomic start, where the core API counts
# them from its end, so a reverse-strand transcript's 3' UTR is read as the reverse complement
# of the span mirrored about the slice. A cache transcript carries its UTR from the dumper and
# never reaches the read. The stop codon a frameshift or extension runs into lies in that UTR,
# so the count after ``fsTer`` or ``extTer`` of the HGVSp differs on a ``STRAND`` -1 row and
# nothing else does.
GXF_REVERSE_STRAND_UTR_SEQUENCE = "gxf_reverse_strand_utr_sequence"
HGVSP_FIELDS = ("HGVSp", "hgvsp")
TER_COUNT_RE = re.compile(r"^(.*(?:fs|ext)Ter)(\d+|\?)$")
# ``organelle_peptide_paths_translation_table``: TranscriptVariationAllele.pm:2333, :2450,
# :2492 and :2555 translate a frameshift's new frame, the reference an inserted peptide is
# tested against, the sequence a stop distance is counted in and an empty-window deletion's
# peptide with ``translate()`` bare, the standard code, against a reference peptide of the
# transcript's codon table; the candidate translates both sides under that table. The rows do
# not state the table, so the class reads it from the candidate JSON cache
# (``OrganelleCodonTables``, ``--json-cache``). The shape: an HGVSp on a row whose transcript
# the cache files under a table that reads some codon unlike the standard one, both values
# naming one protein, the row's consequence set and the two descriptions placing the
# difference in one of the paths, and nothing else of the row's coding change differing
# (``organelle_peptide_path``). Without the cache the class applies from a manifest alone.
ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE = "organelle_peptide_paths_translation_table"
# The genetic codes that read all 64 codons as the standard code does: the standard code and
# table 11, which differs from it in its start codons alone.
STANDARD_READING_CODON_TABLES = (1, 11)
# ``conservation_max_text_order``: Conservation.pm:242-243 sorts a span's covered scores as
# text under MAX and reports the last, so a span whose scores are all negative reports the
# most negative (``-5.2`` sorts after ``-2.17``); the candidate reports the numeric maximum.
# The shape: both values numbers, the reference's below the candidate's, the candidate's below
# zero.
CONSERVATION_MAX_TEXT_ORDER = "conservation_max_text_order"
CONSERVATION_FIELDS = ("Conservation", "conservation")
# ``gwas_haplotype_cross_product``: GWAS.pm:338-340 hands every location of a haplotype or
# interaction row (``CHR_ID`` ``21;21;...``, ``SNPS`` ``rs1; rs2; ...``) to every SNP of the
# row, so a site carries every SNP's risk allele and study; the candidate pairs SNP i with
# location i. The reference annotates a superset of sites, so the shape is a ``GWAS_`` field
# the reference fills and the candidate leaves empty, or one both fill differently. A field
# the candidate alone fills is the class's where both rows are annotated under different
# studies: the reference reports another SNP's catalog row, which leaves the field empty, and
# the candidate the site's own; on any other row it stays a difference. A field the catalog row
# leaves empty is a defined empty string on an annotated reference row (``GWAS_beta_coef=``)
# and no key on a candidate row the plugin never annotated: the pair is the class's where the
# reference row fills another ``GWAS_`` field and a difference where it fills none. The class
# reads the catalog's rows of two or more SNPs, as the interaction class below says. A row of
# several locations whose SNPs differ from them in count is skipped by the candidate and
# written whole by the reference, so every SNP's allele at every location of it is the class's.
GWAS_HAPLOTYPE_CROSS_PRODUCT = "gwas_haplotype_cross_product"
GWAS_PREFIX = "gwas_"
# ``gwas_interaction_first_snp_allele_whitespace``: GWAS.pm:343-348 take each risk allele of a
# ``SNP x SNP`` row as the split on ``x`` leaves it, so every allele but the last keeps the
# space before the `` x `` and matches no variant; the candidate trims it. The shape: a
# ``GWAS_`` field the candidate alone fills on a row at the location of such an SNP, the row's
# allele that SNP's risk allele and the candidate's ``GWAS_study`` naming the catalog row's
# study; an empty field of the candidate's annotated row (``GWAS_odds_ratio=``) is the class's
# with its filled ones. Neither output says which catalog row annotated a site, so this class and
# ``gwas_haplotype_cross_product`` read the catalog's rows of two or more SNPs
# (``GwasCatalogPairs``, ``--gwas-catalog``) and accept nothing without them: the cross-product
# class a row at a location of such a row whose allele is another of its SNPs' risk alleles
# and whose reference ``GWAS_study`` names its study.
GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE = "gwas_interaction_first_snp_allele_whitespace"
GWAS_STUDY_FIELDS = ("GWAS_study",)
GWAS_CATALOG_COLUMNS = ("CHR_ID", "CHR_POS", "STRONGEST SNP-RISK ALLELE", "STUDY ACCESSION")
GWAS_SNPS_COLUMN = "SNPS"
GWAS_LIST_SPLIT_RE = re.compile(r"[;,x]")
GWAS_RS_ID_RE = re.compile(r"rs[0-9]+")
# ``alphamissense_two_digit_transcript_version``: AlphaMissense.pm:259 strips a prediction's
# transcript version with ``s/\.[0-9]//g``, one digit, so under ``transcript_match=1`` a
# two-digit version (``ENST00000359568.10`` read as ``ENST0000035956810``) matches no row and
# the reference writes no ``am_`` field there; the candidate strips the whole version. The
# shape: an ``am_`` field (the JSON ``alphamissense`` object) the candidate alone fills.
ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION = "alphamissense_two_digit_transcript_version"
ALPHAMISSENSE_PREFIXES = ("am_", "alphamissense")
# ``cadd_structural_allele_match``: CADD.pm:222-223 gives a structural allele outside
# ``insertion``, ``deletion`` and ``duplication`` no allele string, and Sequence.pm:1182 and
# :1196 reverse-complement a reverse-strand record's class word, so a structural row of an
# Ensembl-format reverse-strand record matches no CADD-SV row in Ensembl VEP and every row of
# the other classes raises ``Plugin 'CADD' went wrong`` (OutputFactory.pm:2244) in its
# warnings; the candidate matches the row and raises nothing. The shape: a ``CADD_`` field the
# candidate alone fills on a structural allele's row, and in --warnings mode the reference's
# ``Plugin 'CADD' went wrong`` lines.
CADD_STRUCTURAL_ALLELE_MATCH = "cadd_structural_allele_match"
CADD_PREFIX = "cadd_"
# ``splicevault_sibling_allele_rows``: SpliceVault.pm:156 matches a table row against every
# alternate allele of the record (``alt_alleles``), so a row naming one allele of a
# multi-allelic record fills the ``SpliceVault_`` fields of every allele's rows; the candidate
# matches the row's own allele. The shape: a ``SpliceVault_`` field (the JSON ``splicevault``
# object) the reference fills with the value its row of another allele of the record carries
# on the same feature, where the candidate leaves it empty or fills it with its own allele's
# row; a field the candidate alone fills, and any value on a record of one allele, stays a
# difference. ``splicevault_json_later_rows``:
# under ``transcript_match=0`` the plugin rewrites its cached row at :190 and dies at :183 on a
# record's later transcript rows in the JSON format alone (``Plugin 'SpliceVault' went wrong:
# Not an ARRAY reference``), so the reference's ``splicevault`` object is missing there; the
# shape: the JSON ``splicevault`` object the candidate alone fills.
SPLICEVAULT_SIBLING_ALLELE_ROWS = "splicevault_sibling_allele_rows"
SPLICEVAULT_JSON_LATER_ROWS = "splicevault_json_later_rows"
SPLICEVAULT_PREFIXES = ("splicevault_", "splicevault")
SPLICEVAULT_JSON_FIELD = "splicevault"
# ``utrannotator_gained_start_window``: UTRAnnotator.pm:339-345 tests a gained start codon over
# a window of the reference allele's length, so an insertion completing an ``ATG`` whose third
# base lies past the window reports the lost codon alone, and a long deletion's window wraps
# to the UTR's end through negative indices and reports a pre-existing ``ATG`` as gained; the
# candidate sizes the window by the alternate allele. The window opens two bases before the
# variant, so on a substitution at cDNA position 1 or 2 it wraps the same way. The shape
# (``gained_start_window_row``): one engine alone reports the gained start, at a
# ``CapDistanceToStart`` its own window holds and the other's does not, the candidate on a row
# whose allele is longer than the reference allele and the reference on a row whose allele is
# shorter; or the reference reports it at a negative ``CapDistanceToStart``, which the wrapped
# window alone yields. The row's ``5UTR_annotation`` goes with its ``5UTR_consequence``. A
# sequence allele at a two-base Location under an identifier is an insertion's or a two-base
# reference allele's, so such a row stays a difference without the record's own alleles. A
# start codon an insertion only moves on is none the inserted bases make: where the reference
# reports a start codon lost and the candidate's gained start stands the allele's length past
# that codon's index, the row is refused.
UTRANNOTATOR_GAINED_START_WINDOW = "utrannotator_gained_start_window"
UTRANNOTATOR_PREFIX = "5utr_"
UTRANNOTATOR_GAINED_START_TERM = "5_prime_UTR_premature_start_codon_gain_variant"
# The position of a reported start codon inside ``5UTR_annotation``: ``CapDistanceToStart=64``
# in the text formats, a member of the JSON object.
CAP_DISTANCE_RE = re.compile(r'CapDistanceToStart"?\s*[=:]\s*"?(-?\d+)')
CDNA_START_FIELDS = ("cDNA_position", "cdna_start")
CDNA_END_FIELD = "cdna_end"
# ``utrannotator_consequence_precedence``: when several of the plugin's five analyses fire on
# one row, UTRAnnotator.pm:262-267 walks their flags as a hash and overwrites the output with
# every truthy one, so ``5UTR_consequence`` and ``5UTR_annotation`` carry whichever effect the
# process's hash order listed last; the candidate reports one under a fixed precedence. The
# shape: a ``5UTR_consequence`` whose two values are both among the five ``5_prime_UTR_``
# terms (``:426``, ``:562``, ``:704``, ``:838``, ``:1003``) and differ, set aside with the
# row's ``5UTR_annotation``; a ``-`` on either side stays a difference. Two terms on a row do
# not show that two analyses fire there, so the class reads Ensembl VEP's own runs of the
# input under other Perl hash seeds (``UtrannotatorSeedSweep``) and accepts nothing without
# them: both terms are ones those runs write on the row, and each side's ``5UTR_annotation``
# is one they write beside its term. The class follows
# ``utrannotator_gained_start_window`` in ``FIELD_CLASSES``, so a row that class explains is
# counted under it alone.
UTRANNOTATOR_CONSEQUENCE_PRECEDENCE = "utrannotator_consequence_precedence"
UTRANNOTATOR_CONSEQUENCE_FIELDS = ("5UTR_consequence", "5utr_consequence")
UTRANNOTATOR_ANNOTATION_FIELDS = ("5UTR_annotation", "5utr_annotation")
UTRANNOTATOR_TERMS = frozenset(
    (
        "5_prime_UTR_premature_start_codon_gain_variant",
        "5_prime_UTR_uORF_stop_codon_gain_variant",
        "5_prime_UTR_uORF_stop_codon_loss_variant",
        "5_prime_UTR_premature_start_codon_loss_variant",
        "5_prime_UTR_uORF_frameshift_variant",
    )
)
# ``utrannotator_annotation_pair_order``: UTRAnnotator.pm:282-283 walks the reported effect's
# groups and each group's ``key=value`` pairs as hashes, so the groups of ``5UTR_annotation``
# and the pairs inside each come in the process's hash order; the candidate writes each
# group's pairs in byte order of their keys. The shape: a ``5UTR_annotation`` whose two values
# carry the same groups of the same pairs in different orders, the candidate's pairs in key
# order. The class precedes the plugin's other seven in ``FIELD_CLASSES``, so a value differing
# in order alone is counted under it on every row.
UTRANNOTATOR_ANNOTATION_PAIR_ORDER = "utrannotator_annotation_pair_order"
# ``utrannotator_structural_allele_text``: UTRAnnotator.pm:80-82 accept a structural variant,
# whose allele string is the record's ALT text (Parser/VCF.pm:572), and ``run`` (:151-153)
# reads that text as the reference allele beside an empty alternate, so the plugin removes as
# many bases of the 5' UTR as the text has characters (:1158-1164, :1197, :1297) and reports
# the effect of that deletion; the candidate annotates a symbolic deletion whose stated span
# lies inside one exon of the 5' UTR as the sequence deletion of that span, gives a plain
# symbolic insertion the frameshift its stated length decides and leaves every other
# structural allele's two fields empty. Three shapes. Where the candidate fills neither
# field (``structural_allele_text_row``): on a row whose Allele is a structural class word
# the reference's ``5UTR_consequence`` is one of the plugin's five terms; the row's
# ``5UTR_annotation`` goes with its ``5UTR_consequence``. Where the candidate fills them
# under a deletion's class word (``structural_deletion_span_row``): a default or tab row
# whose Location's span is the distance between its two cDNA positions and whose
# Consequence names ``5_prime_UTR_variant`` beside no term but ``feature_truncation`` and
# ``NMD_transcript_variant``, so the span lies inside one exon of the 5' UTR; the
# candidate's ``5UTR_consequence`` is a term a deletion can fire beside a ``5UTR_annotation``.
# Where the candidate fills them under ``insertion`` (``structural_insertion_length_row``):
# a default or tab row whose Consequence names ``5_prime_UTR_variant`` and whose Location's
# span, the stated length, is no multiple of three; the candidate's ``5UTR_consequence`` is
# the frameshift, and each effect of its ``5UTR_annotation`` holds the six keys read from
# the reference UTR, with or without ``alt_type=uORF``, and no ``alt_type_length``. Under
# both, the reference's ``5UTR_consequence`` is one of the plugin's five terms or both of
# its fields are empty. A sequence allele's row, a breakend's own string, a row the
# candidate fills under any other class word or outside those two shapes, a VCF or JSON
# entry the candidate fills, an annotation beside no term and the plugin's ``Existing_``
# counts stay differences. The class follows the pair-order, the gained-start window and
# the precedence classes in ``FIELD_CLASSES``, so a row one of them explains is counted
# under that one.
UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT = "utrannotator_structural_allele_text"
# The Allele of a structural allele's row is the Sequence Ontology term of the record's class,
# a value of ``%SO_TERMS`` (ensembl-variation Utils/Config.pm:1417-1438); a breakend with a
# mate writes its own string there instead. The breakpoint class's term is left out: for
# such a record Parser/VCF.pm:499 puts the REF column before the allele string built from
# the ALT text unless that string opens or ends with a dot, and the plugin's reference
# allele is the part that stands first (UTRAnnotator.pm:151-152), the REF column or the
# first ALT text. The term does not tell the two apart, so its row stays a difference.
STRUCTURAL_CLASS_WORDS = frozenset(
    (
        "insertion",
        "mobile_element_insertion",
        "Alu_insertion",
        "HERV_insertion",
        "LINE1_insertion",
        "SVA_insertion",
        "deletion",
        "mobile_element_deletion",
        "Alu_deletion",
        "HERV_deletion",
        "LINE1_deletion",
        "SVA_deletion",
        "tandem_repeat",
        "tandem_duplication",
        "duplication",
        "copy_number_variation",
        "inversion",
    )
)
# The class words of a deletion's row among them. A deletion of a stated span fires any of
# the plugin's analyses but the gained stop, which is a substitution's alone
# (UTRAnnotator.pm:474), and the frameshift only where the span is no multiple of three
# (:923). The row of a deletion inside one exon of the 5' UTR names ``5_prime_UTR_variant``
# beside no Consequence term but these two.
STRUCTURAL_DELETION_WORDS = frozenset(word for word in STRUCTURAL_CLASS_WORDS if word.endswith("_deletion") or word == "deletion")
UTRANNOTATOR_FRAMESHIFT_TERM = "5_prime_UTR_uORF_frameshift_variant"
UTRANNOTATOR_DELETION_TERMS = UTRANNOTATOR_TERMS - {"5_prime_UTR_uORF_stop_codon_gain_variant"}
STRUCTURAL_DELETION_ROW_TERMS = frozenset({"feature_truncation", "NMD_transcript_variant"})
# A plain insertion's class word; a mobile element's insertion has a word of its own. A
# stated length that is no multiple of three shifts the frame whatever is inserted
# (UTRAnnotator.pm:889, :923), so the candidate writes the frameshift alone on such a row:
# each effect with the six keys that read no inserted base, and ``alt_type`` as ``uORF``
# where a stop codon past the point settles it.
STRUCTURAL_INSERTION_WORD = "insertion"
UTRANNOTATOR_REFERENCE_UTR_KEYS = frozenset(
    ("Evidence", "KozakContext", "KozakStrength", "ref_StartDistanceToCDS", "ref_type", "ref_type_length")
)
UTRANNOTATOR_ALTERNATE_TYPE = ("alt_type", "uORF")
# ``utrannotator_lost_start_after_insertion``: UTRAnnotator.pm:764-769 pass an insertion whose
# index lies at one of a start codon's three bases, and ``:776-778`` read the three bases at
# the codon's own index of the mutated UTR and report the codon lost unless they spell ``ATG``.
# Where the insertion only moves the codon on, the codon stands whole, as many bases further
# from the cap as the allele is long, and the candidate reports no lost start: the insertion
# stands directly before the codon, or after its first or second base with inserted bases that
# end with that base or those two, which is the same mutated sequence. The shape
# (``lost_start_after_insertion_row``): the reference's ``5UTR_consequence`` is the lost-start
# term and its ``5UTR_annotation`` names one ``CapDistanceToStart``; the row's cDNA positions
# are an insertion's two flanking bases, the lower of them the 5' UTR index of the base after
# the insertion; that number is the reported index, or lies one above it on a row whose
# ``Allele``, read on the transcript's strand (``STRAND``), ends with ``A``, or two above it on
# a row whose allele so read ends with ``AT``; the candidate's row carries the plugin's counts
# and its ``5UTR_consequence`` is not the lost-start term. The row's ``5UTR_annotation`` goes
# with its ``5UTR_consequence``. A row that states no strand is read at the codon's own index
# alone. A reference allele's row reports an index below its first position (``:764-769``), so
# the equality is an insertion's alone; a row whose own record states a reference allele that
# is not empty is refused all the same. The class follows the gained-start window, the
# precedence and the pair-order classes in ``FIELD_CLASSES``, so a row the gained-start window
# or the seeded runs explain is counted under that class.
UTRANNOTATOR_LOST_START_AFTER_INSERTION = "utrannotator_lost_start_after_insertion"
UTRANNOTATOR_LOST_START_TERM = "5_prime_UTR_premature_start_codon_loss_variant"
START_CODON = "ATG"
# ``utrannotator_gained_start_remade_in_place``: UTRAnnotator.pm:345 opens the gained-start
# window two bases before an insertion's index and ``:346-353`` report the first ``ATG`` in it,
# whatever the reference holds. An insertion written after the first or second base of a
# start codon whose inserted bases open with the codon's remaining bases, and do not also end
# with the codon's bases before them, leaves an ``ATG`` at the codon's own index: the mutated
# UTR is the reference with the inserted bases directly after the untouched codon, and that
# codon is reported gained. The candidate reports no gained start for it and reads the
# frameshift directly after the codon, so the codon's own ORF is among the frames it lists
# unless a stop codon or the UTR's end follows the codon directly, where ``:921`` reads the
# insertion inside the codon and leaves that ORF out. The shape
# (``gained_start_remade_in_place_row``): the row's cDNA
# positions are an insertion's two flanking bases; the reference's ``5UTR_consequence`` is the
# gained-start term and its ``5UTR_annotation`` names one ``CapDistanceToStart``, one below
# the lower cDNA position on a row whose ``Allele``, read on the transcript's strand
# (``STRAND``), opens with ``TG``, or two below it on a row whose allele so read opens with
# ``G``; the candidate's row carries the plugin's counts and reports no gained start at that
# index. Where Ensembl VEP's hash order keeps another of its terms on the row, the gained
# start is the one its runs under other hash seeds (``UtrannotatorSeedSweep``) write there,
# at one index. The row's ``5UTR_annotation`` goes with its ``5UTR_consequence``. A row that
# states no strand, and a row whose own record states a reference allele that is not empty,
# is refused. The class follows the plugin's other five in ``FIELD_CLASSES``, so a row the
# gained-start window or the seeded runs explain is counted under that class.
UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE = "utrannotator_gained_start_remade_in_place"
# ``utrannotator_frameshift_at_stop_boundary``: UTRAnnotator.pm:911-917 set an upstream ORF's
# check point one index before its stop codon, and ``:921`` admits the ORF to the frameshift
# only where the record's own end index is at or before the check point, reading no base of
# the mutated UTR. An insertion written directly before the first base of a stop codon stands
# inside the ORF that ends there: at a length that is no multiple of three the ORF is shifted
# unless the inserted bases put a stop codon at that index; it then runs on and its length changes.
# ``:921`` leaves that ORF out (the comment at ``:919-920`` gives the boundary to the stop
# codon's analysis), and ``:631`` skips the lost-stop test for an insertion at the stop
# codon's first index, so no term reports the ORF; the candidate lists it. The shape
# (``frameshift_at_stop_boundary_row``): both ``5UTR_consequence`` values are the frameshift
# term; the reference's ``5UTR_annotation`` groups are a proper subset of the candidate's;
# every group the candidate adds is a uORF (``ref_type=uORF``) whose ``alt_type_length``
# differs from its ``ref_type_length`` and whose stop codon's first base, its
# ``ref_StartDistanceToCDS`` less its ``ref_type_length`` plus three bases before the coding
# start, is the base directly after the insertion as the record writes it
# (``utr_base_after_insertion``: the row's HGVSc coordinates moved back by the size of
# ``HGVS_OFFSET``);
# and no group both list is such an ORF. A row whose HGVSc states no insertion point in the
# 5' UTR is refused. The class sets aside ``5UTR_annotation`` alone and follows the plugin's
# other six in ``FIELD_CLASSES``, so a row the seeded runs or a start codon left in place
# explain is counted under that class.
UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY = "utrannotator_frameshift_at_stop_boundary"
UTRANNOTATOR_FRAMESHIFT_TERM = "5_prime_UTR_uORF_frameshift_variant"
# An insertion or a duplication between two 5' UTR positions as an HGVSc body writes it.
UTR_HGVSC_INSERTION_RE = re.compile(r"c\.(-\d+)_(-\d+)ins[A-Za-z]+")
UTR_HGVSC_DUPLICATION_RE = re.compile(r"c\.(-\d+)(?:_(-\d+))?dup")
# ``utrannotator_lost_stop_written_before_stop_codon``: UTRAnnotator.pm:628-633 read the three
# bases at a stop codon's index of the mutated UTR only where the record's own index lies one
# or two bases past the stop codon's first (``:631`` skips an insertion at or before it). An
# insertion written directly before a start codon whose second base opens the stop codon of
# an upstream ORF, with inserted bases that open with the start codon's first two bases, is
# the mutated sequence of the same bases written after the start codon's second base, inside
# the stop codon: the stop codon's index reads no stop codon and the stop codon is lost.
# Ensembl VEP reports the lost stop for the record written inside the stop codon and, for the
# record written before the start codon, the ORF's frameshift (``:921``); the candidate
# reports the lost stop on every way of writing the sequence. The shape
# (``lost_stop_written_before_stop_codon_row``): the reference's ``5UTR_consequence`` is the
# frameshift term and the candidate's the lost-stop term; the row's cDNA positions are an
# insertion's two flanking bases and its record, where it states its alleles, an empty
# reference allele; each group of the candidate's ``5UTR_annotation`` carries the
# ``KozakContext`` and ``KozakStrength`` of exactly one uORF (``ref_type=uORF``) the
# reference lists, a different one for each group, and the alternative stop codon where that
# uORF's ``alt_type_length`` ends it: ``AltStop=True`` with ``AltStopDistanceToCDS`` equal to
# the ORF's ``ref_StartDistanceToCDS`` plus the allele's length less its ``alt_type_length``
# plus three where its ``alt_type`` is ``uORF``, ``AltStop=False`` with ``NA`` otherwise. The
# row's ``5UTR_annotation`` goes with its ``5UTR_consequence``. The class follows the
# plugin's other seven in ``FIELD_CLASSES``, so a row whose two terms the seeded runs write is
# counted under ``utrannotator_consequence_precedence``.
UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON = "utrannotator_lost_stop_written_before_stop_codon"
UTRANNOTATOR_LOST_STOP_TERM = "5_prime_UTR_uORF_stop_codon_loss_variant"
UTRANNOTATOR_COUNT_FIELDS = ("Existing_uORFs",)
# ``nearest_exon_jb_item_order``: NearestExonJB.pm:267-271, and ``:245-263`` under
# ``intronic=1``, build the value by walking the keys of the hash of the row's exons, so two
# or more exon items come in the process's hash order; the candidate lists them in transcript
# order. The shape: a ``NearestExonJB`` whose two values carry the same items the same number
# of times in different orders.
NEAREST_EXON_JB_ITEM_ORDER = "nearest_exon_jb_item_order"
# ``loeuf_zero_minimum``: LOEUF.pm:204 drops a gene's minimum whose string is ``0`` by Perl
# truthiness; the candidate prints it. The shape: a ``LOEUF`` the candidate alone fills, with
# ``0``.
LOEUF_ZERO_MINIMUM = "loeuf_zero_minimum"
LOEUF_FIELDS = ("LOEUF", "loeuf")
# ``maxentscan_overlapping_windows``: a substitution in an exon of at most five bases lies in
# the acceptor window of the intron before it and the donor window of the intron after it;
# MaxEntScan.pm:270 scores the window of whichever intron Ensembl VEP's interval tree returns
# first, the candidate the first window in transcript order, the upstream intron's acceptor. A
# substitution in both windows of one intron, an intron of at most 25 bases, is scored at the
# donor by both. The shape: on a row whose record is one base in the acceptor window of one
# intron of the row's transcript and the donor window of the next, read from the candidate's
# JSON cache (``SpliceWindows``), a ``MaxEntScan_`` field both engines fill with different
# numbers, and a window text of a ``verbose`` run where the reference writes a donor's nine
# bases and the candidate an acceptor's twenty-three. Without the cache nothing is accepted.
MAXENTSCAN_OVERLAPPING_WINDOWS = "maxentscan_overlapping_windows"
MAXENTSCAN_PREFIX = "maxentscan_"
MAXENTSCAN_WINDOW_TEXT_FIELDS = ("MaxEntScan_ref_seq", "MaxEntScan_alt_seq")
MAXENTSCAN_DONOR_TEXT_RE = re.compile(r"[ACGT]{9}")
MAXENTSCAN_ACCEPTOR_TEXT_RE = re.compile(r"[ACGT]{23}")
# ``alphamissense_sibling_allele_rows``: AlphaMissense.pm:245 matches a table row against every
# alternate allele of the record and :250-253 keep the first row whose amino-acid letters equal
# the output row's, so two alleles of one record giving one amino-acid change both take the
# table row that comes first; the candidate matches each allele's own row. The shape: an
# ``am_`` field (the JSON ``alphamissense`` object) the reference fills with the value its row
# of another allele of the record carries on the same feature beside the same ``Amino_acids``,
# where the candidate leaves it empty or fills it differently.
ALPHAMISSENSE_SIBLING_ALLELE_ROWS = "alphamissense_sibling_allele_rows"
AMINO_ACIDS_FIELDS = ("Amino_acids",)
# ``revel_sibling_allele_rows``: REVEL.pm:165 matches a score row against every alternate
# allele of the record and :176-199 return the first row whose alternate amino acid equals
# the output row's peptide, the widening AlphaMissense's match has. The shape: a ``REVEL``
# the reference fills with the value its row of another allele of the record carries on the
# same feature beside the same ``Amino_acids``, where the candidate leaves it empty or fills
# it differently.
REVEL_SIBLING_ALLELE_ROWS = "revel_sibling_allele_rows"
REVEL_FIELDS = ("REVEL",)
# ``revel_other_reference_amino_acid_rows``: a REVEL table row scores one substitution, from
# its ``aaref`` to its ``aaalt``, on the transcripts it lists. REVEL.pm:210 and :229 read the
# ``aaref`` cell and :193 and :197 test the ``aaalt`` alone, so the reference writes the score
# of a row whose reference amino acid is not the output row's; the candidate tests both
# cells. The shape: a ``REVEL`` the reference fills with the score of the table row its walk
# reaches (REVEL.pm:170-199), a row whose ``aaref`` is not the first residue of the row's
# ``Amino_acids``, where the candidate writes what its own walk gives: nothing, or the score
# of the first row of the row's own alternate that names both residues. The class reads the
# table the run's plugin read (``RevelTable``, ``--revel-table``) and accepts nothing
# without it.
REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS = "revel_other_reference_amino_acid_rows"
# ``gnomadmt_info_value_second_equals``: gnomADMt.pm:327-330 reads each INFO item of a file
# record with ``my ($key, $value) = split /=/, $field``, which Perl limits to three parts, so
# a value holding ``=`` is cut there; the candidate splits the item at its first ``=`` alone.
# The shape: a ``gnomAD_`` field both fill, the reference's value the candidate's up to its
# first ``=``, read with the percent-encoding of ``=`` undone.
GNOMADMT_INFO_VALUE_SECOND_EQUALS = "gnomadmt_info_value_second_equals"
GNOMAD_PLUGIN_PREFIX = "gnomad_"
# ``loftee_anc_allele_multiallelic_record``: LoF.pm:531 returns before the ancestral comparison
# when the record's allele string holds ``-`` or is not three characters long, so no allele of
# a multi-allelic record is filtered ``ANC_ALLELE``; the candidate judges each one-base
# substitution by its own base. The shape: on a one-base substitution's stop-gained or
# splice-site row of a record that has a row of another allele on the same feature, a
# ``LoF_filter`` the candidate
# writes as the reference's plus ``ANC_ALLELE``, set aside with the row's ``LoF``, ``LC``
# beside the item and ``HC`` on a side with no item.
LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD = "loftee_anc_allele_multiallelic_record"
LOF_FIELDS = ("LoF",)
LOF_FILTER_FIELDS = ("LoF_filter",)
ANC_ALLELE_FILTER = "ANC_ALLELE"
LOF_SUBSTITUTION_TERMS = frozenset(("stop_gained", "splice_acceptor_variant", "splice_donor_variant"))
# ``loftee_first_alt_donor_filter``: LoF.pm:363-367 test the reference and the first alternate
# allele of the record for every allele's row, so a later allele of a multi-allelic donor
# record is filtered ``GC_TO_GT_DONOR``, or spared, by its sibling's base; the candidate tests
# each allele's own base. The shape: on a ``splice_donor_variant`` row of a record that has a
# row of another allele on the same feature, and that is not the row of the record's first
# alternate, a ``LoF_filter`` whose two values differ by ``GC_TO_GT_DONOR`` alone, either
# side carrying it, set aside with the row's ``LoF`` as the ancestral class's is. The first
# alternate is read where an entry states the record's alleles in order (``first_alternate``);
# where none does the class accepts nothing.
LOFTEE_FIRST_ALT_DONOR_FILTER = "loftee_first_alt_donor_filter"
GC_TO_GT_DONOR_FILTER = "GC_TO_GT_DONOR"
LOF_DONOR_TERM = "splice_donor_variant"
# ``loftee_gerp_zoom_level_summary``: gerp_dist.pl:70-79 total each exon interval's GERP score
# through one summary bin of the bigWig, which the library serves from a zoom level once the
# interval is twice that level's span, 80 bases on the file LoFTEE's GRCh38 bundle ships, so
# ``GERP_DIST`` approximates the per-base sum the candidate writes, and LoF.pm:254 pushes
# ``END_TRUNC`` by that total. The shape: on a stop-gained or frameshift row whose ``BP_DIST``
# is at least 79 (an interval of N bases walks N - 1, so 79 is the 80-base interval), a
# ``LoF_info`` whose two values differ in the ``GERP_DIST`` item alone, both
# numbers; and, where the row's ``DIST_FROM_LAST_EXON`` is at most 50 and the two totals
# straddle the cutoff -58, a ``LoF_filter`` differing by ``END_TRUNC`` alone on the side of
# the lower total, set aside with the row's ``LoF`` as the other LoFTEE classes' is. The
# cutoff is the plugin's default: a run that sets ``gerp_end_trunc_cutoff`` is not read.
LOFTEE_GERP_ZOOM_LEVEL_SUMMARY = "loftee_gerp_zoom_level_summary"
LOF_INFO_FIELDS = ("LoF_info",)
LOF_TRUNCATING_TERMS = frozenset(("stop_gained", "frameshift_variant"))
END_TRUNC_FILTER = "END_TRUNC"
GERP_END_TRUNC_CUTOFF = -58.0
GERP_EXACT_BELOW_BASES = 80
LAST_EXON_REACH = 50
# ``nearest_exon_jb_stretch_priming``: NearestExonJB.pm:184 reads ``_overlapped_exons``, whose
# 12-base stretch (BaseTranscriptVariation.pm:858-865) applies only once a record inside the
# transcript has set its frameshift-intron flag (:1032-1034), so a record within 12 bases
# outside the first or last exon gets a value or none by the records before it; the candidate
# decides the stretch from the transcript's introns. The shape: a ``NearestExonJB`` the
# candidate alone fills on a row whose ``DISTANCE`` is at most 12.
NEAREST_EXON_JB_STRETCH_PRIMING = "nearest_exon_jb_stretch_priming"
NEAREST_EXON_JB_FIELDS = ("NearestExonJB",)
DISTANCE_FIELDS = ("DISTANCE",)
FRAMESHIFT_INTRON_STRETCH = 12
# ``nearest_exon_jb_shifted_span``: under ``--shift_3prime 1`` the consequence tests fill the
# ``_overlapped_exons`` memo at the shifted span (BaseVariationFeatureOverlapAllele.pm:454-463,
# :514) and NearestExonJB.pm:184 reads that memo while :225-226 measure from the record's
# span; the candidate reads the exons at the record's span. The shape: a ``NearestExonJB`` on
# a transcript row with a positive ``SHIFT_LENGTH`` that one side fills and the other leaves
# empty, or whose two values name different exons.
NEAREST_EXON_JB_SHIFTED_SPAN = "nearest_exon_jb_shifted_span"
# One item of a ``NearestExonJB`` value, the exon id, distance, side (``start_end`` on a tie)
# and length: ``|`` joins the parts and ``,`` the items in the default, tab and JSON formats,
# ``+`` and ``&`` in the VCF format.
NEAREST_EXON_ITEM = r"([^|+&,]*)[|+](\d+)[|+](start|end|start_end)[|+](\d+)"
NEAREST_EXON_ITEM_RE = re.compile(NEAREST_EXON_ITEM)
NEAREST_EXON_LIST_RE = re.compile(rf"{NEAREST_EXON_ITEM}(?:[,&]{NEAREST_EXON_ITEM})*")
# ``nearest_exon_jb_id_less_exon_fold``: NearestExonJB.pm:220-237 file each exon's distance,
# side and length under its stable id, so the exons of a row that carry none share the empty
# key and :267-271 (:244-263 under ``intronic=1``) print one item from their parts, the
# nearest exon's distance and side with the length filed last; the candidate prints each
# exon's own item. The shape: a ``NearestExonJB`` both fill with items of an empty exon id,
# the reference's one item whose distance and side are those of one of the candidate's items,
# on a row the candidate gives two or more items, or one whose Consequence carries both
# splice-site terms and whose cDNA position names two positions: a variant that begins in one
# exon and ends in another. A deletion over one whole exon carries both terms and no such
# position, and overlaps that exon alone.
NEAREST_EXON_JB_ID_LESS_EXON_FOLD = "nearest_exon_jb_id_less_exon_fold"
# ``downstream_reverse_strand_rotation``: TranscriptVariationAllele.pm:1330-1347 rotate a
# 3'-shifted allele ``seq_length - shift_length`` places on a reverse-strand transcript, none
# once the shift exceeds the allele, and Downstream.pm:110 splices the bases as they stand;
# the candidate rotates by the shift modulo the allele's length. The shape: a
# ``DownstreamProtein`` or ``ProteinLengthChange`` both fill differently on a ``frameshift_variant`` row
# of a ``STRAND`` -1 transcript whose ``HGVS_OFFSET`` exceeds the number of bases its HGVSc
# inserts, duplicates or deletes and is not a multiple of it.
DOWNSTREAM_REVERSE_STRAND_ROTATION = "downstream_reverse_strand_rotation"
DOWNSTREAM_FIELDS = ("DownstreamProtein", "ProteinLengthChange")
# ``spliceai_chopped_multiallelic``: for a record whose allele string holds ``-``,
# SpliceAI.pm:259-266 rebuild every allele with the anchor base before matching, so a
# substitution beside a deletion is compared as a two-base pair, which equals no single-base
# row (:272); the candidate compares each allele at its own position. The shape: a ``SpliceAI_``
# field (the JSON ``spliceai`` object) the candidate alone fills on a substitution allele's
# row of a record whose rows carry an allele written ``-``.
SPLICEAI_CHOPPED_MULTIALLELIC = "spliceai_chopped_multiallelic"
SPLICEAI_PREFIXES = ("SpliceAI_", "spliceai")
# ``spliceai_multi_entry_cell``: SpliceAI.pm:326-332 remove ``SpliceAI=`` from a score record's
# INFO column and split the rest on ``|`` into the allele, the first symbol and everything
# after it, and :305 file the record under that symbol alone. A record whose key holds several
# comma-joined entries so gives the rows of the first entry's symbol the whole remainder as
# their value (the nine fields of the other formats cut from it, :276-287) and the rows of a
# later entry's symbol nothing (:314-316); the candidate reads each entry under its own
# symbol. Where a record line stands more than once at an allele, as in the SpliceAI golden
# corpora, :305 keep one line a first symbol, the last, so the reference writes a row the
# whole value of the last line its gene leads, or nothing where it leads none, and the
# candidate the gene's entry; where the lines give the gene entries that differ in text, the
# candidate's is that of the last line that names the gene. The shape is read from the score
# files themselves (``SpliceAiTable``,
# ``--spliceai-table``) and nothing is accepted without them. A row's gene is its ``SYMBOL``
# where the row writes one, else its transcript's symbol in the candidate JSON cache
# (``TranscriptSymbols``, ``--json-cache``); a row whose gene neither names is not accepted.
SPLICEAI_MULTI_ENTRY_CELL = "spliceai_multi_entry_cell"
SPLICEAI_KEY = "SpliceAI="
SPLICEAI_VALUE_FIELD = "SpliceAI_pred"
SPLICEAI_CUTOFF_FIELD = "SpliceAI_cutoff"
SPLICEAI_MEMBERS = ("SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL")
SPLICEAI_CUTOFF_FLAGS = ("PASS", "FAIL")
# The grounds on which ``spliceai_multi_entry_cell`` refuses a row it is asked about under the
# score files, in the order ``spliceai_multi_entry_cell_ground`` meets them: a key and what
# the predicate read.
SPLICEAI_GROUNDS = {
    "no_candidate_value": "the candidate's row writes no SpliceAI value",
    "no_span": "the row's Location names no span",
    "other_info": "a record line at the row's variant carries another INFO part beside the key or an entry of another length",
    "several_references": "the record lines at the row's variant spell more than one REF",
    "no_record": "no record line of the score files stands at the row's variant",
    "no_gene": "neither a SYMBOL nor the candidate cache names the row's gene",
    "gene_not_named": "no entry of a record line at the row's variant names the row's gene",
    "lone_entry": "every record line that holds the gene's entry holds it alone",
    "candidate_value": "the candidate's value is not the gene's entry of the last record line that names it",
    "reference_not_silent": "the row's gene leads no record line and the reference's row writes a SpliceAI value or a cutoff flag, or the field is not the candidate's alone",
    "reference_value": "the row's gene leads a record line and the reference's value is not the whole value of the last line it leads",
    "reference_flag": "the row's gene leads a record line and the two cutoff flags differ",
}
# ``plugin_contig_mitochondrial_spelling``: BaseVepTabixPlugin.pm:693-708 map a record's
# chromosome name to a plugin's file by adding or removing ``chr`` alone, so a record spelt
# ``MT``, ``chrM`` or ``chrMT`` finds nothing in a file naming the contig otherwise; the
# candidate queries the file under all four spellings where they name one sequence: on any
# assembly but GRCh37, and there for AlphaMissense alone, whose table counts on the assembly's
# sequence. The shape: a field of a tabix-backed plugin other than gnomADMt that the candidate
# alone fills on a row of a record so spelt and, in a run whose manifest names GRCh37
# (``GRCH37_NAMES``, the names ``Assembly::parse`` reads as that assembly), an AlphaMissense
# field alone; a run that names no assembly is read as one on another assembly than GRCh37.
# The fields are those the plugins name themselves; a column a dbNSFP run selects and a
# TranscriptAnnotator table's fields carry names of the user's and are not read.
PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING = "plugin_contig_mitochondrial_spelling"
MITOCHONDRIAL_SPELLINGS = ("MT", "chrM", "chrMT")
GRCH37_NAMES = ("grch37", "hg19", "37")
TABIX_PLUGIN_FIELDS = (
    *ALPHAMISSENSE_PREFIXES, "AVI_", "BayesDel", CADD_PREFIX, "ClinPred", "ada_score", "rf_score", *LOEUF_FIELDS,
    "OpenTargets_", "opentargets", "PrimateAI", *REVEL_FIELDS, *SPLICEAI_PREFIXES, *SPLICEVAULT_PREFIXES,
    GNOMAD_PLUGIN_PREFIX,
)
# The ``gnomAD_`` fields are two plugins': gnomADc's end in ``_cov``, the others are gnomADMt's.
GNOMADC_FIELD_SUFFIX = "_cov"
# ``dosage_sensitivity_regulatory_row_warnings``: DosageSensitivity.pm's ``run`` reads
# ``$tva->transcript`` on a regulatory or motif row under --regulatory and dies inside
# ``run_plugins``' eval, one ``Plugin 'DosageSensitivity' went wrong`` warning per row
# (OutputFactory.pm:2244); the candidate writes nothing on those rows and no warning. Applied
# in --warnings mode alone: the reference's lines are set aside.
DOSAGE_SENSITIVITY_REGULATORY_ROW_WARNINGS = "dosage_sensitivity_regulatory_row_warnings"
# The reference warnings each warnings-mode class sets aside, by the plugin the line names.
PLUGIN_WARNING_CLASSES = {
    CADD_STRUCTURAL_ALLELE_MATCH: "Plugin 'CADD' went wrong",
    DOSAGE_SENSITIVITY_REGULATORY_ROW_WARNINGS: "Plugin 'DosageSensitivity' went wrong",
}
# ``fork_warning_pid_prefix``: under --fork, Ensembl VEP's parent relays a child's captured
# standard error as one warning headed by the child's process id (Runner.pm:565), and
# ``warning_msg`` (BaseVEP.pm:806) puts ``WARNING: `` before it, so the block's first line reads
# ``WARNING: <pid> : WARNING: <text>`` where the serial run and the candidate write
# ``WARNING: <text>``. Applied in --warnings mode alone: the prefix is removed from the
# reference's lines before they are parsed, each line changed counted under the class.
FORK_WARNING_PID_PREFIX = "fork_warning_pid_prefix"
FORKED_CHILD_RELAY_RE = re.compile(r"^WARNING:\s*\d+ : ")
# The warnings-mode classes that rewrite the reference's lines rather than set them aside.
WARNING_TEXT_CLASSES = (FORK_WARNING_PID_PREFIX,)
# ``refseq_edited_model_hgvsc_offset``: on a transcript model whose alignment to the genome
# inserted or deleted bases (``BAM_EDIT=OK``), Ensembl VEP 116.2 numbers the HGVSc from the
# unedited exons: the cache dumper drops the model's edit flag and mapper (Dumper/Core.pm:164,
# :169), so ``_get_cDNA_position`` rebuilds an unfolded one, and an intronic allele counts no
# edit (TranscriptVariationAllele.pm:1608-1610); the candidate numbers every position against
# the edited model. The shape: an HGVSc on a ``BAM_EDIT`` ``OK`` row whose two values name
# the same reference and differ in their coordinate tokens alone, each pair of coordinates
# covering as many positions in one value as in the other. The numbering moves a notation
# and never stretches it: two values that cover different numbers of positions lie on bases
# the model holds and the genome lacks, where the engines describe different changes, or
# across bases the model lacks, which the two numberings do not count alike. Neither is the
# class's row: the pair is refused and stays a difference, each value printed as its engine
# wrote it. The count is read where a pair states it, from two plain positions or two
# offsets off one exon boundary; a pair that mixes the two, or runs over the stop codon, is
# not measured.
REFSEQ_EDITED_MODEL_HGVSC_OFFSET = "refseq_edited_model_hgvsc_offset"
HGVSC_FIELDS = ("HGVSc", "hgvsc")
BAM_EDIT_FIELDS = ("BAM_EDIT", "bam_edit")
CDNA_COORDINATE_RE = re.compile(r"\*?-?\d+(?:[+-]\d+)?")
CDNA_COORDINATE_PAIR_RE = re.compile(r"(\*?)(-?\d+)([+-]\d+)?_(\*?)(-?\d+)([+-]\d+)?")
# ``refseq_edited_model_hgvs_offset``: on an edited ``NM_`` or ``XM_`` model, ``_return_3prime``
# (TranscriptVariationAllele.pm:164-187) copies the first shift hash attached to the variant,
# whatever transcript or strand made it, because its length test at :182-183 fails for every
# insertion and deletion; a transcript on the other strand inherits that transcript's shift, a
# zero that OutputFactory.pm:1752-1755 does not print, while the 3'-first duplication test of
# Utils/Sequence.pm:578 still names the shifted notation; the candidate shifts each row on its
# own transcript. The shape: an HGVS_OFFSET on a coding ``BAM_EDIT`` ``OK`` row whose HGVSc
# names a duplication, the reference leaving the offset empty and the candidate filling it with
# a non-zero integer, the distance of a 3' shift, negative on the reverse strand. The branch is
# entered only for an allele with CDS coordinates (TranscriptVariationAllele.pm:147-148), and
# only a duplication keeps its notation whatever the shift: a deletion or a plain insertion
# that inherits the foreign shift changes its HGVSc too, another shape.
REFSEQ_EDITED_MODEL_HGVS_OFFSET = "refseq_edited_model_hgvs_offset"
HGVS_OFFSET_FIELDS = ("HGVS_OFFSET", "hgvs_offset")
CDS_POSITION_FIELDS = ("CDS_position", "cds_start", "cds_end")
NON_ZERO_INTEGER_RE = re.compile(r"-?[1-9]\d*")
# ``refseq_edited_model_shifted_insertion_bases``: on a model carrying an ``_rna_edit``
# (``BAM_EDIT=OK``), ``hgvs_transcript`` respells the notation's alleles from ``feature_seq``
# (TranscriptVariationAllele.pm:1510-1513) after ``hgvs_variant_notation`` described an insertion
# at its 3'-shifted coordinates with the rotated allele (:1444, :1474-1483). ``feature_seq`` holds
# the rotated allele only on a coding row under --shift_3prime 1, where ``codon`` (:813-820) has
# run ``shift_feature_seqs`` (:1338-1342) before the notation is built; ``hgvs_protein`` (:1778)
# rotates it too, after ``hgvs_transcript`` has read it (OutputFactory.pm:1751, :1760). On every
# other row the notation spells the inserted bases as written at the shifted position, a sequence
# the variant does not produce. The candidate spells the rotated bases on every row. The shape:
# an HGVSc on a ``BAM_EDIT`` ``OK`` insertion row whose HGVS_OFFSET is one non-zero integer on
# both sides, the two notations insertions of one transcript, numbering and coordinate pair, the
# reference's bases the row's Allele read on the transcript's strand (reverse-complemented where
# STRAND is -1; a row without STRAND admits either reading) and the candidate's those bases turned
# towards the transcript's 3' end by the offset modulo their length; the same bases, another
# length, a duplication and a row without the offset stay differences.
REFSEQ_EDITED_MODEL_SHIFTED_INSERTION_BASES = "refseq_edited_model_shifted_insertion_bases"
# ``refseq_edited_model_seam_insertion_hgvsc``: an insertion standing on the bases an edited
# model holds and the genome lacks. On a RefSeq
# model whose alignment to the genome inserted bases (an ``_rna_edit``; the row carries
# ``BAM_EDIT`` ``OK`` or ``FAILED``) the cached mapper numbers the two genomic bases around
# those bases more than one apart, so an insertion between the two has the model's bases between
# them as its span: a ``CDS_position`` of ``75-77`` for three bases and ``1059`` for one, where
# an ordinary insertion prints two neighbouring positions. ``hgvs_transcript`` slides the
# insertion's two transcript positions by the misalignment offset
# (TranscriptVariationAllele.pm:1517, :1524-1525; ``get_misalignment_offset``, :1594-1624) and
# writes ``ins`` of the allele between them, at plain coordinates (``c.77_78``) or with intron
# offsets (``c.2034+2_2034+3``), naming neither the model's bases nor their span. The candidate
# names the replacement of those bases by the allele, typed as ``hgvs_variant_notation`` types
# an allele replacing reference bases (Utils/Sequence.pm:493), clipped as ``_clip_alleles``
# clips one (TranscriptVariationAllele.pm:1492, :2188) and placed by the HGVS 3' rule on the
# model's own sequence, whose bases beyond the row's codons the row does not carry, so the
# notation is judged by what it says on those codons. The shape: an HGVSc on an insertion's row
# carrying ``BAM_EDIT`` ``OK`` or ``FAILED`` on both sides, the reference's notation an ``ins``
# of the row's Allele read on the transcript's strand (reverse-complemented where STRAND is -1;
# a row without STRAND admits either reading) at any coordinate pair of the ``c.`` numbering,
# the candidate's row with a CDS position of one position or a span inside the codons its
# Protein_position and Codons give, the alternate codons being the reference codons with the
# span's bases replaced by the allele, and the candidate's notation none where the allele is
# those bases and otherwise one that, applied to the reference codons, gives the alternate
# codons: under the reference's sequence name and at plain coding positions, a substitution of
# the base it names, a ``del``, a ``delins``, a ``dup``, an ``inv`` or a repeat count of
# positions inside the codons, or an ``ins`` between two neighbouring positions, of which one
# may lie one base outside them. An ordinary insertion's alternate codons are longer than its
# reference codons by the whole allele, which tells it from a two-base span printing the same
# two positions. A notation that gives another sequence or names a position outside the codons,
# an ``ins`` between positions that are not neighbours, any notation where the allele is the
# model's bases, a row without ``BAM_EDIT`` (a mitochondrial model's), a row without a CDS
# position or Codons (a span in an untranslated region or on a non-coding model, whose row
# does not carry the model's bases), a reference notation that is not an ``ins`` of the allele
# and another reference sequence stay differences
# (``edited_model_seam_insertion_replacement``).
REFSEQ_EDITED_MODEL_SEAM_INSERTION_HGVSC = "refseq_edited_model_seam_insertion_hgvsc"
# ``refseq_edited_model_unminimised_allele_hgvsc``: the HGVSc of an allele that keeps a base it
# shares with its reference on a model whose ``_rna_edit`` the alignment applied; the comment
# above ``GIVEN_REF_FIELDS`` states the shape ``edited_model_clipped_alleles`` reads.
REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC = "refseq_edited_model_unminimised_allele_hgvsc"
EDITED_INSERTION_NOTATION_RE = re.compile(r"([cn])\.(\*?-?\d+(?:[+-]\d+)?_\*?-?\d+(?:[+-]\d+)?)ins([ACGT]+)")
# ``hgvs_shifted_insertion_clamped_end``: under --shift_3prime 1 with --hgvs, an insertion
# anchored on the base before a transcript's first base and shifted k bases into exon 1 maps to
# the transcript slice coordinates (1, 0) in ``_var2transcript_slice_coords``
# (TranscriptVariationAllele.pm:2689-2736, either strand at :2704-2714), which :2729 keeps,
# since only the end is below 1, and :2733 clamps the end to 1; ``hgvs_transcript`` adds k to
# both ends (:1474-1483), so ``hgvs_variant_notation`` reads the one base at n.k+1 against the
# inserted bases and the row prints a substitution or a one-base delins there (``n.15A>T``,
# ``n.15delinsTT``), or, where the inserted sequence ends with that base, ``_clip_alleles``
# (:2231-2242, :2271-2278) strips it and prints an insertion of the other bases at (n.k, n.k+1).
# The candidate writes the notation of the shifted span, as Ensembl VEP does for the allele
# written at the shifted position: a duplication ending at n.k or an insertion between n.k and
# n.k+1. The shape: an HGVSc differing on an insertion's row whose SHIFT_LENGTH k is positive,
# whose cDNA position starts at k (the shift carried the insertion from before the first base to
# after base k) and whose HGVS_OFFSET is set, the reference's notation at one position P
# (``P X>Y``, ``P delins``, or ``P-1_P ins``) and the candidate's a ``dup`` ending at P-1 or an
# ``ins`` spanning P-1 and P, both on the same transcript and numbering, the candidate's
# duplicated or inserted bases as many as the allele's, the reference's as many less the one it
# stripped, and a candidate ``ins`` carrying the reference's bases; an empty HGVSc on either side
# stays a difference.
HGVS_SHIFTED_INSERTION_CLAMPED_END = "hgvs_shifted_insertion_clamped_end"
SHIFT_LENGTH_FIELDS = ("SHIFT_LENGTH", "shift_length")
CLAMPED_NOTATION_RE = re.compile(r"([cn])\.(-?[1-9]\d*)(?:[ACGTN]>([ACGT])|delins([ACGT]+))")
INSERTION_NOTATION_RE = re.compile(r"([cn])\.(-?[1-9]\d*)_(-?[1-9]\d*)ins([ACGT]+)")
SHIFTED_DUPLICATION_RE = re.compile(r"([cn])\.(?:(-?[1-9]\d*)_)?(-?[1-9]\d*)dup")
# ``hgvs_shifted_indel_written_span_gate``: ``_var2transcript_slice_coords``
# (TranscriptVariationAllele.pm:2708-2713) maps the span as written onto the transcript slice and
# refuses one wholly outside it (:2729), and ``hgvs_transcript`` refuses a clamped end the shift's
# offset carries past the last base (:1472) before adding the offset (:1477-1478), so an insertion
# or deletion the shift carries into a transcript from wholly before it, or over its last base,
# has no HGVSc and no HGVS_OFFSET on a row annotated inside the transcript. The candidate
# describes the shifted span, as Ensembl VEP does for the allele written there. The shape: an
# empty reference HGVSc against a candidate insertion or deletion notation of the row's
# transcript, an empty reference HGVS_OFFSET against a candidate non-zero integer, and an empty
# reference HGVSIntronOffset plugin field against a candidate integer, on an indel's row annotated
# inside the transcript (no upstream or downstream term; one whose allele or Location says so, or,
# for an allele of two or more bases at a Location of as many, which the text formats print alike
# for an insertion and a substitution, one whose candidate notation is a ``dup`` or ``ins``) whose
# SHIFT_LENGTH, where the run prints one, is positive; the candidate row carries the notation and
# the offset together, so each pair is read beside the row's other fields. The plugin's two
# fields are the intron offsets of the notation's coordinates, so beside a ``del`` notation,
# whose printed coordinates are the ones the plugin reads, they are held to them.
HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE = "hgvs_shifted_indel_written_span_gate"
HGVS_INTRON_START_OFFSET_FIELDS = ("HGVS_IntronStartOffset",)
HGVS_INTRON_END_OFFSET_FIELDS = ("HGVS_IntronEndOffset",)
HGVS_INTRON_OFFSET_FIELDS = (*HGVS_INTRON_START_OFFSET_FIELDS, *HGVS_INTRON_END_OFFSET_FIELDS)
OUTSIDE_TRANSCRIPT_TERMS = frozenset(("upstream_gene_variant", "downstream_gene_variant"))
SHIFTED_INDEL_GATE_FIELDS = (*HGVSC_FIELDS, *HGVS_OFFSET_FIELDS, *HGVS_INTRON_OFFSET_FIELDS)
# ``hgvs_intron_offset_abutting_insertion``: HGVSIntronOffset.pm:62 asks for the HGVS notation
# on every Transcript row, so on an insertion directly beside a transcript's first or last
# base the reference fills the plugin's two fields from a notation the row never prints, and
# the candidate writes neither where the HGVS 3' rule places no notation on the transcript.
# The shape and the three pairs it takes are the predicate's
# (``intron_offset_abutting_insertion``).
HGVS_INTRON_OFFSET_ABUTTING_INSERTION = "hgvs_intron_offset_abutting_insertion"
INTEGER_RE = re.compile(r"-?\d+")
# ``hgvsp_shifted_indel_written_span_gate``: ``hgvs_protein`` shifts the allele 3' and reads its
# translation coordinates at the shift (TranscriptVariationAllele.pm:1711, 1732-1733), and tests
# the ``coding`` pre-consequence predicate of the allele as written (:1731), computed and cached
# with the row's terms (BaseVariationFeatureOverlapAllele.pm:340, 394, 460-462), where an allele
# written in an intron reads ``non_coding`` (:535-537). So an insertion or deletion written in an
# intron whose shift carries it into the coding sequence has the shifted HGVSc and HGVS_OFFSET
# and no HGVSp, except on a transcript with a frameshift intron, whose exons the predicate reads
# 12 bases wider (BaseTranscriptVariation.pm:861-866). The candidate writes the protein notation
# of the shifted allele, as Ensembl VEP does for the allele written there. The shape: an empty
# reference HGVSp against a candidate protein notation, on an indel's row annotated inside the
# transcript with no cDNA, CDS or protein position (the allele as written touches no exon),
# whose HGVSc and non-zero HGVS_OFFSET both engines write alike, the HGVSc a ``c.`` insertion or
# deletion notation of the row's transcript that places the shifted allele in the coding
# sequence: a ``del`` whose every coordinate is a coding base; a ``dup`` ending on a coding
# base, its first coordinate free to lie back in the intron; an ``ins`` between two coding
# bases. An insertion the shift leaves between the intron's last base and the exon's first
# base is a change at the junction, written as a ``dup`` that ends on the intron's last base
# or as an ``ins`` between that base and the coding base after it: the changed sequence
# reads as a longer intron before the unchanged exon as well as a longer exon, the candidate
# leaves the HGVSp empty there as Ensembl VEP does, and a protein notation on one side alone
# is a difference, whatever bases are inserted and however long the shift. So is an
# insertion before the stop codon's last base, whose intron base is printed ``*-1``,
# without that base's coordinate (TranscriptVariationAllele.pm:2854-2859). The
# intron must lie inside the coding region, since the candidate writes no HGVSp for an allele
# shifted out of an intron of the 5' UTR onto the start codon: the first base of the exon the
# allele is shifted into is coding base 2 or later. A notation with a coordinate in the
# intron names that base (the coordinate counts back from it, or forward from the coding base
# before the intron, which it follows); one of coding bases alone bounds it from below, since
# the allele as written lies before it: the coding coordinate the shifted allele starts at,
# less the bases of the shift, plus 1. A row whose shift is too long for that bound stays a
# difference, as do a row whose SHIFT_LENGTH is positive (annotated at the shifted span), a
# reference notation against an empty candidate, two differing notations, and a candidate
# notation whose first changed residue is the stop codon (``p.Ter``) or that names the stop
# codon as a flank of an insertion (``_Ter85ins``): bases inserted inside the stop codon that
# leave a stop codon after the new residues are written like bases inserted before it. A row
# either engine flags ``cds_start_NF`` or ``cds_end_NF`` (OutputFactory.pm:1478) is no row of the
# class: on a transcript flagged incomplete, and on one whose coding sequence does not end in a
# stop codon, the candidate leaves the HGVSp empty as Ensembl VEP does, so a protein notation on
# one side alone is a difference there. A coding sequence without its closing stop codon shows
# in no field of a row that carries no flag.
HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_GATE = "hgvsp_shifted_indel_written_span_gate"
PROTEIN_POSITION_FIELDS = ("Protein_position", "protein_start", "protein_end")
PROTEIN_NOTATION_RE = re.compile(r"[^:\s]+:p\.\S+")
STOP_CODON_NOTATION_RE = re.compile(r"[^:\s]+:p\.\(?(?:Ter|\*)\d")
STOP_CODON_FLANK_INSERTION_RE = re.compile(r"[^:\s]+:p\.[^_\s]+_(?:Ter|\*)\d+ins")
# A duplicated, deleted or replaced range that ends at the stop codon holds the stop's residue,
# so it changes the stop codon as a notation that starts there does.
STOP_CODON_RANGE_END_RE = re.compile(r"[^:\s]+:p\.[^_\s]+_(?:Ter|\*)\d+(?:dup|del)")
INCOMPLETE_CDS_FLAG_RE = re.compile(r"cds_(?:start|end)_NF")
CODING_COORDINATE_RE = re.compile(r"[1-9]\d*")
ACCEPTOR_INTRON_COORDINATE_RE = re.compile(r"(\*|[1-9]\d*)-([1-9]\d*)")
DONOR_INTRON_COORDINATE_RE = re.compile(r"([1-9]\d*)\+[1-9]\d*")
# ``hgvs_shifted_deletion_clamped_span``: ``_var2transcript_slice_coords`` clamps a written span
# straddling the transcript's first base onto that base (TranscriptVariationAllele.pm:2732-2733)
# before the offset moves it, so a deletion the shift carries inside is described from the base
# after the shift where its shifted span starts earlier: two deleted bases print as one. The
# shape: two ``del`` notations of one transcript and numbering ending at the same cDNA coordinate,
# the reference's starting later than the candidate's, on a deletion's row (allele ``-``) both
# engines give a non-zero HGVS_OFFSET and, where the run prints one, a positive SHIFT_LENGTH; the
# candidate's coordinates agree with the row's cDNA position, which under that run is the shifted
# span's (an ``n.`` coordinate equal to the column's, a same-region ``c.`` span as long as it), so
# a span off by one base stays a difference. The two notations start at different coordinates,
# so the HGVSIntronOffset plugin's start offset differs with them on such a row: each engine's
# is the intron offset of its own notation's start coordinate.
HGVS_SHIFTED_DELETION_CLAMPED_SPAN = "hgvs_shifted_deletion_clamped_span"
SHIFTED_DELETION_SPAN_FIELDS = (*HGVSC_FIELDS, *HGVS_INTRON_START_OFFSET_FIELDS)
DELETION_NOTATION_RE = re.compile(r"([cn])\.(\*?-?\d+(?:[+-]\d+)?)(?:_(\*?-?\d+(?:[+-]\d+)?))?del")
# ``hgvs_clipped_reference_against_whole_alternate``: ``_var2transcript_slice_coords`` clamps a
# written span running over a transcript's edge onto the bases inside
# (TranscriptVariationAllele.pm:2732-2733) and ``hgvs_transcript`` passes that span with the whole
# alternate allele (:1474-1483), which ``hgvs_variant_notation`` reads against the reference bases
# of the span alone (Utils/Sequence.pm:526), so a substitution of several bases running over an
# edge, or a multi-allelic record's allele still carrying a shared outside base, is described as
# a ``delins`` of every alternate base over the bases inside, as the substitution, ``ins`` or
# ``delins`` the clip leaves of them (:2231-2242, :2271-2284), or, where every base inside is
# unchanged, as an ``ins`` of the outside bases at the coordinate ``_get_cDNA_position`` gives
# the base before the transcript's genomic start (:2800-2827, an intron offset off the first exon
# in genomic order). The candidate cuts the alternate as the clamp cut the reference and
# describes the change inside the transcript, none when there is none. The shape: a differing
# HGVSc on a row of a sequence allele whose cDNA position has an open side (``?``), the
# reference's notation naming the row's transcript, or empty where the whole alternate equals
# the bases inside (Utils/Sequence.pm:529-531); the candidate's of the same transcript and
# numbering, every coordinate within the row's span plus one of every coordinate of the
# reference's (both engines clip against the bases inside, so both notations lie within the
# allele's span over the transcript, give or take the base before exon 1 and the base after an
# insertion site), or empty where the change lies outside; against the empty reference, a
# ``del``; against an empty candidate, an ``ins`` or a ``dup``, what the outside bases read as
# once the clip has stripped the unchanged bases inside. An allele shorter than its span that
# deletes bases over the transcript's edge, the deletion written nowhere wholly inside the
# transcript, has no description on it, a deletion extending beyond a transcript being the
# genomic description's alone,
# and the candidate writes none: against that empty candidate the reference's is a ``del``, a
# ``delins`` or a substitution of the bases inside, what the whole alternate reads as against
# them (``deletion_over_the_edge``). The coordinate test reads one scale,
# on which the stop codon or an intron inside a notation's span misstates a distance, and an
# allele whose other end lies in an intron or past the other edge has no open side at all:
# neither end maps to a cDNA base and the position is empty on both engines
# (BaseTranscriptVariation.pm:145-148, Utils.pm:163-165), as it is on every intronic row. On a
# substitution's row the two notations show the clamp against the transcript's exon spans,
# read from the candidate cache (``whole_alternate_over_fewer_bases``, ``TranscriptEdges``):
# each ``delins`` writes a run of the row's allele, and the offset between the two runs plus
# the signed distance between the two start coordinates, and likewise at the end, is the count
# of the allele's bases that lie outside the transcript at that end; at an end inside the
# transcript the two notations name one coordinate, and each coordinate of either notation is
# written as the spans place its base, in an exon or an intron. Such a row is the class whether
# its position is open or empty; without the cache no row is.
# A deletion's row (allele ``-``) and a closed cDNA position stay differences.
# The HGVSIntronOffset plugin's two fields are the intron parts of the notation's coordinates
# (:1530-1533, :1543, :1545), so on a row whose HGVSc the class takes against a candidate that
# prints none, the reference alone writes them and the candidate neither: they are the class
# where they equal those intron parts (``clipped_notation_intron_offsets``). At a
# reverse-strand transcript's first base, its genomic end, the clip leaves Ensembl VEP an
# insertion at slice position 0, for which ``_get_cDNA_position`` returns no coordinate
# (:2834): on such a transcript's row whose cDNA position is open at its start the reference
# prints no notation where the candidate prints an ``ins`` of the bases a longer allele adds
# inside the span, and the candidate alone writes the two fields.
# The class takes the two fields of such a row, the intron parts of the candidate's notation,
# and leaves its HGVSc a difference. Where both print an ``ins`` the class takes as a pair,
# each side's fields are the intron parts of its own notation.
# Before a forward-strand transcript's first base the clip leaves Ensembl VEP an insertion
# that ends at slice position 0, the base before the first exon, which ``_get_cDNA_position``
# writes through its branch for a position between two exons as the first base's coordinate
# with ``-1`` (:2800-2826): the end offset stored beside that ``ins`` is -1 where the
# transcript has no intron. Where the candidate describes the bases such an allele adds as
# the duplication the HGVS 3' rule gives, both coordinates exonic and the bases the rule
# carried them in HGVS_OFFSET, its end offset is 0. The class takes the end offset of such a
# row where the two rows show the clamp at the first base
# (``insertion_before_the_first_base_against_duplication``); the row's HGVS_OFFSET is no
# field of the class and its HGVSc is the class's only by the coordinate test above.
HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE = "hgvs_clipped_reference_against_whole_alternate"
CLIPPED_REFERENCE_FIELDS = (*HGVSC_FIELDS, *HGVS_INTRON_OFFSET_FIELDS)
# The exon coordinate and the intron offset ``hgvs_transcript`` reads off a cDNA coordinate
# (:1530-1531): the first integer and the one after it, the star of a 3' UTR coordinate unread.
COORDINATE_PARTS_AS_READ_RE = re.compile(r"(-?[0-9]+)\+?(-?[0-9]+)?")
CLIPPED_NOTATION_RE = re.compile(
    r"(?P<numbering>[cn])\.(?P<start>\*?-?\d+(?:[+-]\d+)?)(?:_(?P<end>\*?-?\d+(?:[+-]\d+)?))?"
    r"(?:[ACGTN]>(?P<alternate_base>[ACGT])|(?P<kind>delins|ins)(?P<bases>[ACGT]+)|(?P<plain>inv|del|dup))"
)
# ``hgvs_added_bases_left_at_the_span_end``: ``_return_3prime`` shifts a record whose class is
# insertion or deletion alone (TranscriptVariationAllele.pm:120-123), so the bases added by an
# allele of any other record, one that keeps bases it shares with its reference, stay where
# ``_clip_alleles`` leaves them, at the end of the record's span, typed an ``ins`` with no
# duplication lookup (:2204-2284). Where that end is a reverse-strand transcript's last base,
# ``_get_cDNA_position`` numbers the position past it off the last base with an intron offset
# of 1 (:2802, :2820-2827): ``c.*27+1_*27insAC``, beside which the HGVSIntronOffset plugin
# writes 1 and 0. The candidate places such bases by the HGVS 3' rule where the run of points
# at which they give the same sequence reaches the transcript's first or last base: a
# duplication of the transcript's last bases where they repeat them, the notation at the most
# 3' point of a run that starts at or before the first base, with the bases the rule moved
# them in HGVS_OFFSET, and no notation where they have no placement inside the transcript.
# The shapes are the predicate's (``added_bases_left_at_the_span_end``).
HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END = "hgvs_added_bases_left_at_the_span_end"
ADDED_BASES_FIELDS = (*HGVSC_FIELDS, *HGVS_OFFSET_FIELDS, *HGVS_INTRON_OFFSET_FIELDS)
# A coordinate one base off an exonic one, as ``_get_cDNA_position`` writes the position past
# a reverse-strand transcript's last base.
ONE_PAST_COORDINATE_RE = re.compile(r"(\*?-?\d+)\+1")
# A coordinate one base before an exonic one, as ``_get_cDNA_position`` writes the position
# before a forward-strand transcript's first base.
ONE_BEFORE_COORDINATE_RE = re.compile(r"(\*?-?\d+)-1")
# ``hgvs_insertion_left_where_the_shift_stops``: ``perform_shift`` reads 1,000 bases of the
# reference either side of an insertion (TranscriptVariationAllele.pm:427-437) and walks them
# at most 1001 less the allele's length times, 1,000 times for a longer allele (:303-306);
# going 5', on a reverse-strand transcript, the limit of an allele of 1,001 bases is 0, which
# :304 leaves, and nothing is walked. An
# insertion whose run of placements goes further is described where that walk stops, beside
# an HGVS_OFFSET that is the bound itself, none where the bound is 0. The candidate walks to the end of the run, the most
# 3' placement, which the HGVS 3' rule assigns the change to: a ``dup`` or an ``ins`` further
# 3' beside a larger HGVS_OFFSET, or no notation where that placement lies past the
# transcript's last base. The HGVSIntronOffset plugin's two fields follow each notation.
HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS = "hgvs_insertion_left_where_the_shift_stops"
SHIFT_BOUND_FIELDS = (*HGVSC_FIELDS, *HGVS_OFFSET_FIELDS, *HGVS_INTRON_OFFSET_FIELDS)
# The bases of flank ``perform_shift`` reads on either side of an allele (:427-437).
SHIFT_FLANK_BASES = 1000
# ``hgvs_added_bases_written_outside_the_transcript``: OutputFactory.pm:1721 fills the HGVSc
# column only on a row within the feature, and an insertion written directly beside a
# transcript, or before its first base, is a row outside it
# (BaseVariationFeatureOverlapAllele.pm:480), so the reference prints no HGVSc and no
# HGVS_OFFSET there, while HGVSIntronOffset.pm:62 asks for the notation on every Transcript
# row and writes the pair of a notation the row never prints, or none. The description of a
# change does not turn on which of its placements the record writes: where the HGVS 3' rule
# places the inserted bases on the transcript, as a duplication of its last bases or at a
# point between two of its bases that a run from before its first base ends at, the
# candidate writes that notation on the row, with the bases the rule moved it in HGVS_OFFSET
# and the intron offsets of its own coordinates. The shape is the predicate's
# (``added_bases_written_outside_the_transcript``).
HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT = "hgvs_added_bases_written_outside_the_transcript"
# ``hgvsp_stop_codon_protein_comparison``: ``hgvs_protein`` describes a change from the peptides
# of the allele's codon window (``_get_hgvs_protein_type``, TranscriptVariationAllele.pm:2047,
# and ``_get_hgvs_protein_format``, :1898), where the HGVS nomenclature describes a protein
# change by comparing the variant protein with the reference protein. At the stop codon the
# window's description is not the nomenclature's in nine shapes, each on a row whose two
# values name one protein. In the first five both name the stop at one residue number above
# 1; the candidate writes the stop as unchanged where the protein does not change, and an
# extension named at the stop (the first new residue, the count of residues added or ``?``)
# where the whole protein is kept and residues follow it. In the sixth the candidate writes a
# frameshift, named at the first changed residue, where a change reaches from before the stop
# codon into it. In the seventh both write an extension and the candidate counts its residues.
# In the eighth the reference writes the stop as deleted before an extension's count, and the
# candidate names the residue that stands in the stop's place. In the ninth the reference
# writes the last residue and a stop as inserted before that residue, where an insertion
# brings a stop codon of its own to the stop's residue and the protein does not change, and
# the candidate writes the stop as unchanged. The fifth to the ninth shape read the row's
# terms, which are the reference's.
#   The stop replaced by a stop (``Ter<N>delinsTer``) against the stop unchanged (``Ter<N>=``).
#   The stop alone replaced by the first residues that follow the protein, with or without a
#   closing ``Ter`` (``Ter<N>delins<residues>[Ter]``), against an extension with the first of
#   them and ``?`` or a count above the residues the reference names (their number where the
#   reference closes with no ``Ter``).
#   The same where the reference names every added residue and then the stop the inserted
#   bases bring, against an extension whose count is the number of residues named.
#   Residues inserted between the last residue and the stop
#   (``<Aaa><N-1>_Ter<N>ins<residues>``) against an extension with the first of them and
#   their number.
#   The stop written as unchanged (``Ter<N>=``) on a row whose terms hold ``stop_lost``
#   against an extension.
#   A range from a residue to the stop replaced by residues, none of them a stop
#   (``<Aaa><M>_Ter<N>delins<residues>``), on a row whose terms hold ``stop_lost``, against a
#   frameshift named at that residue with the first of them (``<Aaa><M><Bbb>fsTer<count>``),
#   its count ``?``, or above their number and other than ``N - M + 1``, the count that leaves
#   the new stop at the stop's own residue.
#   An extension with an unknown count (``Ter<N><Xaa>extTer?``) on a row whose terms hold
#   ``stop_lost`` against the same extension with a count: ``_stop_loss_extra_AA``
#   (TranscriptVariationAllele.pm:2477) ends at the first stop codon of the translated bases, so a
#   stop codon inside the coding sequence that the transcript's peptide reads as a residue
#   leaves its count unknown where the 3' UTR holds a new stop.
#   The stop written as deleted with an extension count and no residue in its place
#   (``Ter<N>delextTer<count>``) on a row whose terms hold ``stop_lost``, against an
#   extension at the same residue with the same count and the first new residue.
#   One residue and a stop inserted between two neighbouring residues, the residue the right
#   flank's (``<Aaa><M>_<Bbb><M+1>ins<Bbb>Ter``), on a row whose terms hold
#   ``stop_retained_variant``, against the stop unchanged at the residue after the right
#   flank (``Ter<M+2>=``): by the reference's value the alternate protein is the reference
#   protein up to the right flank's residue and then a stop, at the residue where the
#   candidate's value has the reference protein's stop.
# A row whose FLAGS name ``cds_end_NF``, where no stop codon closes the coding sequence, is no
# row of the class. The row alone does not say that the whole reference protein is kept, that
# a change runs to the stop codon's last base or shifts the reading frame, that the coding
# sequence holds a stop codon its peptide reads as a residue, nor what a count is that the
# reference does not spell out; the same notation with the sides exchanged, another
# residue number or first residue, a range ending at the stop on a row whose terms hold no
# ``stop_lost``, and a reference that is empty or writes the stop as deleted with no extension
# count stay differences. For the ninth shape the reference's value names the two flanks
# alone, so the stop's residue number is the candidate's, held to the one after the right
# flank. More than one residue before the inserted stop stays a difference, the row not
# naming the reference residues past the right flank; so do an inserted residue that is not
# the right flank's, flanks that are not two neighbouring residues, and a row whose terms
# hold no ``stop_retained_variant``.
HGVSP_STOP_CODON_PROTEIN_COMPARISON = "hgvsp_stop_codon_protein_comparison"
FLAGS_FIELDS = ("FLAGS", "flags")
RESIDUE_NAME = r"[A-Z][a-z]{2}"
STOP_UNCHANGED_RE = re.compile(r"Ter([1-9]\d*)=")
STOP_EXTENSION_RE = re.compile(rf"Ter([1-9]\d*)({RESIDUE_NAME})extTer([1-9]\d*|\?)")
STOP_REPLACED_RE = re.compile(rf"Ter([1-9]\d*)delins((?:{RESIDUE_NAME})+)")
STOP_INSERTION_RE = re.compile(rf"{RESIDUE_NAME}([1-9]\d*)_Ter([1-9]\d*)ins((?:{RESIDUE_NAME})+)")
STOP_RANGE_REPLACED_RE = re.compile(rf"({RESIDUE_NAME})([1-9]\d*)_Ter([1-9]\d*)delins((?:{RESIDUE_NAME})+)")
STOP_FRAMESHIFT_RE = re.compile(rf"({RESIDUE_NAME})([1-9]\d*)({RESIDUE_NAME})fsTer([1-9]\d*|\?)")
STOP_DELETED_EXTENSION_RE = re.compile(r"Ter([1-9]\d*)delextTer([1-9]\d*|\?)")
RESIDUE_AND_STOP_INSERTED_RE = re.compile(rf"({RESIDUE_NAME})([1-9]\d*)_({RESIDUE_NAME})([1-9]\d*)ins({RESIDUE_NAME})Ter")
# ``hgvsp_stop_named_on_incomplete_coding_sequence``: ``hgvs_protein`` three-letters the ``X`` of
# a partial codon and of an unknown residue as ``Xaa`` and rewrites it ``Ter``
# (``s/Xaa/Ter/g``), so on a coding sequence that no stop codon closes, and at residue 1 of
# one that opens inside its first codon, it names a stop the reference protein does not hold.
# The HGVS nomenclature takes as a reference only a transcript with a complete coding
# sequence, and the candidate writes no HGVSp there. A row of the class: the reference's
# description names the stop alone, as unchanged (``Ter<N>=``), deleted (``Ter<N>del``),
# replaced by one residue (``Ter<N><Xaa>``) or extended (``Ter<N><Xaa>extTer<count>``); the
# candidate's HGVSp is empty; and the row shows the incomplete coding sequence: its FLAGS name
# ``cds_end_NF``, its terms hold ``incomplete_terminal_codon_variant``, or the description is
# ``Ter1=`` and its FLAGS name ``cds_start_NF``. The row alone does not show a coding sequence
# that no flag marks incomplete: such a row stays a difference, and so do any other reference
# description, a candidate that writes a notation and the sides exchanged.
HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE = "hgvsp_stop_named_on_incomplete_coding_sequence"
STOP_NAMED_ALONE_RE = re.compile(rf"Ter([1-9]\d*)(=|del|{RESIDUE_NAME}(?:extTer(?:[1-9]\d*|\?))?)")
# ``hgvsp_intron_end_insertion``: under --hgvs, an insertion the input writes between an intron's
# last base and the first base of the next exon in the transcript's direction, that base inside
# the coding region, which the 3' rule leaves in that gap. Ensembl VEP reads the inserted bases
# as the exon's first bases and writes the HGVSp of that reading. The candidate writes none,
# because the HGVS nomenclature gives a change at an intron's end no predicted protein
# consequence. Both engines write the same HGVSc, which names the form, ``N`` the exon's first
# coding base: a copy of the intron's last base (``c.N-1dup``); a duplication of two or more
# intron positions ending on that base (``c.N-K_N-1dup``, or ``c.M+J_N-1dup`` with ``M`` the
# coding base before ``N`` where the copy starts in the intron's first half); an insertion
# after that base (``c.N-1_Nins``) of one base, or of bases that do not end in one of the two
# pairs that end an intron; or an insertion there of two or more bases that end in such a pair.
# Every ``c.N-1_Nins`` insertion is therefore of a form. The candidate keeps the notation at the
# end of an intron of 13 bases or fewer and on a row the 3' rule moved, which carries an
# HGVS_OFFSET, so an HGVSp that differs there stays a difference. The shape: an HGVSp the
# reference fills and the candidate leaves empty, on a row whose HGVSc is one notation of the
# row's transcript on both sides, of one of the four forms, beside no HGVS_OFFSET, the record
# written in such a gap of that transcript. A record that carries another allele beside the
# insertion holds it with the reference base it shares, the allele that base and the inserted
# bases: both engines write the ``c.N-1_Nins`` notation of the allele clipped of that base, and
# the class takes the row by the gap the clipped insertion lies in (``held_insertion_gaps``).
# A text row gives such a record's Location, not the gap, and spells its allele on the record's
# strand, so the class places it by the Location's span: the gaps its allele's clippings leave
# on either strand, or, where the allele is the insertion alone (a minimised allele), every gap
# of the span (``lone_insertion_gaps``).
# A span of two adjacent positions with a minimised allele is an insertion written alone or a
# reference allele of two bases, so it names its three gaps.
# A VCF or JSON entry spells its allele in more ways than a text row does, so the class places it
# from its input line, a VCF line or an Ensembl-format one, by the notation alone: an ALT that
# clips against the REF to a pure insertion names a gap, which the entry claims where the
# inserted bases are the notation's or as many as the positions a ``dup`` notation duplicates
# (``notation_insertion_gaps``).
# Neither output carries the intron or the coding bounds, so the class reads them from the
# candidate's JSON cache (``IntronEnds``, ``--json-cache``) and accepts nothing without it.
HGVSP_INTRON_END_INSERTION = "hgvsp_intron_end_insertion"
INTRON_END_COPY_RE = re.compile(r"c\.([1-9]\d*)-1dup")
INTRON_END_DUPLICATION_RE = re.compile(r"c\.([1-9]\d*)([+-])([1-9]\d*)_([1-9]\d*)-1dup")
INTRON_END_INSERTION_RE = re.compile(r"c\.([1-9]\d*)-1_([1-9]\d*)ins([A-Z]+)")
# The two pairs of bases that end an intron, of the major kind and of the minor kind, as
# character codes, on the transcript's strand as an HGVSc spells inserted bases.
INTRON_END_PAIRS = frozenset(bytes(pair).decode() for pair in ((65, 71), (65, 67)))
# The longest intron at whose end the candidate keeps the notation: Ensembl VEP reads an intron
# of this many bases or fewer with the exons beside it.
SHORT_INTRON_BASES = 13
# ``hgvsp_intron_start_insertion``: under --hgvs, an insertion the 3' rule leaves between a
# coding exon's last base and the first base of the intron after it, written in that gap or
# inside the exon and carried to it, whose bases are no copy of the bases before the gap.
# Ensembl VEP reads the inserted bases as the exon's last bases and writes the HGVSp of that
# reading. The candidate writes none, because the HGVS nomenclature gives a change at an
# intron's start no predicted protein consequence. A copy of the bases before the gap is a
# duplication of coding positions and keeps its notation on both engines. Both engines write
# the same HGVSc, ``c.N_N+1ins`` with ``N`` the exon's last coding base in both coordinates,
# and, on a record the 3' rule moved, the same HGVS_OFFSET: the bases moved, positive on a
# forward-strand transcript and negative on a reverse-strand one, so the gap the insertion is
# left in is the gap the record is written in moved by the offset as signed. A run that prints
# the moved location (``--shift_genomic 1``) names the gap the insertion is left in, so the
# class accepts a Location's gap as named or moved by the offset. The intron's length plays no
# part. The shape: an HGVSp the reference fills and the candidate leaves empty, on a row whose
# HGVSc is one notation of the row's transcript on both sides, of that form, beside no
# HGVS_OFFSET or one non-zero integer on both sides, the gap the insertion is left in one in
# which an exon of that transcript ends inside the coding region before an intron. A
# notation on both sides, an HGVSc of another form and offsets that differ stay differences.
# A record that holds the insertion with a reference base it shares carries no offset and is
# taken by the gap its clipped allele lies in (``held_insertion_gaps``).
# A text row of such a record is placed by its Location's span, as in the intron-end class; a
# minimised allele the 3' rule moved is the notation's bases turned back by the offset
# (``lone_insertion_gaps``).
# A run under ``--shift_3prime 1`` prints a moved record's allele as the rule moved it, the
# notation's bases as they stand, so on a row that carries the offset the class reads the bases
# either way, a text row's minimised allele against the notation's with its span's gaps as named
# or moved and, for a VCF or JSON entry, the bases an ALT of its input line inserts against the
# notation's turned back by the offset (``notation_insertion_gaps``).
# The class reads the exons and the coding bounds from the candidate's JSON cache
# (``IntronEnds``, ``--json-cache``) and accepts nothing without it.
HGVSP_INTRON_START_INSERTION = "hgvsp_intron_start_insertion"
INTRON_START_INSERTION_RE = re.compile(r"c\.([1-9]\d*)_([1-9]\d*)\+1ins([A-Z]+)")
# The classes that read the candidate cache's intron gaps (``IntronEnds``).
INTRON_GAP_CLASSES = (HGVSP_INTRON_END_INSERTION, HGVSP_INTRON_START_INSERTION)
# The version a stable id may end in, which a row's Feature carries under --transcript_version
# and a cache's stable id of an Ensembl transcript does not.
VERSION_SUFFIX_RE = re.compile(r"\.\d+$")
# A cache shard's file name: the first and last position of its window.
CACHE_SHARD_NAME_RE = re.compile(r"(\d+)-(\d+)\.json(?:\.gz)?")
# ``shift_intron_term_unshifted_span``: under --shift_3prime 1, ``_intron_effects``
# (BaseTranscriptVariationAllele.pm:99-150) sets ``intronic`` when the shifted span or the
# unshifted span (:123-124, tested at :146-147) lies in an intron's interior, while every other
# key of the method reads the shifted span alone; an insertion or deletion in an intron that
# the shift carries onto the two bases of the acceptor site, or past them onto the exon's first
# bases, keeps ``intron_variant`` beside the shifted span's terms (an intron is tested for a
# span beginning within three bases past its end, the margin of the intron interval tree,
# BaseTranscriptVariation.pm:1027). The candidate reads the shifted span alone, as Ensembl VEP
# does for the allele written at the shifted position. The shape: a Consequence set the
# reference writes as the candidate's plus ``intron_variant``, on a row whose SHIFT_LENGTH is
# positive and whose HGVSc places the shifted span at the acceptor's two bases or in the exon
# (every coordinate exonic or at the intron offset -1 or -2, a duplication's first excepted,
# since its bases run back into the intron; a notation without an offset on a row without a
# splice-site term, since a deletion over a whole intron is numbered over the exons alone and
# covers both sites), so that no intron base away from the splice site is touched; a deletion
# reaching such a base carries the term from its shifted span on both engines. A row without
# an HGVSc is read by its allele and terms: an insertion's shifted span lies between the two
# acceptor bases whenever the row carries ``splice_acceptor_variant``.
SHIFT_INTRON_TERM_UNSHIFTED_SPAN = "shift_intron_term_unshifted_span"
CONSEQUENCE_FIELDS = ("Consequence", "consequence_terms")
INTRON_TERM = "intron_variant"
SPLICE_ACCEPTOR_TERM = "splice_acceptor_variant"
SPLICE_SITE_TERMS = frozenset((SPLICE_ACCEPTOR_TERM, "splice_donor_variant"))
INDEL_NOTATION_RE = re.compile(r"[cn]\.(\*?-?\d+(?:[+-]\d+)?)(?:_(\*?-?\d+(?:[+-]\d+)?))?(dup|del|ins[ACGT]+)")
# The intron offset a coordinate ends with: the sign follows the position's digits, so the
# leading minus of a 5' UTR position (``c.-3``) is none.
INTRON_OFFSET_RE = re.compile(r"(?<=\d)[+-]\d+$")
ACCEPTOR_OFFSETS = ("-1", "-2")
# ``mature_mirna_unshifted_span``: under --shift_3prime 1, the predicates read the shifted span
# from ``$bvf->{start}`` and ``{end}``, which BaseVariationFeatureOverlapAllele.pm:269-282 offset
# for each predicate call; ``within_mature_miRNA`` (Utils/VariationEffect.pm:532-558) overlaps
# the mature miRNA's genomic span with ``$bvf->seq_region_start`` and ``seq_region_end``
# (:545-549) instead, which BaseVariationFeature.pm:279-282 and :304-307 cache at their first
# read, before any shift, so the test reads the unshifted span. An insertion or deletion the
# shift carries off the mature sequence keeps ``mature_miRNA_variant`` where the shifted span
# is ``non_coding_transcript_exon_variant`` (``non_coding_exon_variant``, :502-519, reads the
# shifted span and yields to the mature test), and one carried onto it keeps the exon term. The
# candidate reads the shifted span alone, as Ensembl VEP does for the allele written at the
# shifted position. The shape: two Consequence sets equal but for one engine's
# ``mature_miRNA_variant`` against the other's ``non_coding_transcript_exon_variant``, on a row
# whose SHIFT_LENGTH is positive, which only an insertion or deletion is given.
MATURE_MIRNA_UNSHIFTED_SPAN = "mature_mirna_unshifted_span"
MATURE_MIRNA_TERM = "mature_miRNA_variant"
NON_CODING_EXON_TERM = "non_coding_transcript_exon_variant"
# ``stop_codon_unshifted_span``: ``_overlaps_stop_codon_cil`` (Utils/VariationEffect.pm
# :1389-1430) overlaps the stop codon's genomic span with ``$bvf->seq_region_start`` and
# ``seq_region_end`` (:1403, :1424-1427), the cached unshifted coordinates, where every other
# test of ``stop_retained`` (:1290-1323) reads the shifted span. A deletion the shift carries
# onto the stop codon from the bases before it therefore fails the overlap: ``stop_retained`` is
# false, and ``coding_unknown`` (:1616-1649) writes ``coding_sequence_variant`` for a span whose
# shifted cDNA start is undefined. The candidate reads the shifted span, as Ensembl VEP does for
# the allele written at the shifted position. The shape: two Consequence sets equal but for the
# reference's ``coding_sequence_variant`` against the candidate's ``stop_retained_variant``, on
# a row whose SHIFT_LENGTH is positive. The other direction cannot arise: a span the 3' shift
# carries off the stop codon lies in the 3' UTR, where the term is not a candidate.
STOP_CODON_UNSHIFTED_SPAN = "stop_codon_unshifted_span"
# The classes of the stop codon tests behind ``stop_retained`` (Utils/VariationEffect.pm
# :1290-1324); the comment above each one's predicate states its shape.
STOP_RETAINED_WITHOUT_STOP_CODON_TEST = "stop_retained_without_stop_codon_test"
STOP_CODON_OVERLAP_GENOMIC_SPAN = "stop_codon_overlap_genomic_span"
STOP_RETAINED_THREE_BASE_STOP_CODON_TEST = "stop_retained_three_base_stop_codon_test"
STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST = "stop_codon_removed_without_stop_codon_test"
STOP_CODON_OVERLAP_INSERTION_RANGE = "stop_codon_overlap_insertion_range"
# ``splice_term_past_the_record_span``, the class of the six terms an intron gives a row: the
# five splice terms whose windows lie inside it and the intron term, whose window it is
# (VariationFeatureOverlapAllele.pm:405-460, BaseTranscriptVariationAllele.pm:117-118 and
# :143-204). On a record whose alternate allele is longer than its reference allele, Ensembl
# VEP writes a splice donor, splice acceptor, donor fifth base, donor region, polypyrimidine
# tract or intron term from positions past the end of the record's span; where no base of the
# span lies in a window of that term, vep-rs does not write it. The comment above the
# predicate in ``compare_vep_outputs.py`` states the shape with the record's alleles and
# without them.
SPLICE_TERM_PAST_THE_RECORD_SPAN = "splice_term_past_the_record_span"
# ``refseq_edited_model_deleted_span_term``: on a model whose alignment to the genome deleted
# bases (``BAM_EDIT=OK`` with a deletion edit), TranscriptMapper.pm:168-221 pairs no cDNA
# position with the deleted bases (:213-218), so a sequence allele wholly inside them maps to
# a Bio::EnsEMBL::Mapper::Gap: ``cds_coords`` and ``cdna_coords`` (BaseTranscriptVariation.pm
# :478-528) hold the gap alone, ``within_cds`` and ``within_cdna`` (Utils/VariationEffect.pm
# :636-693) need a Coordinate or a frameshift intron and find neither, no predicate is true,
# and BaseVariationFeatureOverlapAllele.pm:285 assigns the default term, ``intergenic_variant``
# (Utils/Constants.pm:636-649, whose ``within_feature => 0`` condition the default never
# consults), on the transcript's row although ``_bvfo_preds`` (:454-558) marked the allele
# within the feature, in an exon and coding. The candidate writes the term of the transcript
# region whose genomic span holds the allele: ``coding_sequence_variant`` inside the coding
# sequence, the UTR term inside a UTR. The shape: a Transcript row of a sequence allele
# carrying ``BAM_EDIT`` ``OK`` on both sides whose Consequence is ``intergenic_variant`` in the
# reference and exactly one of the three region terms in the candidate, neither row with a
# cDNA position and both with one IMPACT. A structural allele's row (a breakend mate within
# --distance takes the same default), a row with a cDNA position (an allele reaching a mapped
# base is ``coding_sequence_variant`` with a position on both engines) and any other candidate
# term stay differences.
REFSEQ_EDITED_MODEL_DELETED_SPAN_TERM = "refseq_edited_model_deleted_span_term"
IMPACT_FIELDS = ("IMPACT", "impact")
CDNA_POSITION_FIELDS = ("cDNA_position", "cdna_start", "cdna_end")
INTERGENIC_TERM = "intergenic_variant"
DELETED_SPAN_REGION_TERMS = frozenset({"coding_sequence_variant", "5_prime_UTR_variant", "3_prime_UTR_variant"})
# ``refseq_edited_model_tail_term``: on a model whose edits inserted bases (``BAM_EDIT=OK``
# with an insertion edit), the cached mapper (TranscriptMapper.pm:168-221, :207-211) numbers
# every base after an insertion in the edited frame, so the model's last mapped bases carry
# cDNA coordinates above the sum of its exon lengths; ``within_cdna`` (Utils/VariationEffect.pm
# :670-694, line 680) requires the coordinate to be at most ``$feat->length``, and
# ``Transcript::length`` (the core API's Transcript.pm:1666-1682) adds the edits' lengths only
# under ``edits_enabled``, which the cached transcript object does not carry, so it returns the
# unedited sum; ``within_3_prime_utr`` (:744-756) is false for as many final bases as the edits
# inserted, no predicate is true, and BaseVariationFeatureOverlapAllele.pm:285 assigns the
# default term, ``intergenic_variant``, on a row that prints the cDNA position the same
# coordinate gave. The candidate writes ``3_prime_UTR_variant``. The shape: a Transcript row of
# a sequence allele carrying ``BAM_EDIT`` ``OK`` on both sides whose Consequence is
# ``intergenic_variant`` in the reference and ``3_prime_UTR_variant`` in the candidate, a cDNA
# position on both and one IMPACT. A row without a cDNA position (the deleted-span class's
# shape), a structural allele's row and any other candidate term stay differences; a
# non-coding model's exon term reads the exon's genomic span, not the mapper, so it has no
# such row.
REFSEQ_EDITED_MODEL_TAIL_TERM = "refseq_edited_model_tail_term"
THREE_PRIME_UTR_TERM = "3_prime_UTR_variant"
# ``refseq_edited_model_deleted_span_codons``: TranscriptVariationAllele.pm:872-874 splices the
# reference allele's bases over the CDS span of a model carrying ``_rna_edit`` before reading
# the reference codons (:877), ``substr($cds, $cds_start - 1, $vf_nt_len) = $seq``, which
# assumes the allele has the span's length. For an allele spanning bases the model's alignment
# deleted, TranscriptVariation.pm:459 refuses the model's own bases as the reference allele
# (the span maps to three segments) and the genomic bases stand, more of them than the CDS
# span holds; the splice lengthens the CDS and the codon window reads the genomic bases from
# the span's start, among them bases the model lacks, while the model's own final bases of the
# span fall out of the window. The candidate reads the model's codons. The shape: the Codons
# of a Transcript row of a sequence allele carrying ``BAM_EDIT`` ``OK`` on both sides whose
# USED_REF is longer than its CDS span, the two values' alternate codons equal, their reference
# codons of one length agreeing at every lower-case (unchanged) base and differing in an
# upper-case one, the reference's upper-case bases spelling USED_REF's leading bases on the
# transcript's strand. A row whose allele has the span's length (the splice replaces the
# model's bases with themselves), a differing alternate codon, a differing amino acid and
# upper-case bases that are not the allele's stay differences.
REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS = "refseq_edited_model_deleted_span_codons"
CODONS_FIELDS = ("Codons", "codons")
USED_REF_FIELDS = ("USED_REF", "used_ref")
# The fields that state a row's coding change beside its consequence set; a row of
# ``organelle_peptide_paths_translation_table`` carries each alike on both sides.
CODING_CHANGE_FIELDS = (
    *CDNA_POSITION_FIELDS, *CDS_POSITION_FIELDS, "Protein_position", "protein_start", "protein_end",
    "Amino_acids", "amino_acids", *CODONS_FIELDS, *HGVSC_FIELDS, *HGVS_OFFSET_FIELDS,
)
# ``frequency_filter_zero_value``: BaseCacheVariation.pm:419 at release/116.2 reads an allele's
# frequency with ``|| 'NA'``, so a frequency of exactly 0, Perl-false, reads ``NA`` and fails
# --check_frequency and --filter_common; the candidate judges the zero, as Ensembl VEP 116.3
# does. Two shapes: a ``FREQS`` whose items agree but for a value the reference reads ``NA``
# and the candidate ``0``, and under ``--freq_filter include`` the rows of a record the
# candidate alone keeps, every ``FREQS`` item of which carries the value ``0``, with the
# record's VCF line and JSON object. A JSON consequence entry carries no ``FREQS``, so in that
# format the row shape applies from a manifest's candidate-only tuples alone.
FREQUENCY_FILTER_ZERO_VALUE = "frequency_filter_zero_value"
FREQS_FIELDS = ("FREQS", "freqs")
# ``merged_same_stable_id_collapse``: under a RefSeq or merged source AnnotationType/Transcript.pm
# :313-343 at release/116.2 keeps one transcript per stable identifier among the models a
# buffer's regions load (the ``ensembl``-sourced copy, else the lowest dbID), where the method's
# first pass and its description uniquify by dbID; a model the cache carries at two loci keeps
# the copy listed first and loses every row at the other, whether or not a record overlaps the
# kept copy. The candidate annotates both copies. The shape: a Transcript entry the candidate
# alone writes on a record the reference wrote a Transcript row for, whose SOURCE, where the
# run writes one, is RefSeq (a merged cache's Ensembl-sourced models are unique per locus; a
# RefSeq cache writes no SOURCE), of a model the reference wrote no row for on that record; the
# JSON record's most_severe_consequence is left uncompared as for a documented key. A record
# whose every overlapping model is a lost copy gets an intergenic row from the reference, a
# shape outside the class, as are another allele's row of a kept model and a duplicated row.
MERGED_SAME_STABLE_ID_COLLAPSE = "merged_same_stable_id_collapse"
SOURCE_FIELDS = ("SOURCE", "source")
REFSEQ_SOURCE = "RefSeq"
# ``hgnc_id_copy``: AnnotationType/Transcript.pm:290 at release/116.2 records the HGNC id of
# every transcript carrying one under its gene symbol and :305 lends it to every same-symbol
# transcript among the features fetched for one buffer, so a transcript whose gene is named
# outside HGNC (EntrezGene, RFAM; every RefSeq transcript of a merged cache) prints an id in a
# buffer that also loaded the same-named HGNC gene and nothing in one that did not. The
# candidate lends the id once, as the cache is built. The shape: the HGNC_ID of a Transcript
# entry whose SYMBOL_SOURCE is not HGNC, one side empty and the other an id (``HGNC:10528``,
# or ``10528`` in a cache before the prefix), in either direction. Only a run that writes the
# field compares it: ``--symbol``, or the VCF format, whose option set implies it
# (``VCF_IMPLIED_COLUMNS``).
HGNC_ID_COPY = "hgnc_id_copy"
HGNC_ID_FIELDS = ("HGNC_ID", "hgnc_id")
SYMBOL_SOURCE_FIELDS = ("SYMBOL_SOURCE", "gene_symbol_source")
HGNC_SOURCE = "HGNC"
HGNC_ID_RE = re.compile(r"(?:HGNC:)?\d+")
# The two field classes of the Downstream plugin's ``DownstreamProtein``. Downstream.pm:108-110
# write the allele over the change's first to last coding base, :112-113 round the first base's
# coordinate, counted from 1, down to a multiple of three, :115-119 read the changed coding
# sequence from that offset, counted from 0, and follow it with the 3' UTR, and :139-142
# translate the result and cut it at its first stop. For a deletion, with or without inserted
# bases, whose first base is a codon's third base the offset is the base after it, so the
# translation starts one codon after the first codon the change touches; the candidate starts
# at the touched codon. The rows do not state the sequences, so both classes read them from the
# candidate JSON cache (``DownstreamSequences``, ``--json-cache``) and recompute both values
# (``downstream_class``); without the cache they accept nothing. The shape the two share: a
# DownstreamProtein both sides carry (an empty value is a value, which the plugin writes when
# the first codon read is a stop; ``-`` or no key is none) on a row whose two entries have one
# consequence set, or two that differ by the reference's ``start_lost`` as the registry's
# ``start_cooccurrence_swap`` admits (``registry_class``), with ``frameshift_variant`` and no
# splice term (:95-99), one CDS_position (or the two of a later allele's row,
# ``later_allele_prints_shifted_pair``) and one ProteinLengthChange; the cache holds the
# row's transcript with a
# ``translateable_seq``; and on the record's own coordinates, the CDS_position as printed, the
# allele read on the transcript's strand (``-`` as no bases, any other reverse-complemented
# where the cache's ``strand``, or the row's STRAND for a transcript the cache gives none, is
# -1), the first base is a codon's third base, the reference's value is the changed sequence
# translated under the transcript's codon table (:130-133, ``GENETIC_CODES``) from the offset
# with the UTR and the candidate's the same translated from three bases before it. Under
# --hgvs the plugin holds the coordinates the protein notation shifted
# (``plugin_hgvs_shifted_cds_coordinates`` below): a row is one of these two classes' where
# the plugin's reading of the shifted coordinates gives the value its reading of the record's
# own gives, and that class's where the two differ. A value one side lacks, another value on
# either side, a first base that is a codon's first or second, a differing
# ProteinLengthChange, a transcript under a codon table ``GENETIC_CODES`` lacks and a
# transcript the cache lacks stay differences.
#
# ``downstream_utr_frame_past_coding_end``: the offset lies past the changed sequence's end, a
# deletion without inserted bases (allele ``-``) that runs to the last base of the
# ``translateable_seq``. The changed sequence is one base shorter than the offset, the
# ``substr`` of :115 returns nothing, and the 3' UTR is translated from its first base, one
# base out of the changed sequence's frame, where the candidate translates the two coding
# bases before the deletion with the UTR. The transcript has a 3' UTR: without one both
# values are empty.
#
# ``downstream_start_past_own_codon``: the offset lies inside the changed sequence or at its
# end, every other deletion and every deletion with inserted bases of the shared shape; the
# reference's value is the candidate's less the touched codon's residue, or the residues past
# a stop the touched codon spells where the candidate's value is empty. A transcript without
# a 3' UTR is read with an empty one. An insertion between two bases the model keeps adjacent
# is never the class's: its CDS_position prints the two bases around it, and read as a
# replacement of those two bases its recomputed values are neither engine's. An insertion on
# a seam of the model, two mapper pairs that touch on the genome with model bases between
# them, prints those bases as a span of its own (Mapper.pm:504-558), which is read as any
# other span. The two classes share no row.
DOWNSTREAM_UTR_FRAME_PAST_CODING_END = "downstream_utr_frame_past_coding_end"
DOWNSTREAM_START_PAST_OWN_CODON = "downstream_start_past_own_codon"
DOWNSTREAM_CLASSES = (DOWNSTREAM_UTR_FRAME_PAST_CODING_END, DOWNSTREAM_START_PAST_OWN_CODON)
DOWNSTREAM_PROTEIN_FIELDS = ("DownstreamProtein", "downstreamprotein")
PROTEIN_LENGTH_CHANGE_FIELDS = ("ProteinLengthChange", "proteinlengthchange")
# ``plugin_hgvs_shifted_cds_coordinates``: under --hgvs a plugin reads the CDS coordinates the
# protein notation left shifted on the row. OutputFactory.pm:1737 prints CDS_position from the
# record's own span and :1760 calls ``hgvs_protein``, whose ``codon`` deletes the row's kept
# CDS pair and reads it again at the 3'-shifted position (TranscriptVariationAllele.pm:844-858);
# ``_get_alternate_cds`` reads it there once more and leaves it there (:2406-2410), and
# ``shift_feature_seqs`` rotates the allele's sequence for the shifted position (:1330-1345,
# called at :1778). ``hgvs_protein`` returns before all of it when the shifted span has no
# translation coordinates (:1730-1738). Nothing clears the pair before OutputFactory.pm:398
# runs the plugins, so Downstream.pm:108 and :110 read the shifted pair and the rotated
# allele, and NMD.pm:88, :98 and :103 the shifted pair (:137 reads the record's own genomic
# end). Without --hgvs the same plugins read the record's own pair, which the candidate reads
# with and without the flag. The class reads the candidate JSON cache as the two Downstream
# classes do (``DownstreamSequences``, ``--json-cache``), with the transcript's strand, exons,
# coding bounds and mapper pairs, and accepts nothing without it.
#
# The two pairs of a row whose reference entry carries a non-zero HGVS_OFFSET: the own pair
# is the CDS_position, on a later allele's row whose reference entry prints the shifted pair
# the candidate's (``later_allele_prints_shifted_pair``: the alleles of a record share one
# TranscriptVariation, and OutputFactory.pm:1737 prints on a later allele's row the pair an
# earlier allele's notation left there, the span moved by the row's own offset or by the
# offset of the last earlier allele's row that carries one, ``earlier_allele_shift``); the
# shifted pair is the record's genomic span moved by the offset as
# printed (signed along the genome) and mapped through the mapper pairs
# (``TranscriptModel.cds_pair``: BaseTranscriptVariation.pm:252-264, and Mapper.pm:504-558 for
# an insertion, whose pair runs from the base after its 5' flank to the base before its 3'
# flank). The record's span is the Location of a text entry, which must map to the row's
# CDS_position; a VCF or JSON entry carries none, so its span is the one the CDS pair maps
# back to, and it has none where the pair holds a gap or a base outside the mapper pairs or
# is an insertion's between two bases the genome does not keep adjacent. The shift reaches
# the plugin where both ends of the shifted span are CDS coordinates; the allele it reads
# there is the Allele turned left one base as many times as the offset's absolute value on a
# forward-strand transcript, as its length less that value on a reverse-strand one (not at
# all when that is not positive), then read on the transcript's strand. An insertion and a
# two-base replacement print the same two adjacent numbers in Location and CDS_position, so a
# sequence allele at such a span is read both ways, a reading counting only where it
# recomputes both entries' values.
#
# The Downstream shape: a DownstreamProtein or a ProteinLengthChange both entries carry, on a
# row of the two Downstream classes' consequence set and CDS_position, where the candidate's
# value is the candidate's rule on the own pair with the Allele, the reference's is
# Downstream.pm:104-149 on the shifted pair with the rotated allele, and that differs from
# Downstream.pm:104-149 on the own pair, the value Ensembl VEP writes without the flag. A row
# whose two readings of Downstream.pm give one value is left to the two classes above, so the
# three share no row. Where the candidate's value is not that own-pair value either, the own
# pair has the shape of one of the two classes above and the start rule changes the value as
# well as the flag: the row is this class's only where that class is declared beside it, and
# with either of the two declared alone it stays a difference. The NMD shape: an ``NMD_escaping_variant`` one entry carries and the
# other lacks, on a row whose two entries have one consequence set, or two that differ as
# ``stop_retained_without_stop_codon_test_pair`` admits, each naming one of the plugin's four
# terms (NMD.pm:65), and one CDS_position or a later allele's two, where the reference carries it exactly when the rule of
# NMD.pm:98-109 holds on the shifted pair (both ends coordinates, and the CDS end at most 101
# or the record's genomic end in the last exon or within 51 bases of the 3' end of the second
# to last, ``variant_exon_check``, :132-176), the candidate exactly when it holds on the own
# pair, and the two differ. :90-91 count the one array reference ``get_all_Introns`` returns,
# so the rule of :109 for a transcript without introns never holds, and a value the candidate
# writes by it alone is outside the class. A later allele's row without an HGVS_OFFSET whose
# earlier allele's row carries one is the Downstream shape with the pair that allele's
# notation left, the pair the reference prints, in the shifted pair's place and the Allele
# as written in the rotated allele's. The same values on any other row without an HGVS_OFFSET,
# an NMD on a row without one, with a shift that does not reach the plugin or with a shifted pair that gives the own
# pair's value, a value planted on either side, a Downstream value one side lacks, an NMD
# both entries carry or both lack, and a transcript the cache lacks or holds without its
# mapper pairs stay differences.
PLUGIN_HGVS_SHIFTED_CDS_COORDINATES = "plugin_hgvs_shifted_cds_coordinates"
NMD_FIELDS = ("NMD", "nmd")
NMD_ESCAPING = "NMD_escaping_variant"
NMD_GATE_TERMS = frozenset(("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"))
# The classes that read the candidate cache's transcripts (``DownstreamSequences``).
TRANSCRIPT_READING_CLASSES = (*DOWNSTREAM_CLASSES, PLUGIN_HGVS_SHIFTED_CDS_COORDINATES)
# ``plugin_value_on_reference_skipped_record``: a field of DosageSensitivity or UTRAnnotator
# that the candidate alone writes on a row of a record the reference dropped
# (``reference_skipped_record``). Ensembl VEP writes no row for a record it marks skipped
# (Parser.pm:491-501, :821-827), so no plugin runs on it; the candidate annotates the record
# and its plugins write on those rows as on any structural row. The shape is the base class's
# own: a Location at which the candidate wrote a structural row and the reference none, no
# candidate row there naming a record the reference wrote (``ReferenceRecordScope`` and
# ``reference_named_locations`` of ``compare_vep_outputs.py``). A row of such a record has no
# reference entry, so this comparison sets the row aside whole under the structural-variant
# class and ``field_class_accepts`` takes no pair for the class; the plugin tuple scorer,
# which scores a value and not a row, decides it from the rows of the value's Location. A
# value on the reference's side, a value at a Location the reference wrote a row for and a
# value of a tuple-file engine, which states no rows, stay differences.
PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD = "plugin_value_on_reference_skipped_record"
DOSAGE_SENSITIVITY_FIELDS = ("pHaplo", "pTriplo")
UTRANNOTATOR_EXISTING_FIELDS = ("Existing_InFrame_oORFs", "Existing_OutOfFrame_oORFs", "Existing_uORFs")
# The keys of one pair of a cache transcript's ``mapper``: a run of cDNA coordinates, the
# genomic bases it maps to and their orientation.
MAPPER_PAIR_KEYS = ("from_start", "from_end", "to_start", "to_end", "ori")
# The genetic codes Bio::Tools::CodonTable translates with, by table id, as
# ``crates/vep-core/src/codon.rs`` ``CODES`` lists them: each codon's amino acid, the codons in
# the order TTT TTC TTA TTG TCT ... GGG. Ids 7, 8 and 17 to 20 hold no table there and ids from
# 25 on are unknown to it, so a transcript under one of them has no code here.
CODON_TABLE_AMINO_ACIDS = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    24: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
CODON_ORDER = tuple(first + second + third for first in "TCAG" for second in "TCAG" for third in "TCAG")
# Each table's amino acid by codon, and the bases each nucleotide letter stands for.
GENETIC_CODES = {table: dict(zip(CODON_ORDER, amino_acids)) for table, amino_acids in CODON_TABLE_AMINO_ACIDS.items()}
NUCLEOTIDE_READINGS = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT",
    "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT", "X": "ACGT",
}
# The letter Bio::Tools::CodonTable writes for a codon whose readings code two amino acids.
AMBIGUOUS_RESIDUE_PAIRS = {frozenset("DN"): "B", frozenset("EQ"): "Z"}
# The three classes of the structural-variant comparator (``scripts/validation/
# compare_sv_concordance.py``), built in so that a structural-variant case sets aside what the
# adjusted structural-variant figures set aside: each is decided by that comparator's own
# predicate over the entries read as its tuples (Location, Allele, Feature, Feature_type,
# Consequence set), a VCF or JSON entry's Location the one its input record gives
# (``vep_rs_record_rows``). ``reference_skipped_record`` (``filter_reference_skipped_records``):
# a row the candidate alone writes on a record Ensembl VEP drops before its input buffer, an
# unsupported type (Parser/VCF.pm:477-481: ``<NON_REF>``, ``<CPX>``, ``<CTX>``), a span over
# ``--max_sv_size`` (Parser.pm:491-501), an incomplete deletion or an invalid span
# (``perl_skip_route``), and wrote no row for; the record's VCF line and JSON object are
# counted, not compared, as a manifest's are. It reads the input VCF in every format.
# ``breakend_mate_local_read`` (``filter_breakend_mate_local_read_pairs``): a key both engines
# write under a bracket-mate allele whose consequence sets are the local-coordinate and the
# mate-coordinate derivations over the cache's transcript model
# (BaseVariationFeatureOverlapAllele.pm:257,273 evaluate the mate row's predicates on the local
# variation feature); the pair is skipped whole, as a registry pair is, so its DISTANCE, cDNA,
# CDS and Protein positions follow. ``bnd_own_chromosome_rows_lost``
# (``filter_bnd_own_rows_lost``): a Transcript row the candidate alone writes under an own-end
# allele at a Location where the reference wrote mate rows naming another chromosome and no
# own-end row (AnnotationType/Transcript.pm:122 names the record's slice by the first
# transcript reached, :125-143 then fail every own-chromosome transcript), the transcript on
# the Location's chromosome within 5,000 bases of the breakpoint. The two breakend classes
# read the candidate's JSON cache (``--json-cache``), as the structural-variant comparator does.
BREAKEND_MATE_LOCAL_READ = "breakend_mate_local_read"
BND_OWN_CHROMOSOME_ROWS_LOST = "bnd_own_chromosome_rows_lost"
STRUCTURAL_VARIANT_CLASSES = (REFERENCE_SKIPPED_RECORD, BREAKEND_MATE_LOCAL_READ, BND_OWN_CHROMOSOME_ROWS_LOST)
CACHE_READING_CLASSES = (BREAKEND_MATE_LOCAL_READ, BND_OWN_CHROMOSOME_ROWS_LOST)
# The classes documented on the divergence pages that a manifest alone applies (a corpus's
# ``divergences``, ``vep_rs_only_tuples`` or ``field_divergences`` under the class), each
# a shape no predicate decides from the two outputs: the fields a TranscriptAnnotator table
# names carry the user's prefix; a later allele of a multi-allelic insertion record annotated
# at an earlier allele's shifted span (OutputFactory.pm:306-322, TranscriptVariationAllele.pm
# :502-527) and an indel within 1,000 bases of a sequence region's edge shifted short of the
# bases the sequence admits (TranscriptVariationAllele.pm:427-437, :462) change the notation
# of a row whose coordinates the outputs do not carry; and a matched allele's frequency read
# after the input allele's own key (BaseCacheVariation.pm:419-429 at release/116.2) decides
# which records --check_frequency keeps by a cache value the outputs do not print. The others
# are shapes of a flag run:
# under --minimal a structural record with two symbolic ALTs is trimmed as two sequences
# (Parser.pm:852, Utils/Sequence.pm:965-1038), so its rows move to a span that starts after
# the characters the two names share and ends before the ``>`` they share, marked MINIMISED,
# and a record with three or more keeps its span and loses the rows of every allele past the
# second (InputBuffer.pm:638-718, OutputFactory.pm:2311-2318): the reference's keys and the
# candidate's are both documented, the reference's with no candidate set; under
# --shift_3prime 1 the alleles of one record share one transcript variation whose span
# predicates the first allele stores at its own shifted span and whose coding coordinates
# are read back cached (BaseVariationFeatureOverlapAllele.pm:363, BaseTranscriptVariation.pm
# :142), so a later allele shifting by another length, or the unshifted reference row of
# --process_ref_homs, takes another span's terms; and the start-codon test of the same flag
# adds the shift to the unshifted cDNA start a deletion beginning in an intron lacks
# (Utils/VariationEffect.pm:977-981), so such a deletion reaching past the start codon gains
# start_retained_variant, or start_lost with a lost-start HGVSp on a transcript without a 5'
# UTR (:1069). The outputs carry neither the record's other alleles nor the unshifted span.
# Two of the classes are decided without a manifest as well, each on the rows its own comment
# below names: ``hgvs_shift_shared_across_alleles`` and ``hgvs_shift_sequence_edge_window``.
TRANSCRIPT_ANNOTATOR_SIBLING_ALLELE_ROWS = "transcript_annotator_sibling_allele_rows"
HGVS_SHIFT_SHARED_ACROSS_ALLELES = "hgvs_shift_shared_across_alleles"
# ``blosum62_hgvs_shift_shared_across_alleles``: the BLOSUM62 the reference alone writes on a
# later allele's row of ``hgvs_shift_shared_across_alleles``. Blosum62.pm:110 scores a row
# whose peptide pair is one residue over one residue, with the matrix cell of the two
# (:50-71, :73, :84-88). On a later allele's row the reference reads the reference residue at
# the span an earlier allele's notation moved, so an insertion between two codons, whose own
# pair has no reference residue and gets no score, carries one residue over one and a score.
# The shape: the candidate writes no BLOSUM62 and its own ``Amino_acids`` fails the plugin's
# test; the row's CDS and protein positions are the base class's
# (``later_allele_shared_shift``); the reference's ``Amino_acids`` is one residue over one,
# the first the translation of the moved pair's codon and the second the candidate's own
# alternate residue; and the reference's value is the matrix cell of that pair. A value that
# is not that cell, a value on the candidate's side, a candidate pair the plugin scores, a
# first allele's row, a record whose first allele carries no offset and a run without the
# cache stay differences.
BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES = "blosum62_hgvs_shift_shared_across_alleles"
BLOSUM62_FIELDS = ("BLOSUM62",)
BLOSUM62_RESIDUES = "ARNDCQEGHILKMFPSTWYV"
BLOSUM62_ROWS = (
    "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0",
    "-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3",
    "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3",
    "-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3",
    "0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1",
    "-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2",
    "-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2",
    "0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3",
    "-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3",
    "-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3",
    "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1",
    "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2",
    "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1",
    "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1",
    "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2",
    "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2",
    "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0",
    "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3",
    "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1",
    "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4",
)
BLOSUM62_MATRIX = {
    (first, second): cell
    for first, row in zip(BLOSUM62_RESIDUES, BLOSUM62_ROWS)
    for second, cell in zip(BLOSUM62_RESIDUES, row.split())
}
PEPTIDE_PAIR_RE = re.compile(r"([A-Z])/([A-Z])")
HGVS_SHIFT_SEQUENCE_EDGE_WINDOW = "hgvs_shift_sequence_edge_window"
FREQUENCY_FILTER_ALLELE_KEY_ORDER = "frequency_filter_allele_key_order"
MINIMAL_SYMBOLIC_ALLELES = "minimal_symbolic_alleles"
SHIFT_3PRIME_SHARED_PREDICATES = "shift_3prime_shared_predicates"
SHIFT_START_CODON_TEST = "shift_start_codon_test"
# ``transcript_annotator_insertion_window``: TranscriptAnnotator.pm:218 reads the table rows
# from ``start - 2`` to ``end``, a two-base window for an insertion (``end`` = ``start`` - 1)
# that never reaches a row written at the insertion's start; the candidate reads from two
# bases before the lower coordinate to the higher. The shape, the table's fields the candidate
# alone fills on an insertion's row, is documented from a corpus manifest.
TRANSCRIPT_ANNOTATOR_INSERTION_WINDOW = "transcript_annotator_insertion_window"
MANIFEST_CLASSES = (
    TRANSCRIPT_ANNOTATOR_SIBLING_ALLELE_ROWS,
    TRANSCRIPT_ANNOTATOR_INSERTION_WINDOW,
    HGVS_SHIFT_SHARED_ACROSS_ALLELES,
    HGVS_SHIFT_SEQUENCE_EDGE_WINDOW,
    FREQUENCY_FILTER_ALLELE_KEY_ORDER,
    MINIMAL_SYMBOLIC_ALLELES,
    SHIFT_3PRIME_SHARED_PREDICATES,
    SHIFT_START_CODON_TEST,
)
# The classes under which a record the reference writes nothing for is documented by its
# candidate-only tuples: the record's VCF line and JSON object are counted, not compared, as
# a ``reference_skipped_record``'s are.
RECORD_DROPPING_CLASSES = (REFERENCE_SKIPPED_RECORD, FREQUENCY_FILTER_ZERO_VALUE, FREQUENCY_FILTER_ALLELE_KEY_ORDER)
# ``hgvs_shift_shared_across_alleles`` beside the candidate JSON cache: three columns of a
# later allele's row are decided without a manifest. The alleles of a record share one
# TranscriptVariation per transcript, which keeps the CDS pair and the translation pair once
# read (BaseTranscriptVariation.pm:252-267, :365-380). The first allele's protein notation
# reads both again at its 3'-shifted position and leaves them there
# (TranscriptVariationAllele.pm:1724-1733, :844-858, :2404-2410), and reads the shared
# reference allele's peptide again at those positions (:1780-1783, :1793). A later allele's
# row prints the kept pairs (OutputFactory.pm:1737, :1741) and that peptide before its own
# (:1733, TranscriptVariationAllele.pm:618-630). Its own peptide stands at the record's own
# position: a run that keeps statistics reads every allele's consequences, and with them its
# peptide, before the record's first row (OutputFactory.pm:347, :1022; Stats.pm:241;
# VariationFeature.pm:1000-1002; Utils/VariationEffect.pm:785). The shape
# (``later_allele_shared_shift``): a differing CDS_position, Protein_position or Amino_acids
# on the row of a second or later alternate allele of its record (``later_alternate``) whose
# first allele's reference row on the feature carries a non-zero HGVS_OFFSET
# (``first_allele_shift``). The candidate prints the record's own CDS pair and the reference
# that span moved by the first allele's offset, both mapped through the cache's mapper pairs
# (``shifted_pair_readings``). Each side's Protein_position is the codons of its pair, and
# the two Amino_acids share their alternate part, each side's reference part the coding
# sequence's translation over its codons under the transcript's table. A first allele's row,
# a first allele without an offset, a move of another size, a candidate value that is not the
# record's own, any other column, a row whose record's alleles no entry states in order and a
# transcript the cache holds without its mapper pairs stay differences; without the cache the
# class decides nothing from the two outputs.
SHARED_SHIFT_FIELDS = (*CDS_POSITION_FIELDS, *PROTEIN_POSITION_FIELDS, *AMINO_ACIDS_FIELDS)
# ``hgvs_shift_sequence_edge_window`` beside a reference FASTA: ``_genomic_shift`` cuts the two
# flanks ``perform_shift`` walks as the first and the last 1,000 bases of a slice it expanded by
# 1,000 bases each side and clipped to the region (TranscriptVariationAllele.pm:427-437), so
# within 1,000 bases of an edge the flank on the clipped side is not the allele's. The region's
# length and its bases are in neither output, so the predicate reads them from the reference
# and repeats both walks (:303-349): over the allele's own flanks, which is the shift the
# sequence admits within the bound two whole flanks give, and over the clipped slice's two
# cuts. Where each walk leaves the notation depends on the transcript's span, which no column
# carries, so the predicate reads the candidate JSON cache as well. The shape: a differing
# HGVSc or HGVS_OFFSET on an insertion's or deletion's Transcript row whose strand walks
# towards an edge less than 1,000 bases away, the two walks differing, and each side's HGVSc
# and HGVS_OFFSET those its walk gives on the row's transcript. Ensembl VEP clamps the written
# span onto the transcript and adds its walk, writing no notation where the clamped end then
# passes the last base (:2729-2733, :1472, :1477-1478); the candidate moves the span by its
# walk and clamps what is left inside, writing none where nothing is. A deletion containing
# the transcript whole so keeps one notation and differs in its offset alone, and a row whose
# shift crosses the transcript's edge is read with the clamp each engine applies. An
# insertion duplicates the bases before its point where they repeat the allele as the walk
# turned it and inserts that allele where they do not, Ensembl VEP duplicating the bases
# after the point first (Utils/Sequence.pm:573-596). An offset is written beside a notation or
# an HGVSp (OutputFactory.pm:1755, :1767). A value other than those, the same values further
# from the edge or on the other strand, a coordinate inside an intron, an insertion written
# at or over the transcript's edge and a transcript the cache lacks stay differences.
SEQUENCE_EDGE_WINDOW = 1000
SEQUENCE_EDGE_FIELDS = (*HGVSC_FIELDS, *HGVS_OFFSET_FIELDS)
# The built-in field classes, each with the field names or name prefixes it reads, matched
# without regard to case (the JSON format lower-cases every key). A differing field two
# classes accept is credited to the first of them in this order.
FIELD_CLASSES = {
    GXF_REVERSE_STRAND_UTR_SEQUENCE: HGVSP_FIELDS,
    ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE: HGVSP_FIELDS,
    CONSERVATION_MAX_TEXT_ORDER: CONSERVATION_FIELDS,
    GWAS_HAPLOTYPE_CROSS_PRODUCT: (GWAS_PREFIX,),
    GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE: (GWAS_PREFIX,),
    ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION: ALPHAMISSENSE_PREFIXES,
    CADD_STRUCTURAL_ALLELE_MATCH: (CADD_PREFIX,),
    SPLICEVAULT_SIBLING_ALLELE_ROWS: SPLICEVAULT_PREFIXES,
    SPLICEVAULT_JSON_LATER_ROWS: (SPLICEVAULT_JSON_FIELD,),
    UTRANNOTATOR_ANNOTATION_PAIR_ORDER: UTRANNOTATOR_ANNOTATION_FIELDS,
    UTRANNOTATOR_GAINED_START_WINDOW: (UTRANNOTATOR_PREFIX,),
    UTRANNOTATOR_CONSEQUENCE_PRECEDENCE: UTRANNOTATOR_CONSEQUENCE_FIELDS + UTRANNOTATOR_ANNOTATION_FIELDS,
    UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT: UTRANNOTATOR_CONSEQUENCE_FIELDS + UTRANNOTATOR_ANNOTATION_FIELDS,
    UTRANNOTATOR_LOST_START_AFTER_INSERTION: UTRANNOTATOR_CONSEQUENCE_FIELDS + UTRANNOTATOR_ANNOTATION_FIELDS,
    UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE: UTRANNOTATOR_CONSEQUENCE_FIELDS + UTRANNOTATOR_ANNOTATION_FIELDS,
    UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY: UTRANNOTATOR_ANNOTATION_FIELDS,
    UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON: UTRANNOTATOR_CONSEQUENCE_FIELDS + UTRANNOTATOR_ANNOTATION_FIELDS,
    NEAREST_EXON_JB_ITEM_ORDER: NEAREST_EXON_JB_FIELDS,
    LOEUF_ZERO_MINIMUM: LOEUF_FIELDS,
    MAXENTSCAN_OVERLAPPING_WINDOWS: (MAXENTSCAN_PREFIX,),
    ALPHAMISSENSE_SIBLING_ALLELE_ROWS: ALPHAMISSENSE_PREFIXES,
    REVEL_SIBLING_ALLELE_ROWS: REVEL_FIELDS,
    REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS: REVEL_FIELDS,
    GNOMADMT_INFO_VALUE_SECOND_EQUALS: (GNOMAD_PLUGIN_PREFIX,),
    LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD: LOF_FIELDS + LOF_FILTER_FIELDS,
    LOFTEE_FIRST_ALT_DONOR_FILTER: LOF_FIELDS + LOF_FILTER_FIELDS,
    LOFTEE_GERP_ZOOM_LEVEL_SUMMARY: LOF_INFO_FIELDS + LOF_FIELDS + LOF_FILTER_FIELDS,
    NEAREST_EXON_JB_STRETCH_PRIMING: NEAREST_EXON_JB_FIELDS,
    NEAREST_EXON_JB_SHIFTED_SPAN: NEAREST_EXON_JB_FIELDS,
    NEAREST_EXON_JB_ID_LESS_EXON_FOLD: NEAREST_EXON_JB_FIELDS,
    DOWNSTREAM_REVERSE_STRAND_ROTATION: DOWNSTREAM_FIELDS,
    # The registry pair class on the fields the Downstream plugin's gate withholds:
    # Downstream.pm:95-99 write on a frameshift row none of whose terms names a splice site, so
    # on a key whose Consequence sets pair by the class's rule the reference, whose terms pass
    # the gate, alone fills them.
    SPLICE_FAMILY_SWAP: DOWNSTREAM_FIELDS,
    SPLICEAI_CHOPPED_MULTIALLELIC: SPLICEAI_PREFIXES,
    SPLICEAI_MULTI_ENTRY_CELL: SPLICEAI_PREFIXES,
    PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING: TABIX_PLUGIN_FIELDS,
    REFSEQ_EDITED_MODEL_HGVSC_OFFSET: HGVSC_FIELDS,
    REFSEQ_EDITED_MODEL_SHIFTED_INSERTION_BASES: HGVSC_FIELDS,
    REFSEQ_EDITED_MODEL_SEAM_INSERTION_HGVSC: HGVSC_FIELDS,
    REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC: HGVSC_FIELDS,
    REFSEQ_EDITED_MODEL_HGVS_OFFSET: HGVS_OFFSET_FIELDS,
    HGVS_INTRON_OFFSET_ABUTTING_INSERTION: HGVS_INTRON_OFFSET_FIELDS,
    HGVS_SHIFTED_INSERTION_CLAMPED_END: HGVSC_FIELDS,
    HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE: SHIFTED_INDEL_GATE_FIELDS,
    HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_GATE: HGVSP_FIELDS,
    HGVS_SHIFTED_DELETION_CLAMPED_SPAN: SHIFTED_DELETION_SPAN_FIELDS,
    HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE: CLIPPED_REFERENCE_FIELDS,
    HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END: ADDED_BASES_FIELDS,
    HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS: SHIFT_BOUND_FIELDS,
    HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT: ADDED_BASES_FIELDS,
    HGVS_SHIFT_SEQUENCE_EDGE_WINDOW: SEQUENCE_EDGE_FIELDS,
    HGVS_SHIFT_SHARED_ACROSS_ALLELES: SHARED_SHIFT_FIELDS,
    HGVSP_STOP_CODON_PROTEIN_COMPARISON: HGVSP_FIELDS,
    HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE: HGVSP_FIELDS,
    HGVSP_INTRON_END_INSERTION: HGVSP_FIELDS,
    HGVSP_INTRON_START_INSERTION: HGVSP_FIELDS,
    SHIFT_INTRON_TERM_UNSHIFTED_SPAN: CONSEQUENCE_FIELDS,
    REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS: CODONS_FIELDS,
    FREQUENCY_FILTER_ZERO_VALUE: FREQS_FIELDS,
    HGNC_ID_COPY: HGNC_ID_FIELDS,
    DOWNSTREAM_UTR_FRAME_PAST_CODING_END: DOWNSTREAM_PROTEIN_FIELDS,
    DOWNSTREAM_START_PAST_OWN_CODON: DOWNSTREAM_PROTEIN_FIELDS,
    PLUGIN_HGVS_SHIFTED_CDS_COORDINATES: DOWNSTREAM_PROTEIN_FIELDS + PROTEIN_LENGTH_CHANGE_FIELDS + NMD_FIELDS,
    PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD: DOSAGE_SENSITIVITY_FIELDS + UTRANNOTATOR_EXISTING_FIELDS,
    BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES: BLOSUM62_FIELDS,
}
# The field classes whose shape spans a record's rows: ``compare_entries`` groups the
# reference's entries by record (``RecordRows``) when one of them is requested. Those that
# read a record's alleles read them where an entry names one record (``RecordRows.one_record``),
# so a default or tab row under an identifier stays a difference under them unless the input
# VCF names its record: the one line that carries the identifier, or of several the one that
# writes the row at its Location with its allele (``RecordRows.line_of``); the gained-start
# class also reads such a row's own text where it states the reference allele
# (``record_allele_lengths``).
RECORD_READING_CLASSES = (
    GWAS_HAPLOTYPE_CROSS_PRODUCT,
    GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE,
    SPLICEVAULT_SIBLING_ALLELE_ROWS,
    ALPHAMISSENSE_SIBLING_ALLELE_ROWS,
    REVEL_SIBLING_ALLELE_ROWS,
    REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS,
    LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD,
    LOFTEE_FIRST_ALT_DONOR_FILTER,
    SPLICEAI_CHOPPED_MULTIALLELIC,
    PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING,
    UTRANNOTATOR_GAINED_START_WINDOW,
    HGVS_SHIFT_SHARED_ACROSS_ALLELES,
    PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
    BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES,
)
# Those of them that read a record's alleles: under one of them the comparison reads the input
# VCF in the default and tab formats too, to tell which identifier names one record.
ALLELE_READING_CLASSES = (
    SPLICEVAULT_SIBLING_ALLELE_ROWS,
    ALPHAMISSENSE_SIBLING_ALLELE_ROWS,
    REVEL_SIBLING_ALLELE_ROWS,
    REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS,
    LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD,
    LOFTEE_FIRST_ALT_DONOR_FILTER,
    SPLICEAI_CHOPPED_MULTIALLELIC,
    UTRANNOTATOR_GAINED_START_WINDOW,
    HGVS_SHIFT_SHARED_ACROSS_ALLELES,
    PLUGIN_HGVS_SHIFTED_CDS_COORDINATES,
    BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES,
)
# The classes of a value the reference copies from a sibling allele's row: they read the
# candidate's rows of the record as well, to tell the allele the table names, whose row both
# engines fill, from the sibling's.
SIBLING_ALLELE_CLASSES = (SPLICEVAULT_SIBLING_ALLELE_ROWS, ALPHAMISSENSE_SIBLING_ALLELE_ROWS, REVEL_SIBLING_ALLELE_ROWS)
# The notes of a floor a format cannot count: the JSON format's row-less columns, a custom
# source's columns in the JSON format (JSON.pm:249-258 at release/116.2 put a record's custom
# annotations, their statistics first, on the record object alone), and a column neither
# header of a tab or VCF output declares.
FLOOR_NOT_APPLICABLE = (
    "a column the JSON format has no row count for (a record-level column is carried once "
    "per record or not at all, Feature_type names the consequence list, Extra is the text "
    "formats' container), so the floor is checked in the text formats of the same case"
)
FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE = (
    "a custom source's column, which the JSON format carries once per record under "
    "custom_annotations and not on its consequence entries, so the floor is checked in the "
    "text formats of the same case"
)
FLOOR_UNDECLARED_NOT_APPLICABLE = (
    "a column neither output's header declares in this format (the tab header line and the "
    "CSQ Format list carry the fields the run's options name, where the default format's Extra "
    "carries every key a row has), so the floor is checked in the default format of the same case"
)
# The record field holding a JSON record's custom annotations, an object keyed by source.
CUSTOM_ANNOTATIONS_FIELD = "custom_annotations"
# The columns Ensembl VEP's own options write (the FLAG_FIELDS list of Constants.pm at
# release/116.2 with the default output columns), which the JSON factory keys at the top of a
# consequence entry wherever they are populated: a floor on one of them is never read inside a
# plugin's object, whose member of the same name (SpliceAI's SYMBOL) is the plugin's copy, not
# the column.
CORE_COLUMNS = {
    name.lower()
    for name in [
        "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type", "Consequence",
        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation",
        "Extra", "AF", "AFR_AF", "ALLELE_NUM", "AMBIGUITY", "AMR_AF", "APPRIS", "BAM_EDIT", "BIOTYPE",
        "CANONICAL", "CCDS", "CELL_TYPE", "CHECK_REF", "CLINVAR_SOMATIC_CLASSIFICATION", "CLIN_SIG", "DISTANCE",
        "DOMAINS", "EAS_AF", "ENSP", "EUR_AF", "EXON", "FLAGS", "FREQS", "GA4GH_VRS", "GENCODE_PRIMARY",
        "GENE_PHENO", "GIVEN_REF", "HGNC_ID", "HGVS_OFFSET", "HGVSc", "HGVSg", "HGVSp", "HIGH_INF_POS", "IMPACT",
        "IND", "INTRON", "MANE", "MANE_PLUS_CLINICAL", "MANE_SELECT", "MAX_AF", "MAX_AF_POPS", "MINIMISED",
        "MOTIF_NAME", "MOTIF_POS", "MOTIF_SCORE_CHANGE", "NEAREST", "OverlapBP", "OverlapPC", "PHENO", "PICK",
        "PUBMED", "PolyPhen", "REFSEQ_MATCH", "REFSEQ_OFFSET", "REF_ALLELE", "RefSeq", "SAS_AF", "SHIFT_LENGTH",
        "SIFT", "SOMATIC", "SOURCE", "SPDI", "STRAND", "SV", "SWISSPROT", "SYMBOL", "SYMBOL_SOURCE",
        "TRANSCRIPTION_FACTORS", "TREMBL", "TSL", "UNIPARC", "UNIPROT_ISOFORM", "UPLOADED_ALLELE", "USED_REF",
        "VARIANT_CLASS", "VAR_SYNONYMS", "ZYG", "miRNA", "gnomADe_AF", "gnomADe_AFR_AF", "gnomADe_AMR_AF",
        "gnomADe_ASJ_AF", "gnomADe_EAS_AF", "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF",
        "gnomADe_REMAINING_AF", "gnomADe_SAS_AF", "gnomADg_AF", "gnomADg_AFR_AF", "gnomADg_AMI_AF",
        "gnomADg_AMR_AF", "gnomADg_ASJ_AF", "gnomADg_EAS_AF", "gnomADg_FIN_AF", "gnomADg_MID_AF",
        "gnomADg_NFE_AF", "gnomADg_REMAINING_AF", "gnomADg_SAS_AF",
    ]
}
# Config.pm:501-508 at release/116.2: --vcf implies --symbol, --biotype and --numbers, so a
# VCF run carries these columns where a default run without those flags does not.
VCF_IMPLIED_COLUMNS = {"SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "EXON", "INTRON"}

# ensembl-variation Utils/Constants.pm %OVERLAP_CONSEQUENCES at release/116, SO term to rank.
CONSEQUENCE_RANK = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "feature_elongation": 9,
    "feature_truncation": 10,
    "inframe_insertion": 11,
    "inframe_deletion": 12,
    "missense_variant": 13,
    "protein_altering_variant": 14,
    "splice_donor_5th_base_variant": 15,
    "splice_region_variant": 16,
    "splice_donor_region_variant": 17,
    "splice_polypyrimidine_tract_variant": 18,
    "incomplete_terminal_codon_variant": 19,
    "start_retained_variant": 20,
    "stop_retained_variant": 21,
    "synonymous_variant": 22,
    "coding_sequence_variant": 23,
    "mature_miRNA_variant": 24,
    "5_prime_UTR_variant": 25,
    "3_prime_UTR_variant": 26,
    "non_coding_transcript_exon_variant": 27,
    "intron_variant": 28,
    "NMD_transcript_variant": 29,
    "non_coding_transcript_variant": 30,
    "coding_transcript_variant": 31,
    "upstream_gene_variant": 32,
    "downstream_gene_variant": 33,
    "TFBS_ablation": 34,
    "TFBS_amplification": 35,
    "TF_binding_site_variant": 36,
    "regulatory_region_ablation": 37,
    "regulatory_region_amplification": 38,
    "regulatory_region_variant": 39,
    "intergenic_variant": 40,
    "sequence_variant": 41,
}


def read_text(path: Path) -> str:
    """The text of ``path``, gunzipped when its name ends in ``.gz``."""
    try:
        if str(path).endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as fh:
                return fh.read()
        return Path(path).read_text(encoding="utf-8")
    except (OSError, EOFError, UnicodeDecodeError, zlib.error) as exc:
        raise ComparatorError(f"cannot read {path}: {getattr(exc, 'strerror', None) or exc}") from exc


def is_directory(path: Path) -> bool:
    """Whether ``path`` is a directory. A path the system cannot look up (a name longer than a
    file name can be) is a malformed input and ends the run with the comparator's error,
    whatever the interpreter's own ``is_dir`` makes of it."""
    try:
        return stat.S_ISDIR(os.stat(path).st_mode)
    except (FileNotFoundError, NotADirectoryError):
        return False
    except (OSError, ValueError) as exc:
        raise ComparatorError(f"cannot read {str(path)[:200]}: {getattr(exc, 'strerror', None) or exc}") from exc


def text_lines(text: str) -> list[str]:
    """The lines of ``text`` without their terminators, no trailing empty line."""
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return [line.removesuffix("\r") for line in lines]


def input_records(path: Path) -> list[str]:
    """The data lines of an input VCF, in order."""
    return [line for line in text_lines(read_text(path)) if not line.startswith("#")]


@dataclass
class Entry:
    """One row of the default or tab format, one CSQ entry or record line of the VCF
    format, or one consequence object or record object of the JSON format.

    ``key`` is the entry's identity within its format: the Location (default, tab) or the
    input record ordinal (VCF, JSON), the allele, the feature and the feature type, plus a
    ``#n`` counter on exact duplicates.
    """

    key: str
    record: int | None
    location: str | None
    allele: str
    feature: str
    consequence: str
    fields: dict[str, str] = field(default_factory=dict)
    # Whether the allele is an insertion's or a deletion's, read from the record's own alleles
    # (``allele_is_indel`` on the VCF columns, the JSON allele string or an identifier-less
    # Uploaded_variation); ``None`` where they say nothing, and ``indel_row`` falls back to
    # the Location.
    indel: bool | None = None
    # Whether the same alleles make every sequence allele an insertion, the reference allele
    # empty once the shared first base is dropped (``reference_allele_is_empty``); ``None``
    # where the entry's record spells no alleles.
    insertion: bool | None = None
    # The reference bases the allele stands in for at the entry's Location, read by position
    # from the record's own alleles (``replaced_reference``): an identifier-less
    # Uploaded_variation, the VCF line's columns or the JSON allele string, else the records
    # of the input VCF where a class reads them (``read_entry_references``). Empty for a pure
    # insertion's entry; ``None`` where no record was read.
    reference: str | None = None
    # The ordinal of the input line a VCF entry's line is written from: the record's, or, where
    # the line's own columns name no input record, the Ensembl-format line an engine writes as
    # a VCF line of its own (``ensembl_line_ordinal``). ``None`` for any other entry.
    line: int | None = None
    # The Location the default format writes for a VCF or JSON entry's allele, where the
    # record's own position and alleles give the allele there; ``None`` for any other entry.
    # ``at_its_location`` hands it to a predicate that reads a text row's Location.
    allele_location: str | None = None


def dedup_keys(entries: list[Entry]) -> None:
    """Suffixes the second and later occurrences of a key with ``#n``."""
    seen: Counter[str] = Counter()
    for e in entries:
        seen[e.key] += 1
        n = seen[e.key]
        if n > 1:
            e.key = f"{e.key}#{n}"


def canonical_json_text(text: str) -> str:
    """A JSON object's text with its members sorted (``json_compact``), or ``text`` as it is
    when it is not a JSON object. Ensembl VEP encodes the ``GA4GH_VRS`` object with its keys
    in hash order, which differs from row to row, so the text formats compare the parsed
    object, as the JSON format's nested object does. A text nested deeper than the
    interpreter reads or walks is a malformed input and ends the run with the comparator's
    error."""
    try:
        value = json.loads(text)
    except RecursionError as exc:
        raise ComparatorError(f"a JSON value is nested too deep to read: {text[:24]}...") from exc
    except ValueError:
        return text
    if not isinstance(value, dict):
        return text
    try:
        return json_compact(value)
    except RecursionError as exc:
        raise ComparatorError(f"a JSON object is nested too deep to read: {text[:24]}...") from exc


def numberified_json_text(text: str) -> str:
    """``canonical_json_text`` with every string member that is a whole number made a number:
    the ``GA4GH_VRS`` object of the text formats carries its ``start`` as a string, and the
    JSON format's nested object has both bounds as numbers (``OutputFactory/JSON.pm:288``,
    ``numberify``)."""
    try:
        value = json.loads(text)
    except ValueError:
        return text
    if not isinstance(value, dict):
        return text

    def numberify(v):
        if isinstance(v, dict):
            return {k: numberify(x) for k, x in v.items()}
        if isinstance(v, list):
            return [numberify(x) for x in v]
        return int(v) if isinstance(v, str) and v.isdecimal() else v

    return json_compact(numberify(value))


# The fields whose text is a JSON object, compared in canonical form.
JSON_OBJECT_FIELDS = ("GA4GH_VRS",)


def parse_default(text: str) -> list[Entry]:
    """The default format: the columns of the ``#`` header line, the ``Extra`` column's
    ``key=value`` pairs becoming fields of their own while ``Extra`` keeps its raw text, a
    ``JSON_OBJECT_FIELDS`` value in canonical form in both."""
    cols: list[str] = []
    entries: list[Entry] = []
    for line in text_lines(text):
        if line.startswith("#"):
            if not line.startswith("##"):
                cols = line[1:].split("\t")
            continue
        vals = line.split("\t")
        if len(vals) != len(cols):
            raise ComparatorError(f"column count on line: {line}")
        fields = dict(zip(cols, vals))
        for name in JSON_OBJECT_FIELDS:
            if name in fields:
                fields[name] = canonical_json_text(fields[name])
        extra = fields.get("Extra")
        if extra is not None and extra != "-":
            pairs: list[str] = []
            for kv in extra.split(";"):
                k, sep, v = kv.partition("=")
                if sep and k in JSON_OBJECT_FIELDS:
                    v = canonical_json_text(v)
                    kv = f"{k}={v}"
                fields[k] = v
                pairs.append(kv)
            fields["Extra"] = ";".join(pairs)
        location = fields.get("Location", "")
        allele = fields.get("Allele", "")
        feature = fields.get("Feature", "")
        indel = insertion = reference = None
        uploaded = fields.get("Uploaded_variation", "")
        record_alleles = uploaded_alleles(uploaded)
        if record_alleles is not None:
            reference_allele, alternate_alleles = record_alleles
            indel = allele_is_indel(reference_allele, alternate_alleles, allele)
            insertion = reference_allele_is_empty(reference_allele, alternate_alleles)
            position = uploaded_position(uploaded)
            if position is not None:
                reference = replaced_reference(position, reference_allele, alternate_alleles, location, allele)
        entries.append(
            Entry(
                key="|".join([location, allele, feature, fields.get("Feature_type", "")]),
                record=None,
                location=location,
                allele=allele,
                feature=feature,
                consequence=fields.get("Consequence", ""),
                fields=fields,
                indel=indel,
                insertion=insertion,
                reference=reference,
            )
        )
    dedup_keys(entries)
    return entries


def parse_tab(text: str) -> list[Entry]:
    """The tab format: every column named in the ``#`` header line."""
    return parse_default(text)


class InputIndex:
    """Input records by identity, for mapping an output record back to its ordinal.
    Identical input lines are matched in order of appearance."""

    def __init__(self, inputs: list[str], key: Callable[[str], str]):
        self.unmatched: dict[str, deque[int]] = defaultdict(deque)
        for i, line in enumerate(inputs):
            self.unmatched[key(line)].append(i)

    def record(self, key: str) -> int | None:
        queue = self.unmatched.get(key)
        return queue.popleft() if queue else None


def vcf_line_key(line: str) -> str:
    """A VCF data line's identity: CHROM, POS, ID, REF and ALT."""
    return "\t".join(line.split("\t")[:5])


# The strand an Ensembl-format input line's fifth field names.
ENSEMBL_FORMAT_STRANDS = {"+": 1, "1": 1, "+1": 1, "-": -1, "-1": -1}


def ensembl_format_record(line: str) -> tuple[str, int, str, list[str]] | None:
    """The chromosome, first position, reference allele and alternate alleles of an input line
    in the Ensembl default format: the chromosome, an integer start and end, the allele string
    ``REF/ALT[/ALT...]`` with ``-`` for an empty allele, and the strand, separated by white
    space. The alleles come in upper case on the forward strand, those of a line on the reverse
    strand reverse-complemented. An insertion written alone has its start one above its end and
    an empty reference allele. ``None`` for a line of any other shape."""
    fields = line.split()
    if len(fields) < 5 or not (fields[1].isdigit() and fields[2].isdigit()):
        return None
    if "/" not in fields[3] or fields[4] not in ENSEMBL_FORMAT_STRANDS:
        return None
    alleles = ["" if allele == "-" else allele.upper() for allele in fields[3].split("/")]
    if ENSEMBL_FORMAT_STRANDS[fields[4]] < 0:
        alleles = [reverse_complement(allele) for allele in alleles]
    return fields[0], int(fields[1]), alleles[0], alleles[1:]


def input_line_record(line: str) -> tuple[str, int, str, list[str]] | None:
    """The chromosome, first position, reference allele and alternate alleles an input line
    gives, the alleles in upper case on the forward strand. The line's fourth field tells its
    format: a VCF line's REF is letters alone, beside an integer POS and an ALT column, and an
    Ensembl-format line's allele string holds a ``/`` (``ensembl_format_record``). ``None``
    for a line of neither shape."""
    columns = line.split("\t")
    if len(columns) >= 5 and columns[1].isdigit() and re.fullmatch(r"[A-Za-z]+", columns[3]) is not None:
        return columns[0], int(columns[1]), columns[3].upper(), columns[4].upper().split(",")
    return ensembl_format_record(line)


def ensembl_line_ordinal(inputs: list[str], ordinal: int, chromosome: str, position: str) -> int | None:
    """The ordinal of the Ensembl-format input line a VCF output line is written from, for a
    line whose own columns name no input record: an engine writes one VCF line per input line,
    in order, at the input line's start or, where it pads the alleles with the base before
    them, one position below. ``ordinal`` is the output line's own, ``chromosome`` and
    ``position`` its CHROM and POS; ``None`` where the input line there is no such line."""
    record = ensembl_format_record(inputs[ordinal]) if ordinal < len(inputs) else None
    if record is None or record[0] != chromosome or not position.isdigit():
        return None
    return ordinal if int(position) in (record[1], record[1] - 1) else None


def csq_format_fields(line: str) -> list[str]:
    """The CSQ field names of a ``##INFO=<ID=CSQ`` header line's ``Format:`` list."""
    parts = line.split("Format: ", 1)
    if len(parts) < 2:
        raise ComparatorError(f"CSQ header without a Format list: {line}")
    fmt = parts[1]
    while fmt.endswith('">'):
        fmt = fmt[:-2]
    return fmt.split("|")


def header_columns(text: str) -> list[str]:
    """The columns of a default or tab output's ``#`` header line, empty without one."""
    for line in text_lines(text):
        if line.startswith("#") and not line.startswith("##"):
            return line[1:].split("\t")
    return []


def vcf_line_span(chromosome: str, position: str, reference: str, alternates: list[str]) -> tuple[int, str] | None:
    """The POS of a VCF line and the Location the default format writes for its rows: the
    line's span, without the first base where Parser/VCF.pm:300-336 drop it
    (``vcf_alleles_as_spelt``), an emptied reference allele's two neighbouring positions.
    ``None`` for a POS that is not a run of digits."""
    if re.fullmatch(r"[0-9]+", position) is None:
        return None
    first = integer_of(position, "VCF line POS")
    dropped = vcf_alleles_as_spelt(reference, alternates)[0] != reference
    return first, sequence_location(chromosome, first + dropped, first + len(reference) - 1)


def parse_vcf(text: str, inputs: list[str]) -> list[Entry]:
    """The VCF format: one entry per CSQ entry, keyed by the input record ordinal, the
    allele, the feature and the feature type, plus one record entry per line carrying the
    seven fixed columns, the INFO without CSQ and whether a CSQ was present. A line
    matching no input record is keyed by its output ordinal; its CSQ entries carry the ordinal
    of the Ensembl-format input line it is written from, where it is (``Entry.line``). A CSQ
    entry carries the reference bases its allele stands in for and the Location of the
    default format (``Entry.reference``, ``Entry.allele_location``) where the line's POS, REF
    and ALT give the allele at the line's span."""
    csq_fields: list[str] = []
    entries: list[Entry] = []
    index = InputIndex(inputs, vcf_line_key)
    lines_seen = 0
    for line in text_lines(text):
        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = csq_format_fields(line)
            continue
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 8:
            raise ComparatorError(f"VCF line has fewer than 8 columns: {line}")
        record = index.record(vcf_line_key(line))
        rec_label = f"rec{record}" if record is not None else f"line{lines_seen}"
        input_line = record if record is not None else ensembl_line_ordinal(inputs, lines_seen, cols[0], cols[1])
        lines_seen += 1
        info_rest: list[str] = []
        csq: str | None = None
        for kv in cols[7].split(";"):
            if kv.startswith("CSQ="):
                csq = kv[4:]
            else:
                info_rest.append(kv)
        fields = dict(zip(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"], cols[:7]))
        fields["INFO_without_CSQ"] = ";".join(info_rest)
        fields["CSQ_present"] = "yes" if csq is not None else "no"
        entries.append(
            Entry(
                key=f"{rec_label}|record",
                record=record,
                location=None,
                allele="",
                feature="",
                consequence="",
                fields=fields,
            )
        )
        if csq is None:
            continue
        reference_allele, alternate_alleles = fields["REF"], fields["ALT"].split(",")
        span = vcf_line_span(fields["CHROM"], fields["POS"], reference_allele, alternate_alleles)
        for chunk in csq.split(","):
            vals = chunk.split("|")
            if len(vals) != len(csq_fields):
                raise ComparatorError(f"CSQ chunk width on {rec_label}: {chunk}")
            cfields = dict(zip(csq_fields, vals))
            allele = cfields.get("Allele", "")
            feature = cfields.get("Feature", "")
            reference = None
            if span is not None:
                reference = replaced_reference(span[0], reference_allele, alternate_alleles, span[1], allele)
            entries.append(
                Entry(
                    key=f"{rec_label}|{allele}|{feature}|{cfields.get('Feature_type', '')}",
                    record=record,
                    location=None,
                    allele=allele,
                    feature=feature,
                    consequence=cfields.get("Consequence", ""),
                    fields=cfields,
                    indel=allele_is_indel(reference_allele, alternate_alleles, allele),
                    insertion=reference_allele_is_empty(reference_allele, alternate_alleles),
                    reference=reference,
                    line=input_line,
                    allele_location=span[1] if span is not None and reference is not None else None,
                )
            )
    dedup_keys(entries)
    return entries


def _json_float(f: float) -> str:
    """A float as a JSON text, shortest round trip, the exponent without sign padding."""
    if math.isnan(f) or math.isinf(f):
        return "null"
    s = repr(f)
    if "e" in s:
        mantissa, exponent = s.split("e")
        return f"{mantissa}e{int(exponent)}"
    return s


def json_compact(v) -> str:
    """A JSON value as compact text with object keys sorted."""
    if v is None or isinstance(v, bool):
        return json.dumps(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return _json_float(v)
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, list):
        return "[" + ",".join(json_compact(x) for x in v) + "]"
    members = ",".join(f"{json.dumps(k, ensure_ascii=False)}:{json_compact(x)}" for k, x in sorted(v.items()))
    return "{" + members + "}"


def json_scalar(v) -> str:
    """A JSON value as one comparable string: a string itself, an array its items joined
    with ``,``, anything else its compact JSON text."""
    if isinstance(v, str):
        return v
    if isinstance(v, list):
        return ",".join(json_scalar(x) for x in v)
    return json_compact(v)


def json_objects(text: str) -> list[tuple[int, dict]]:
    """The (line ordinal, object) of every non-empty line of a JSON Lines output."""
    out = []
    for n, line in enumerate(text_lines(text)):
        if not line.strip():
            continue
        try:
            obj = json.loads(line)
        except ValueError as exc:
            raise ComparatorError(f"JSON line {n + 1} does not parse: {exc}") from exc
        if not isinstance(obj, dict):
            raise ComparatorError(f"JSON line {n + 1} is not an object")
        out.append((n, obj))
    return out


def json_record_label(index: InputIndex, n: int, obj: dict) -> tuple[int | None, str]:
    """The input record ordinal of a JSON object and its label. An object names its record
    by the first five columns of its ``input``, as a VCF line does: Ensembl VEP stores a
    structural variant's line with an empty or ``0`` QUAL, FILTER or INFO written as ``.``
    (``Parser/VCF.pm:464-466``), so the whole ``input`` may differ from the line as read."""
    inp = obj.get("input")
    record = index.record(vcf_line_key(inp) if isinstance(inp, str) else "")
    return record, (f"rec{record}" if record is not None else f"obj{n}")


def json_object_span(obj: dict) -> tuple[int, str] | None:
    """A JSON object's ``start`` and the Location the default format writes for its rows
    (``record_location``); ``None`` where ``start`` or ``end`` is not a whole number."""
    bounds = obj.get("start"), obj.get("end")
    if not all(isinstance(bound, int) and not isinstance(bound, bool) for bound in bounds):
        return None
    return bounds[0], record_location(obj)


def parse_json(text: str, inputs: list[str]) -> list[Entry]:
    """The JSON format: one record entry per object (the top-level scalars other than the
    engine identity, keyed by the input record the ``input`` line names) and one entry per
    element of each consequence list, keyed by the record, ``variant_allele``, the feature id
    and the list's feature type. An object whose ``input`` matches no record is keyed by its
    output ordinal. A consequence entry carries the reference bases its allele stands in for
    and the Location of the default format (``Entry.reference``, ``Entry.allele_location``)
    where the object's ``start``, ``end`` and ``allele_string`` give the allele at that
    span."""
    index = InputIndex(inputs, vcf_line_key)
    entries: list[Entry] = []
    set_aside = {name for name, _ in JSON_CONSEQUENCE_LISTS} | {"input", JSON_IDENTITY_KEY}
    for n, obj in json_objects(text):
        record, rec_label = json_record_label(index, n, obj)
        fields = {k: json_scalar(v) for k, v in obj.items() if k not in set_aside}
        entries.append(
            Entry(
                key=f"{rec_label}|record",
                record=record,
                location=None,
                allele="",
                feature="",
                consequence="",
                fields=fields,
            )
        )
        reference_allele, *alternate_alleles = json_scalar(obj.get("allele_string", "")).split("/")
        span = json_object_span(obj)
        # The allele string is spelt without the base before ``start``, on which an insertion's
        # Location begins: ``N``, which stands for any base, takes that base's place in the
        # alleles ``replaced_reference`` reads, an allele of ``-`` being empty.
        anchored = ["N" + ("" if spelt == "-" else spelt) for spelt in (reference_allele, *alternate_alleles)]
        for list_name, kind in JSON_CONSEQUENCE_LISTS:
            items = obj.get(list_name)
            if not isinstance(items, list):
                continue
            for item in items:
                if not isinstance(item, dict):
                    raise ComparatorError(f"{rec_label}: a {list_name} element is not an object")
                cfields = {k: json_scalar(v) for k, v in item.items()}
                allele = cfields.get("variant_allele", "")
                if kind == "intergenic":
                    feature = ""
                else:
                    feature = cfields.get(f"{kind}_id", "")
                reference = None
                if span is not None:
                    reference = replaced_reference(span[0] - 1, anchored[0], anchored[1:], span[1], allele)
                entries.append(
                    Entry(
                        key=f"{rec_label}|{allele}|{feature}|{kind}",
                        record=record,
                        location=None,
                        allele=allele,
                        feature=feature,
                        consequence=cfields.get("consequence_terms", ""),
                        fields=cfields,
                        indel=allele_is_indel(reference_allele, alternate_alleles, allele, chopped=True),
                        insertion=reference_allele_is_empty(reference_allele, alternate_alleles, chopped=True),
                        reference=reference,
                        allele_location=span[1] if span is not None and reference is not None else None,
                    )
                )
    dedup_keys(entries)
    return entries


def recoder_value(key: str, value, documented: Documented) -> str | None:
    """One allele key's value in comparable form: a list as the sorted set of its items'
    canonical texts, a ``var_synonyms`` entry with its identifiers sorted and a ``mane_select``
    triple as its canonical JSON; a scalar as its canonical text. An item an applied class
    masks is dropped and counted under the class: under a ``RECODER_MASKS`` class, a string
    under one of ``RECODER_NOTATION_KEYS`` starting with one of its prefixes, or a ``mane_select``
    triple whose ``hgvsc`` or ``hgvsp`` does; under ``archive_dbsnp_synonyms``, a ``var_synonyms``
    entry of the ``Archive dbSNP`` source. A list the classes emptied is ``None``, the key
    compared as absent, since the other side writes no key for an empty list."""

    def masked_text(text: str) -> bool:
        for klass, prefixes in RECODER_MASKS.items():
            if klass in documented.masks and text.startswith(prefixes):
                documented.applied[klass] += 1
                return True
        return False

    if not isinstance(value, list):
        return json_compact(value)
    items: set[str] = set()
    for item in value:
        if key in RECODER_NOTATION_KEYS and isinstance(item, str) and masked_text(item):
            continue
        if key == "mane_select" and isinstance(item, dict):
            named = [item.get(k) for k in RECODER_TRANSCRIPT_KEYS if isinstance(item.get(k), str)]
            if any(masked_text(text) for text in named):
                continue
        if key == "var_synonyms" and isinstance(item, str) and ": " in item:
            if ARCHIVE_DBSNP_SYNONYMS in documented.masks and item.startswith(ARCHIVE_DBSNP_SOURCE):
                documented.applied[ARCHIVE_DBSNP_SYNONYMS] += 1
                continue
            source, ids = item.split(": ", 1)
            item = f"{source}: {','.join(sorted(ids.split(',')))}"
        items.add(json_compact(item))
    if value and not items:
        return None
    return "[" + ",".join(sorted(items)) + "]"


def recoder_identifier_indel(label: str, allele: str, fields: dict) -> bool:
    """Whether an allele object is an insertion or deletion of an identifier input: the input
    has no notation separator and no space, and the allele is ``-``, or an ``hgvsg`` names a
    deletion, duplication or insertion, or an ``spdi`` deletes and inserts sequences of
    different lengths, or a ``vcf_string`` has a reference and alternate of different lengths."""
    if ":" in label or any(c.isspace() for c in label):
        return False
    if allele == "-":
        return True
    hgvsg = fields.get("hgvsg")
    if isinstance(hgvsg, list) and any(isinstance(h, str) and INDEL_NOTATION.search(h) for h in hgvsg):
        return True

    def unbalanced(values, separator: str) -> bool:
        """Whether a list holds a text whose third and fourth ``separator``-joined fields, the
        deleted and inserted sequences of an SPDI or the reference and alternate of a VCF form,
        differ in length; a value that is not a list holds none."""
        if not isinstance(values, list):
            return False
        for text in values:
            parts = text.split(separator) if isinstance(text, str) else []
            if len(parts) == 4 and len(parts[2]) != len(parts[3]):
                return True
        return False

    return unbalanced(fields.get("spdi"), ":") or unbalanced(fields.get("vcf_string"), "-")


def parse_recoder(text: str, documented: Documented) -> list[Entry]:
    """The recoder format: one entry per allele of each object of the JSON array, in allele
    order, keyed by the object's ``input`` and the allele, plus one ``|warnings`` entry for an
    object's warnings; a ``null`` object is one entry keyed ``null``, an object of warnings alone
    one keyed ``warnings``. The keys ``identifier_indel_transcript_notations`` drops are left out
    first; every other value is read through ``recoder_value``."""
    try:
        objects = json.loads(text)
    except ValueError as exc:
        raise ComparatorError(f"the recoder output does not parse: {exc}") from exc
    if not isinstance(objects, list):
        raise ComparatorError("the recoder output is not a JSON array")
    entries: list[Entry] = []

    def entry(key: str, allele: str, fields: dict[str, str]) -> None:
        entries.append(Entry(key=key, record=None, location=None, allele=allele, feature="", consequence="", fields=fields))

    def comparable_fields(fields: dict) -> dict[str, str]:
        """The fields read through ``recoder_value``, a key it returns ``None`` for left out."""
        values = {k: recoder_value(k, v, documented) for k, v in fields.items()}
        return {k: v for k, v in values.items() if v is not None}

    for n, obj in enumerate(objects):
        if obj is None:
            entry("null", "", {})
            continue
        if not isinstance(obj, dict):
            raise ComparatorError(f"recoder object {n + 1} is neither an object nor null")
        alleles = {k: v for k, v in obj.items() if k != "warnings"}
        inputs = {v.get("input") for v in alleles.values() if isinstance(v, dict) and isinstance(v.get("input"), str)}
        if len(inputs) > 1:
            raise ComparatorError(f"recoder object {n + 1} names more than one input: {sorted(inputs)}")
        label = next(iter(inputs)) if inputs else "warnings"
        # An object's allele keys come out of a Perl hash, in no order.
        for allele, fields in sorted(alleles.items()):
            if not isinstance(fields, dict):
                raise ComparatorError(f"recoder object {n + 1}: allele {allele!r} is not an object")
            dropped: set[str] = set()
            if IDENTIFIER_INDEL_NOTATIONS in documented.masks and recoder_identifier_indel(label, allele, fields):
                dropped = {k for k in IDENTIFIER_INDEL_KEYS if k in fields}
                documented.applied[IDENTIFIER_INDEL_NOTATIONS] += len(dropped)
            kept = {k: v for k, v in fields.items() if k not in dropped}
            entry(f"{label}|{allele}", allele, comparable_fields(kept))
        if "warnings" in obj:
            entry(f"{label}|warnings", "", comparable_fields({"warnings": obj["warnings"]}))
    dedup_keys(entries)
    return entries


def parquet_as_tab(path: Path, canonical: bool = True) -> str:
    """A Parquet output (a file or a directory of part files), or with ``canonical`` a
    ``--tab`` output too, read through the adapter as tab text. The canonical form is the
    header line, then the rows sorted, every DOUBLE cell as its shortest decimal: a Parquet
    file stores a frequency as a number, so the ``0.0010`` and ``1.47e-05`` spellings of a
    ``--tab`` output come back as ``0.001`` and ``0.0000147``, and reading the tab side
    through the same form compares the two as numbers and every other cell verbatim. Without
    ``canonical`` the rows come in file order, which self-parity needs to pair the rows of a
    record that share a key with the default format's rows in their order."""
    options = ["--canonical"] if canonical else []
    command = [sys.executable, str(ADAPTER), *options, str(path)]
    proc = subprocess.run(command, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise ComparatorError(f"the Parquet adapter failed on {path}: {' '.join(proc.stderr.split())}")
    return proc.stdout


def header_key(line: str) -> str | None:
    """The key a header line describes: ``## <key> : <text>`` in the default and tab formats
    (OutputFactory.pm:2435-2503), ``##INFO=<ID=<key>,...>`` in the VCF; ``None`` for any
    other line."""
    if line.startswith("## "):
        key, sep, _ = line[3:].partition(" : ")
        return key if sep else None
    if line.startswith("##INFO=<ID="):
        return line[len("##INFO=<ID=") :].split(",", 1)[0]
    return None


def reference_text_masked(text: str, documented: Documented) -> str:
    """The reference output with the text classes of ``documented.masks`` applied before it is
    parsed: under ``conservation_short_name_spelling`` every ``phlyoP<N>way`` key spelt
    ``phyloP<N>way``, under ``pc_description_text`` every header line describing a ``<key>_PC``
    field with one of the two ``PC_DESCRIPTIONS`` carrying the candidate's text. A rewrite is
    counted under its class (a key per occurrence, a header line per line); the entries read
    no header line, so the second class changes what a header comparison would read and
    nothing the entries compare."""
    if CONSERVATION_SHORT_NAME_SPELLING in documented.masks:
        text, count = MISSPELT_CONSERVATION_KEY.subn(CONSERVATION_KEY_SPELLING, text)
        if count:
            documented.applied[CONSERVATION_SHORT_NAME_SPELLING] += count
    if PC_DESCRIPTION_TEXT in documented.masks:
        lines = text.split("\n")
        for i, line in enumerate(lines):
            key = header_key(line)
            if key is None or not key.endswith(PC_FIELD_SUFFIX):
                continue
            for perl_text, candidate_text in PC_DESCRIPTIONS.items():
                if perl_text in line:
                    lines[i] = line.replace(perl_text, candidate_text, 1)
                    documented.applied[PC_DESCRIPTION_TEXT] += 1
                    break
        text = "\n".join(lines)
    return text


def parse_format(
    fmt: str, path: Path, inputs: list[str] | None, rewrite: Callable[[str], str] | None = None
) -> list[Entry]:
    """Parses one output file in ``fmt``; VCF and JSON need the input records. ``parquet``
    reads a Parquet output, and ``tab-canonical`` a ``--tab`` output in the canonical form a
    Parquet output is compared against. ``rewrite`` is applied to the text before it is parsed
    (the reference side's text classes, ``reference_text_masked``), to the adapter's tab text
    for the two Parquet forms."""
    if fmt in ("parquet", "tab-canonical"):
        text = parquet_as_tab(path)
    else:
        text = read_text(path)
    if rewrite is not None:
        text = rewrite(text)
    if fmt == "default":
        return parse_default(text)
    if fmt in ("tab", "parquet", "tab-canonical"):
        return parse_tab(text)
    if inputs is None:
        raise ComparatorError(f"--format {fmt} needs --input, the input VCF")
    if fmt == "vcf":
        return parse_vcf(text, inputs)
    return parse_json(text, inputs)


def decode_extra(s: str) -> str:
    """Undoes the percent-encoding of ``=`` and ``;`` in Extra and CSQ values."""
    return s.replace("%3D", "=").replace("%3B", ";")


def documented_form(s: str) -> str:
    """A field value on the default format's spelling: ``-`` for an empty or absent
    value, percent-encoding undone, and the ``&`` the VCF format writes for every comma
    of a value (``OutputFactory/VCF.pm:401``), between the items of a list and inside a
    name alike, read as the default format's ``,``."""
    return "-" if s == ABSENT else dash(decode_extra(s).replace("&", ","))


def _record_indices(item: dict) -> list[int]:
    out = []
    for r in item.get("record_indices") or []:
        if isinstance(r, int) and not isinstance(r, bool) and r >= 0:
            out.append(r)
    return out


def _s(item: dict, key: str) -> str:
    v = item.get(key)
    return v if isinstance(v, str) else ""


def with_reference_sets(item: dict) -> list[str] | None:
    """The consequence sets a ``divergences`` item states for a run that has a reference
    sequence (``WITH_REFERENCE_SETS``), or ``None`` for an item without the key. The key is a
    non-empty list of distinct sets in sorted order, each spelt as ``consequence_set`` spells
    it and other than the item's ``vep_rs_consequence_sets`` as a whole; a list that holds
    Ensembl VEP's set (``vep_consequence_set``) holds nothing else, and states that the two
    engines agree on the key in such a run. A key of any other shape is refused, whatever the
    comparison's setting."""
    if WITH_REFERENCE_SETS not in item:
        return None
    sets = item[WITH_REFERENCE_SETS]
    reason = None
    if not isinstance(sets, list) or not sets:
        reason = "is not a non-empty list"
    elif not all(isinstance(s, str) and s and consequence_set(s) == s for s in sets):
        reason = "holds an item that is not a sorted, comma-joined set of terms"
    elif sets != sorted(set(sets)):
        reason = "is not a sorted list of distinct sets"
    elif sets == item.get("vep_rs_consequence_sets"):
        reason = "repeats vep_rs_consequence_sets and states no second outcome"
    elif consequence_set(_s(item, "vep_consequence_set")) in sets and len(sets) > 1:
        reason = "names Ensembl VEP's set beside another"
    if reason:
        key = " ".join(_s(item, name) for name in ("location", "allele", "feature"))
        raise ComparatorError(f"the manifest's item on {key}: {WITH_REFERENCE_SETS} {reason}")
    return sets


def manifest_run(manifest: dict, run: str, input_name: str | None) -> tuple[dict, str]:
    """The object whose three lists document the run ``run`` of the manifest, and the name of
    the run's input: the manifest itself for ``OWN_RUN``; for a flag set its entry where the
    entry carries ``divergences`` of its own, else the manifest, as the golden harness reads
    them. A name the manifest gives no run, and a compared input (``input_name``) other than
    the run's, are refused: a pair is of one run."""
    if run == OWN_RUN:
        lists, run_input = manifest, CORPUS_INPUT
    else:
        flag_sets = manifest.get("flag_sets")
        entry = next(
            (fs for fs in flag_sets if isinstance(fs, dict) and fs.get("name") == run), None
        ) if isinstance(flag_sets, list) else None
        if entry is None:
            raise ComparatorError(f"--run {run}: the manifest has no flag set of that name ({OWN_RUN} names its own run)")
        lists = entry if "divergences" in entry else manifest
        run_input = entry.get("input", CORPUS_INPUT)
    if input_name is not None and input_name != run_input:
        raise ComparatorError(f"--run {run} annotates {run_input}, and --input names {input_name}")
    return lists, run_input


T = TypeVar("T")


class Documented:
    """The keys a corpus manifest documents as divergent, restricted to the classes applied:
    consequence sets per key, keys only the candidate emits, field values per key, and the
    records the reference skipped. Keys are addressable by Location or by record ordinal.
    ``masks`` is the set of classes requested, the built-in ones included.
    ``reference_sequence`` is the comparison's setting where the caller states it, and
    ``run`` the manifest's run whose lists alone were read where the caller names one.
    ``assembly`` is the assembly the manifest names for its runs, ``None`` without one."""

    def __init__(self) -> None:
        self.masks: set[str] = set()
        self.run: str | None = None
        self.assembly: str | None = None
        self.reference_sequence: str | None = None
        # The keys an item states as equal on both engines in a run with a reference sequence,
        # filled under ``REFERENCE_PRESENT`` alone: such a key is compared like any other and
        # no class accepts a difference of its Consequence.
        self.equal_by_location: dict[tuple[str, str, str], bool] = {}
        self.equal_by_record: dict[tuple[int, str, str], bool] = {}
        # Per class, the manifest items counted in ``applied`` as they were read.
        self.read: Counter[str] = Counter()
        # Whether the lists read number the records of another input than the compared one: an
        # input no run of the manifest names reads the manifest's own lists, whose ordinals
        # count the records of ``CORPUS_INPUT``. No key is then looked up by ordinal.
        self.foreign_ordinals = False
        self.by_location: dict[tuple[str, str, str], list[str]] = defaultdict(list)
        self.by_record: dict[tuple[int, str, str], list[str]] = defaultdict(list)
        self.extra_by_location: set[tuple[str, str, str]] = set()
        self.extra_by_record: set[tuple[int, str, str]] = set()
        # A field divergence's key carries the row's consequence set when the manifest states
        # one ("" otherwise), as the golden harness keys it: two records sharing a Location and
        # Allele (an insertion and a substitution of the same bases) each keep their own rows.
        self.fields_by_location: dict[tuple[str, str, str, str], dict[str, str]] = defaultdict(dict)
        self.fields_by_record: dict[tuple[int, str, str, str], dict[str, str]] = defaultdict(dict)
        self.skipped_records: set[int] = set()
        self.available: Counter[str] = Counter()
        self.applied: Counter[str] = Counter()
        # The classes a value set aside under another class needs named beside it (the second
        # class of a value two causes change): used, though nothing is counted under them.
        self.needed: set[str] = set()

    @classmethod
    def from_manifest(
        cls,
        manifest: dict | None,
        masks: set[str],
        input_name: str | None = None,
        run: str | None = None,
        reference_sequence: str | None = None,
    ) -> Documented:
        """``input_name`` is the compared input's file name when the caller knows it; the runs
        documented over another input are then left out, since their record ordinals count
        that input's records. ``run`` names the one run of the manifest the pair belongs to
        (``manifest_run``): its lists alone are read, so an item of another run over the same
        input neither is expected on this pair's rows nor sets one aside; without it every
        run over the input is read as one list. ``reference_sequence`` states whether the
        candidate's run had a reference sequence: under ``REFERENCE_PRESENT`` an item that
        carries ``WITH_REFERENCE_SETS`` documents those sets, and one whose sets there are
        Ensembl VEP's own documents nothing, is not counted as applied and leaves its key to
        ``equal_with_reference``; under ``REFERENCE_ABSENT``, and where the caller states no
        setting, every item documents its ``vep_rs_consequence_sets``."""
        d = cls()
        d.masks = set(masks)
        d.reference_sequence = reference_sequence
        if run is not None and not manifest:
            raise ComparatorError(f"--run {run} names a run of a manifest, and none is given")
        if not manifest:
            return d
        assembly = manifest.get("assembly")
        d.assembly = assembly if isinstance(assembly, str) else None
        if run is not None:
            _, input_name = manifest_run(manifest, run, input_name)
            d.run = run
        elif input_name is not None and input_name != CORPUS_INPUT:
            flag_sets = manifest.get("flag_sets")
            d.foreign_ordinals = not any(
                isinstance(fs, dict) and fs.get("input", CORPUS_INPUT) == input_name
                for fs in (flag_sets if isinstance(flag_sets, list) else [])
            )

        def ordinals(item: dict) -> list[int]:
            return [] if d.foreign_ordinals else _record_indices(item)

        for _, item in d._applied_items(manifest, "divergences", masks, input_name):
            with_reference = with_reference_sets(item)
            if reference_sequence == REFERENCE_PRESENT and with_reference is not None:
                sets = with_reference
            else:
                sets = [x for x in item.get("vep_rs_consequence_sets") or [] if isinstance(x, str)]
            loc, allele, feature = _s(item, "location"), _s(item, "allele"), _s(item, "feature")
            d.by_location[(loc, allele, feature)].extend(sets)
            for r in ordinals(item):
                d.by_record[(r, allele, feature)].extend(sets)
        for klass, item in d._applied_items(manifest, "vep_rs_only_tuples", masks, input_name):
            loc, allele, feature = _s(item, "location"), _s(item, "allele"), _s(item, "feature")
            d.extra_by_location.add((loc, allele, feature))
            for r in ordinals(item):
                d.extra_by_record.add((r, allele, feature))
                if klass in RECORD_DROPPING_CLASSES:
                    d.skipped_records.add(r)
        for _, item in d._applied_items(manifest, "field_divergences", masks, input_name):
            fname, value = _s(item, "field"), _s(item, "vep_rs")
            loc, allele, feature = _s(item, "location"), _s(item, "allele"), _s(item, "feature")
            terms = consequence_set(_s(item, "consequence_set"))
            if loc:
                d._document_field(d.fields_by_location, (loc, allele, feature, terms), fname, value)
            for r in ordinals(item):
                # Without the input named, the runs over every input are read and an ordinal
                # names a different record in each; the text formats look nothing up by ordinal.
                d._document_field(d.fields_by_record, (r, allele, feature, terms), fname, value, refuse=input_name is not None)
        return d

    @staticmethod
    def _document_field(table: dict, key: tuple, fname: str, value: str, refuse: bool = True) -> None:
        """Records one documented field value; with ``refuse``, two runs over one input
        documenting the same key's field differently are refused, since a comparison that does
        not name its flag set cannot choose between them."""
        recorded = table[key].get(fname)
        if refuse and recorded is not None and recorded != value:
            raise ComparatorError(
                f"the manifest documents field {fname} on {key} as both {recorded!r} and {value!r}"
            )
        table[key][fname] = value

    @staticmethod
    def _documented_keys(item: dict, by_ordinal: bool) -> Iterator[tuple]:
        """The keys a manifest entry's value is looked up by: its Location and, with
        ``by_ordinal``, each of its record ordinals, each with the Allele, the Feature and the
        field, lower-cased as the JSON format spells it and empty where the entry names none."""
        rest = (_s(item, "allele"), _s(item, "feature"), _s(item, "field").lower())
        if _s(item, "location"):
            yield (_s(item, "location"), *rest)
        if by_ordinal:
            for r in _record_indices(item):
                yield (r, *rest)

    def _applied_items(
        self, manifest: dict, section: str, masks: set[str], input_name: str | None = None
    ) -> Iterator[tuple[str, dict]]:
        """The (class, item) pairs of one section of the manifest's own run and of each of its
        ``flag_sets`` entries whose class ``masks`` names, an item documented by several runs
        read once; every item's class is counted in ``available``, each yielded one's in
        ``applied``. With ``input_name`` known, the manifest's own run counts as a run over
        ``CORPUS_INPUT`` and a flag set over another input is skipped; a selected flag set
        carrying no ``divergences`` of its own, and an input no flag set names, inherit the
        manifest's lists, as the golden harness reads them. An item naming a second class
        (``SECOND_EXPECTED_DIVERGENCE``), a value two causes change, is yielded only where
        ``masks`` names both. Where it does not, no item of the runs read that a lookup
        would find for the same value is yielded either: one sharing a key of
        ``_documented_keys`` with it (an ordinal names one record only where the compared
        input is named) under its consequence set, or with no set on either item. A run
        documenting the value under one class cannot be told from the run both causes change.
        A yielded item's second class is recorded in ``needed``. Where ``self.run`` names a
        run, its lists alone are read (``manifest_run``). Under ``REFERENCE_PRESENT`` a
        ``divergences`` item whose ``WITH_REFERENCE_SETS`` is Ensembl VEP's own set is not
        yielded and not counted in ``applied``: its key goes to ``equal_by_location`` and
        ``equal_by_record``."""

        def over_compared_input(run_input: str) -> bool:
            return input_name is None or run_input == input_name

        runs: list[dict] = []
        if self.run is not None:
            runs.append(manifest_run(manifest, self.run, input_name)[0])
        else:
            flag_sets = manifest.get("flag_sets")
            if not isinstance(flag_sets, list):
                flag_sets = []
            selected = [
                fs for fs in flag_sets if isinstance(fs, dict) and over_compared_input(fs.get("input", CORPUS_INPUT))
            ]
            inherits = not selected or any("divergences" not in fs for fs in selected)
            if over_compared_input(CORPUS_INPUT) or inherits:
                runs.append(manifest)
            runs.extend(selected)
        by_ordinal = input_name is not None
        seen: set[str] = set()
        named: list[tuple[str, dict]] = []
        withheld: dict[tuple, set[str]] = defaultdict(set)
        for run in runs:
            for item in run.get(section) or []:
                if not isinstance(item, dict):
                    continue
                text = json.dumps(item, sort_keys=True)
                if text in seen:
                    continue
                seen.add(text)
                klass, second = _s(item, "expected_divergence"), _s(item, SECOND_EXPECTED_DIVERGENCE)
                self.available[klass] += 1
                if second and not {klass, second} <= masks:
                    for key in self._documented_keys(item, by_ordinal):
                        withheld[key].add(consequence_set(_s(item, "consequence_set")))
                elif klass in masks:
                    named.append((klass, item))
        for klass, item in named:
            terms = consequence_set(_s(item, "consequence_set"))
            held = {t for key in self._documented_keys(item, by_ordinal) for t in withheld.get(key, ())}
            if held and (not terms or "" in held or terms in held):
                continue
            if (
                section == "divergences"
                and self.reference_sequence == REFERENCE_PRESENT
                and with_reference_sets(item) == [consequence_set(_s(item, "vep_consequence_set"))]
            ):
                allele, feature = _s(item, "allele"), _s(item, "feature")
                self.equal_by_location[(_s(item, "location"), allele, feature)] = True
                for r in [] if self.foreign_ordinals else _record_indices(item):
                    self.equal_by_record[(r, allele, feature)] = True
                continue
            self.applied[klass] += 1
            self.read[klass] += 1
            if _s(item, SECOND_EXPECTED_DIVERGENCE):
                self.needed.add(_s(item, SECOND_EXPECTED_DIVERGENCE))
            yield klass, item

    def leave_unmet_items_uncounted(self) -> None:
        """Takes the manifest items counted as they were read out of ``applied``, for a
        comparison none of whose entries names a Location or an input record the lists number
        (the VCF and JSON formats of an input that is not a VCF, or that no run of the
        manifest names): a manifest item is looked up by one or the other, so none meets a
        row there, and a row is counted where a class's own predicate takes it, once."""
        for klass, count in self.read.items():
            self.applied[klass] -= count
            if self.applied[klass] <= 0:
                del self.applied[klass]
        self.read.clear()

    @staticmethod
    def _lookup(
        by_location: dict[tuple[str, str, str], T], by_record: dict[tuple[int, str, str], T], e: Entry
    ) -> T | None:
        """The value documented on the entry's key: by Location when the entry has one, then
        by record ordinal when it has one."""
        allele, feature = dash(e.allele), dash(e.feature)
        if e.location is not None:
            value = by_location.get((e.location, allele, feature))
            if value is not None:
                return value
        if e.record is not None:
            return by_record.get((e.record, allele, feature))
        return None

    def _divergent_fields(self, e: Entry) -> dict[str, str] | None:
        """The field values documented on the entry's key under its consequence set, else
        under no set (an item the manifest documents without a ``consequence_set``)."""
        allele, feature = dash(e.allele), dash(e.feature)
        for terms in (consequence_set(e.consequence), ""):
            if e.location is not None:
                value = self.fields_by_location.get((e.location, allele, feature, terms))
                if value is not None:
                    return value
            if e.record is not None:
                value = self.fields_by_record.get((e.record, allele, feature, terms))
                if value is not None:
                    return value
        return None

    def divergent_field(self, e: Entry, fname: str) -> str | None:
        """The documented candidate value of ``fname`` on the entry's key; the manifest
        spells the default-format name, the JSON format its lower-case form."""
        m = self._divergent_fields(e)
        if m is None:
            return None
        if fname in m:
            return m[fname]
        for k, v in m.items():
            if k.lower() == fname.lower():
                return v
        return None

    def has_divergent_fields(self, e: Entry) -> bool:
        return self._divergent_fields(e) is not None

    def divergent(self, e: Entry) -> list[str] | None:
        """The documented candidate consequence sets of the entry's key, if documented."""
        return self._lookup(self.by_location, self.by_record, e)

    def equal_with_reference(self, e: Entry) -> bool:
        """Whether an item states the entry's key as equal on both engines in a run with a
        reference sequence, the comparison's setting: the entry is compared like any other
        and no class accepts a difference of its Consequence."""
        return self._lookup(self.equal_by_location, self.equal_by_record, e) is not None

    def extra(self, e: Entry) -> bool:
        """Whether the entry's key is documented as one only the candidate emits."""
        allele, feature = dash(e.allele), dash(e.feature)
        if e.location is not None and (e.location, allele, feature) in self.extra_by_location:
            return True
        return e.record is not None and (e.record, allele, feature) in self.extra_by_record

    def record_divergent(self, record: int) -> bool:
        """Whether a documented key belongs to the record, so its summary values may differ."""
        return any(r == record for r, _, _ in self.by_record) or any(
            r == record for r, _, _ in self.extra_by_record
        )

    def reference_skipped(self, record: int) -> bool:
        return record in self.skipped_records


def populated_counts(entries: list[Entry]) -> Counter:
    """Per field, the entries carrying a value other than ``-``."""
    counts: Counter[str] = Counter()
    for e in entries:
        for name, value in e.fields.items():
            if dash(value) != "-":
                counts[name] += 1
    return counts


def extra_pairs(extra: str) -> dict[str, str]:
    """The ``key=value`` pairs of a default-format ``Extra`` value."""
    if extra in ("-", ABSENT):
        return {}
    return dict(kv.partition("=")[::2] for kv in extra.split(";"))


NEAREST_FIELDS = ("NEAREST", "nearest")


def nearest_list(value: str) -> list[str]:
    """The values of a ``NEAREST`` field as every format joins them (``,`` or ``&``)."""
    if value in ("-", ABSENT):
        return []
    return [v for v in re.split(r"[,&]", value) if v]


class NearestSites:
    """The protein-coding transcription start sites of a candidate JSON cache, by transcript
    stable id, gene stable id and gene symbol: a transcript's site is its start on the forward
    strand and its end on the reverse (``AnnotationType/Transcript.pm:749``), a gene's or
    symbol's the sites of its transcripts."""

    def __init__(self, sites: dict[str, set[int]]) -> None:
        self.sites = sites

    @classmethod
    def from_cache(cls, cache: Path) -> NearestSites:
        transcripts = cache / "transcripts"
        if not transcripts.is_dir():
            raise ComparatorError(f"--nearest-sites: {cache} has no transcripts/ directory")
        sites: dict[str, set[int]] = defaultdict(set)
        for shard in sorted(transcripts.glob("*/*.json")) + sorted(transcripts.glob("*/*.json.gz")):
            try:
                parsed = json.loads(read_text(shard))
            except ValueError as exc:
                raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
            for t in cls.shard_transcripts(shard, parsed):
                if t.get("biotype") != "protein_coding":
                    continue
                forward = str(t.get("strand", "1")) == "1"
                site = cls.coordinate(shard, t, "start" if forward else "end")
                for key in ("stable_id", "gene_stable_id", "gene_symbol"):
                    value = t.get(key)
                    if value:
                        sites[value].add(site)
        return cls(dict(sites))

    @staticmethod
    def shard_transcripts(shard: Path, parsed: object) -> list[dict]:
        """The transcript objects of a parsed shard: a bare array, or a converted cache's
        range-keyed map of arrays; a ``null`` in either is a placeholder and is dropped."""
        if isinstance(parsed, dict):
            if not all(isinstance(wrapped, list) for wrapped in parsed.values()):
                raise ComparatorError(f"{shard}: a range-keyed map value is not a transcript array")
            parsed = [t for wrapped in parsed.values() for t in wrapped]
        if not isinstance(parsed, list):
            raise ComparatorError(f"{shard} is neither a transcript array nor a range-keyed map")
        transcripts = [t for t in parsed if t is not None]
        if not all(isinstance(t, dict) for t in transcripts):
            raise ComparatorError(f"{shard}: a transcript is not an object")
        return transcripts

    @staticmethod
    def coordinate(shard: Path, transcript: dict, key: str) -> int:
        """The ``key`` coordinate of a transcript, an integer or a string of digits."""
        value = transcript.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, str)) or not str(value).isdecimal():
            raise ComparatorError(f"{shard}: transcript {transcript.get('stable_id')!r} has no integer {key}")
        return int(value)

    def agree(self, expected: str, actual: str) -> bool:
        """Whether two ``NEAREST`` values name the same winning sites: every expected value is
        among the actual values and every actual value starts at a site of an expected value;
        a value the cache lacks starts nowhere, so it never agrees."""
        expected_values = nearest_list(expected)
        actual_values = nearest_list(actual)
        if not expected_values or not actual_values:
            return expected == actual
        if any(e not in actual_values for e in expected_values):
            return False
        expected_sites = set().union(*(self.sites.get(e, set()) for e in expected_values))
        return all(self.sites.get(a, set()) & expected_sites for a in actual_values)


class OrganelleCodonTables:
    """The transcripts a candidate JSON cache files under a codon table that reads some codon
    unlike the standard one, by stable id: a transcript's
    ``variation_effect_feature_cache.codon_table``, a number in a native cache and the dump's
    string in a converted one, outside ``STANDARD_READING_CODON_TABLES``. A chromosome's shards are read
    on the first lookup of an entry whose Location names its directory, every shard on the
    first lookup of an entry without one (a VCF or JSON entry)."""

    def __init__(self, cache: Path) -> None:
        self.root = cache / "transcripts"
        if not self.root.is_dir():
            raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
        self._tables: dict[str | None, dict[str, int]] = {}

    def _read(self, chromosome: str | None) -> dict[str, int]:
        directories = [self.root / chromosome] if chromosome is not None else sorted(self.root.iterdir())
        tables: dict[str, int] = {}
        for directory in directories:
            for shard in sorted(directory.glob("*.json")) + sorted(directory.glob("*.json.gz")):
                try:
                    parsed = json.loads(read_text(shard))
                except ValueError as exc:
                    raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
                for t in NearestSites.shard_transcripts(shard, parsed):
                    table = str((t.get("variation_effect_feature_cache") or {}).get("codon_table", ""))
                    if t.get("stable_id") and table.isdecimal() and int(table) not in (0, *STANDARD_READING_CODON_TABLES):
                        tables[t["stable_id"]] = int(table)
        return tables

    def table(self, e: Entry) -> int | None:
        """The codon table of the entry's transcript, or ``None`` for a transcript the cache
        files under a table reading as the standard one or does not hold; a Feature carrying a version the
        cache's stable id lacks (``--transcript_version``) is looked up without it."""
        chromosome = e.location.split(":", 1)[0] if e.location else None
        if chromosome is not None and not (self.root / chromosome).is_dir():
            chromosome = None
        if chromosome not in self._tables:
            self._tables[chromosome] = self._read(chromosome)
        tables = self._tables[chromosome]
        return tables.get(e.feature, tables.get(re.sub(r"\.\d+$", "", e.feature)))


def extra_accepted(expected: str, actual: str, accept: Callable[[str, str, str], bool]) -> bool:
    """Whether ``actual`` is ``expected`` with at least one differing ``key=value`` pair that
    ``accept(key, expected value, actual value)`` allows substituted, the keys and their
    order kept; a difference in the order of the pairs is not accepted."""
    ap = extra_pairs(actual)
    rewritten = []
    substituted = False
    for kv in expected.split(";"):
        k, sep, v = kv.partition("=")
        if sep and k in ap and v != ap[k] and accept(k, v, ap[k]):
            kv = f"{k}={ap[k]}"
            substituted = True
        rewritten.append(kv)
    return substituted and ";".join(rewritten) == actual


def extra_nearest_accepted(expected: str, actual: str, sites: NearestSites) -> bool:
    """``extra_accepted`` for an agreeing ``NEAREST`` value."""
    return extra_accepted(expected, actual, lambda k, v, a: k == "NEAREST" and sites.agree(v, a))


def extra_pairs_accepted(expected: str, actual: str, field_pairs: dict[str, set[tuple[str, str]]]) -> bool:
    """``extra_accepted`` for the accepted field pairs."""
    return extra_accepted(expected, actual, lambda k, v, a: (v, a) in field_pairs.get(k, ()))


def extra_order_only_accepted(expected: str, actual: str, order_only: dict[str, str]) -> bool:
    """``extra_accepted`` for the set-valued columns, their items agreeing as a multiset."""
    return extra_accepted(expected, actual, lambda k, v, a: k in order_only and same_items(v, a, order_only[k]))


def extra_class_accepted(expected: str, actual: str, accept: Callable[[str, str, str], bool]) -> bool:
    """Whether the two ``Extra`` texts differ only in pairs ``accept(key, expected value,
    actual value)`` allows: a value that differs, a pair one side alone carries (its other
    value read as ``-``, since Ensembl VEP and the candidate both leave an undefined key out of
    the text), the keys both sides carry in the same order. Unlike ``extra_accepted`` it admits
    a pair the reference lacks, the shape of a field the candidate alone fills."""
    expected_map = extra_pairs(expected)
    actual_map = extra_pairs(actual)
    changed = False
    for key in sorted(set(expected_map) | set(actual_map)):
        ev, av = expected_map.get(key, "-"), actual_map.get(key, "-")
        if ev == av:
            continue
        if not accept(key, ev, av):
            return False
        changed = True
    common_expected = [k for k in expected_map if k in actual_map]
    common_actual = [k for k in actual_map if k in expected_map]
    return changed and common_expected == common_actual


def as_number(value: str) -> float | None:
    """The value as a finite number, or ``None`` when it is not one."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) or math.isinf(number) else number


def field_named(name: str, patterns: Collection[str]) -> bool:
    """Whether a field name is one of ``patterns``, each a full name or a prefix ending in
    ``_``, the case of both set aside (the JSON format lower-cases every key)."""
    lowered = name.lower()
    return any(lowered.startswith(p.lower()) if p.endswith("_") else lowered == p.lower() for p in patterns)


def record_label(e: Entry) -> str:
    """The input record an entry belongs to: the record ordinal of a VCF or JSON entry, and
    in the text formats the Location with the ``Uploaded_variation``, which tells two records
    at one Location apart (two records sharing both are read as one, which only narrows the
    rules that group a record's rows)."""
    label = e.key.split("|", 1)[0]
    if e.location is None:
        return label
    uploaded = dash(e.fields.get("Uploaded_variation", "-"))
    return label if uploaded == "-" else f"{label} {uploaded}"


def entry_kind(e: Entry) -> str:
    """The last segment of an entry's key without its duplicate counter: the Feature_type
    (default, tab, VCF), the JSON list's kind, or ``record``."""
    return e.key.rsplit("|", 1)[-1].partition("#")[0]


class RecordRows:
    """The entries of one output grouped by input record (``record_label``): the row group
    ``compare_entries`` hands a field class whose shape spans a record's rows
    (``RECORD_READING_CLASSES``). ``input_lines`` are the input VCF's data lines, read for the
    default and tab rows that name their record by an identifier. A line carries the first
    ``;``-joined element of its ID column, the name both engines print for its rows, and no
    name where that element is a dot, empty or ``0``."""

    def __init__(self, entries: list[Entry], input_lines: list[str] | None = None) -> None:
        self.by_record: dict[str, list[Entry]] = defaultdict(list)
        for e in entries:
            self.by_record[record_label(e)].append(e)
        # The input lines by identifier, each stored with the identifier as its third column.
        self.lines_by_name: dict[str, list[list[str]]] = defaultdict(list)
        for columns in (line.split("\t") for line in input_lines or ()):
            if len(columns) >= 5:
                columns[2] = columns[2].split(";", 1)[0]
                if columns[2] not in (".", "", "0"):
                    self.lines_by_name[columns[2]].append(columns)
        # Per label, the line of each of its entries that several lines were asked about, kept
        # under the entry itself (``_unasked`` until then), and the label's entries under each
        # line in their order: asking a label's entries for their lines again at every look-up
        # costs the square of its entries times the lines.
        self._line_by_entry: dict[str, dict[int, list[str] | None]] = {}
        self._entries_by_line: dict[str, dict[int, list[Entry]]] = {}
        self._unasked: list[str] = []

    def line_of(self, e: Entry) -> list[str] | None:
        """The input line of the record a default or tab row under an identifier belongs to:
        the one line that carries the identifier, or of several the one that writes the row
        at its Location with its allele (``line_writes_row``). ``None`` for an identifier no
        line carries, and where none or more than one of several lines writes the row. Several
        lines are asked about an entry of this output once, when its own line is wanted."""
        lines = self.lines_by_name.get(e.fields.get("Uploaded_variation", ""), ())
        if len(lines) < 2:
            return lines[0] if lines else None
        label = record_label(e)
        asked = self._line_by_entry.get(label)
        if asked is None:
            asked = self._line_by_entry[label] = dict.fromkeys(map(id, self.by_record.get(label, ())), self._unasked)
        line = asked.get(id(e), self._unasked)
        if line is self._unasked:
            writing = [columns for columns in lines if line_writes_row(columns, e)]
            line = writing[0] if len(writing) == 1 else None
            if id(e) in asked:
                asked[id(e)] = line
        return line

    def _written_by(self, label: str, line: list[str]) -> list[Entry]:
        """The entries of ``label`` whose line (``line_of``) is ``line``, in the label's order."""
        by_line = self._entries_by_line.get(label)
        if by_line is None:
            by_line = {}
            for row in self.by_record.get(label, ()):
                written = self.line_of(row)
                if written is not None:
                    by_line.setdefault(id(written), []).append(row)
            self._entries_by_line[label] = by_line
        return by_line.get(id(line), [])

    def one_record(self, e: Entry) -> bool:
        """Whether ``e``'s label names one input record: a VCF or JSON entry's ordinal, a
        default or tab row whose identifier-less ``Uploaded_variation`` spells the record's
        alleles, or one whose identifier names one line of the input (``line_of``). Under an
        identifier two records at one Location share a label, so without the input their
        rows do not tell a multi-allelic record from two one-allele records."""
        if e.location is None:
            return True
        return uploaded_alleles(e.fields.get("Uploaded_variation", "")) is not None or self.line_of(e) is not None

    def of(self, e: Entry) -> list[Entry] | None:
        """The entries of ``e``'s input record, or ``None`` where the label may name two
        records (``one_record``). The record's own entry is among them in the VCF and JSON
        formats, and for a default or tab row under an identifier one built from its input
        line (``CHROM``, ``POS``, ``REF``, ``ALT``) stands in for it. Where several lines
        carry the identifier, the rows are those of the label that the same line writes."""
        if not self.one_record(e):
            return None
        label = record_label(e)
        rows = self.by_record.get(label, [])
        line = self.line_of(e) if e.location is not None else None
        if line is None:
            return rows
        if len(self.lines_by_name[line[2]]) > 1:
            rows = self._written_by(label, line)
        fields = {"Uploaded_variation": line[2], "CHROM": line[0], "POS": line[1], "REF": line[3], "ALT": line[4]}
        return [*rows, Entry(f"{e.location}|record", None, e.location, "", "", "", fields)]


def line_writes_row(columns: list[str], e: Entry) -> bool:
    """Whether the input line ``columns`` (``CHROM``, ``POS``, ``ID``, ``REF``, ``ALT``)
    writes the default or tab row ``e``: it stands on the row's chromosome, a ``chr`` prefix
    set aside, and one of its alternates is its reference allele with the row's allele in
    place of the Location's bases, or between the Location's two positions
    (``replaced_reference``). Each alternate is asked by itself, so a line two of whose
    alternates give the row in two ways writes it too."""
    if e.location is None:
        return False
    position = vep_rs_pos(columns[1])
    if position is None or columns[0].removeprefix("chr") != e.location.partition(":")[0].removeprefix("chr"):
        return False
    alternates = perl_split(columns[4], ",")
    return any(replaced_reference(position, columns[3], [alternate], e.location, e.allele) is not None for alternate in alternates)


class GwasCatalogPairs:
    """The places where the two engines read a GWAS catalog's rows of two or more SNPs
    differently, each a (chromosome without ``chr``, position, risk allele, study accession):
    ``reference_only`` where Ensembl VEP hands a location another SNP's risk allele
    (GWAS.pm:338-340), ``candidate_only`` where it keeps the space after an allele written
    before `` x `` and so never matches the SNP's own (GWAS.pm:343-348)."""

    def __init__(self) -> None:
        self.reference_only: set[tuple[str, str, str, str]] = set()
        self.candidate_only: set[tuple[str, str, str, str]] = set()

    @classmethod
    def from_catalog(cls, path: Path) -> GwasCatalogPairs:
        """The pairs of the associations file ``path``, read by its header's column names."""
        lines = text_lines(read_text(path))
        header = lines[0].split("\t") if lines else []
        if not all(name in header for name in GWAS_CATALOG_COLUMNS):
            raise ComparatorError(f"{path}: not a GWAS catalog associations file (columns {', '.join(GWAS_CATALOG_COLUMNS)})")
        columns = [header.index(name) for name in GWAS_CATALOG_COLUMNS]
        snps = header.index(GWAS_SNPS_COLUMN) if GWAS_SNPS_COLUMN in header else None
        pairs = cls()
        for line in lines[1:]:
            cells = line.split("\t")
            if len(cells) > max(columns):
                pairs.add_row(*(cells[i] for i in columns), snps=cells[snps] if snps is not None and snps < len(cells) else "")
        return pairs

    def add_row(self, chromosomes: str, positions: str, risk_alleles: str, study: str, snps: str = "") -> None:
        """One catalog row: its lists split as the plugin splits them, the locations
        stripped, each risk allele the text after its token's first ``-`` as written. A row
        of several locations whose SNPs, the ``rs`` ids of ``snps`` or else its risk alleles,
        differ from them in count is one the candidate skips, while the reference writes
        every SNP at every location (GWAS.pm:237 compares the chromosomes with the positions
        alone), so each of its alleles is the reference's alone at each location. A row
        whose chromosomes and positions differ in count, or whose risk alleles differ in
        count from as many SNPs as locations, pairs nothing. Neither does a row whose
        risk-allele cell holds a ``?``: both engines read that cell as empty (GWAS.pm:261)
        and write nothing from the row."""
        chromosome_list = [c.strip().removeprefix("chr") for c in GWAS_LIST_SPLIT_RE.split(chromosomes)]
        position_list = [p.strip() for p in GWAS_LIST_SPLIT_RE.split(positions)]
        written = [token.split("-")[1] if "-" in token else "" for token in GWAS_LIST_SPLIT_RE.split(risk_alleles)]
        if len(position_list) < 2 or len(chromosome_list) != len(position_list) or "?" in risk_alleles:
            return
        places = list(zip(chromosome_list, position_list))
        matched_as_written = {allele for allele in written if re.fullmatch(r"[ACGT]+", allele)}
        if (len(GWAS_RS_ID_RE.findall(snps)) or len(written)) != len(places):
            self.reference_only.update((*place, allele, study) for place in places for allele in matched_as_written)
            return
        if len(written) != len(places):
            return
        for (chromosome, position), allele in zip(places, written):
            own = {allele.strip()} if re.fullmatch(r"[ACGT]+", allele.strip()) else set()
            self.reference_only.update((chromosome, position, other, study) for other in matched_as_written - own)
            self.candidate_only.update((chromosome, position, other, study) for other in own - matched_as_written)


def gwas_site(e: Entry, record_rows: Collection[Entry] | None) -> tuple[str, str] | None:
    """The chromosome, without ``chr``, and position at which the plugin matches a row's
    allele (GWAS.pm:176-195 compare the record's alleles at its start with a risk allele at
    the catalog's position): in the default and tab formats a one-base Location's position, or
    the first position of a Location spanning as many bases as the row's sequence allele, a
    substitution of several bases, which an insertion the parser read from the record's
    alleles is not; else the ``CHROM`` of the VCF format's record entry with the start the
    parser gives the row's allele from the line's ``POS``, ``REF`` and ``ALT``
    (``compare_vep_outputs.py`` ``gvcf_site_rows``; the ``POS`` itself for an allele it does
    not place), or the ``seq_region_name`` and ``start`` of the JSON format's among
    ``record_rows``. A VCF or JSON row is held to the Location's condition where its record
    states the span the allele replaces, the parser's start and end of the row's allele or the
    JSON object's ``start`` and ``end``: a span of another length than one base is read only
    for a sequence allele of as many bases that is no insertion's or deletion's, so an
    insertion's row, whose span is empty, states no site. ``None`` where the row states none."""
    if e.location is not None:
        chromosome, _, span = e.location.rpartition(":")
        start, dash_between, end = span.partition("-")
        if not chromosome or not start.isdigit():
            return None
        if dash_between:
            bases = re.fullmatch(r"[ACGTN]+", e.allele) is not None
            if not (bases and e.indel is not True and end.isdigit() and int(end) - int(start) + 1 == len(e.allele)):
                return None
        return chromosome.removeprefix("chr"), start
    label = record_label(e)
    for r in record_rows or ():
        if entry_kind(r) == "record" and record_label(r) == label:
            chromosome = next((r.fields[f] for f in ("CHROM", "seq_region_name") if f in r.fields), "")
            position = next((r.fields[f] for f in ("POS", "start") if f in r.fields), "")
            span: int | None = None
            if "POS" in r.fields and position.isdigit():
                script_dir = str(Path(__file__).resolve().parent)
                if script_dir not in sys.path:
                    sys.path.insert(0, script_dir)
                import compare_vep_outputs as cvo

                # Parser/VCF.pm:300-336 drop the first base an indel record's alleles share and
                # move its start one base on, so the row's start is not the line's POS there.
                rows = cvo.gvcf_site_rows(int(position), r.fields.get("REF", ""), r.fields.get("ALT", "").split(","))
                starts = {str(start) for allele, start, _ in rows if allele == e.allele.upper()}
                if len(starts) == 1:
                    position = starts.pop()
                    spans = {end - start + 1 for allele, start, end in rows if allele == e.allele.upper()}
                    span = spans.pop() if len(spans) == 1 else None
            elif position.isdigit() and r.fields.get("end", "").isdigit():
                span = int(r.fields["end"]) - int(position) + 1
            if span is not None and span != 1:
                bases = re.fullmatch(r"[ACGTN]+", e.allele) is not None
                if not (bases and e.indel is not True and span == len(e.allele)):
                    return None
            return (chromosome.removeprefix("chr"), position) if chromosome and position.isdigit() else None
    return None


def gwas_pair_named(e: Entry, studies: str, site: tuple[str, str], pairs: set[tuple[str, str, str, str]]) -> bool:
    """Whether ``pairs`` holds the row's site and allele under one of the studies the
    ``GWAS_study`` value ``studies`` names (joined by ``,``; ``&`` in the VCF format)."""
    return not blank(studies) and any((*site, e.allele, study) in pairs for study in re.split(r"[,&]", studies))


def record_alleles(e: Entry, record_rows: Collection[Entry]) -> set[str]:
    """The alleles the rows of ``e``'s record carry, its own among them; an entry of
    ``record_rows`` under another record's label (``record_label``) is not read."""
    label = record_label(e)
    return {e.allele, *(r.allele for r in record_rows if entry_kind(r) != "record" and record_label(r) == label)}


def sibling_rows(e: Entry, record_rows: Collection[Entry]) -> list[Entry]:
    """The rows of ``e``'s record on its feature under another allele; an entry of
    ``record_rows`` under another record's label (``record_label``) is none of them."""
    kind, label = entry_kind(e), record_label(e)
    return [
        r for r in record_rows
        if r.allele != e.allele and r.feature == e.feature and entry_kind(r) == kind and record_label(r) == label
    ]


def record_chromosome(e: Entry, record_rows: Collection[Entry] | None) -> str | None:
    """The chromosome of ``e``'s record as its output spells it: the Location's in the default
    and tab formats, else the ``CHROM`` of the VCF format's record entry or the
    ``seq_region_name`` of the JSON format's among ``record_rows``, under ``e``'s own label
    (``record_label``); ``None`` where none is carried."""
    if e.location is not None:
        return e.location.rpartition(":")[0] or None
    label = record_label(e)
    for r in record_rows or ():
        if entry_kind(r) == "record" and record_label(r) == label:
            return next((r.fields[f] for f in ("CHROM", "seq_region_name") if f in r.fields), None)
    return None


def record_indel_row(e: Entry, record_rows: Collection[Entry]) -> bool:
    """``indel_row`` read from the record's own alleles where the parser did not read them and
    ``record_rows`` carries them: the ``REF`` and ``ALT`` of the record entry built from the
    input line of a default or tab row under an identifier. The Location alone reads an
    insertion written at a one-base Location, the ``CA`` of ``C/-/CA``, as a substitution."""
    if e.allele != "-" and e.indel is None:
        label = record_label(e)
        for r in record_rows:
            if entry_kind(r) == "record" and record_label(r) == label and "REF" in r.fields and "ALT" in r.fields:
                shape = allele_is_indel(r.fields["REF"], r.fields["ALT"].split(","), e.allele)
                if shape is not None:
                    return shape
    return indel_row(e)


def record_allele_lengths(e: Entry, record_rows: Collection[Entry] | None) -> tuple[int, int] | None:
    """The lengths of the reference allele and of ``e``'s allele as Ensembl VEP hands them to
    a plugin, an empty allele 0: read from the record's own alleles (an identifier-less
    ``Uploaded_variation``, else the record entry among ``record_rows``: a VCF line's ``REF``
    and ``ALT``, the JSON ``allele_string``), else from the row's text where it states them.
    A deletion's Location spans the deleted bases and a Location of one base or of three or
    more is the reference allele's; a two-base Location is an insertion's two flanking bases
    or a two-base reference allele, so ``None`` there, as for an allele that is no sequence
    and one that no alternate, or more than one, spells."""
    if re.fullmatch(r"[ACGTN]+|-", e.allele) is None:
        return None
    alleles: tuple[str, list[str]] | None = None
    spelt = uploaded_alleles(e.fields.get("Uploaded_variation", ""))
    if spelt is not None:
        alleles = vcf_alleles_as_spelt(*spelt)
    else:
        label = record_label(e)
        record = next((r for r in record_rows or () if entry_kind(r) == "record" and record_label(r) == label), None)
        if record is not None and "REF" in record.fields and "ALT" in record.fields:
            alleles = vcf_alleles_as_spelt(record.fields["REF"], record.fields["ALT"].split(","))
        elif record is not None and "allele_string" in record.fields:
            reference, *alternates = record.fields["allele_string"].split("/")
            alleles = reference, alternates
    if alleles is not None:
        reference, alternates = alleles
        return (len(reference.strip("-")), len(e.allele.strip("-"))) if alternates.count(e.allele) == 1 else None
    if e.location is None:
        return None
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    if not (start.isdigit() and (end.isdigit() or not end)):
        return None
    span = int(end or start) - int(start) + 1
    if e.allele == "-":
        return span, 0
    return (span, len(e.allele)) if span != 2 else None


def cdna_span(e: Entry) -> tuple[int, int] | None:
    """The row's first and last cDNA positions, the lower first: the ``a-b`` or the single
    position of the text formats' ``cDNA_position``, read up to a ``/total`` suffix, or the
    JSON ``cdna_start`` and ``cdna_end``, which an insertion writes the higher first.
    ``None`` where either is not a number."""
    position = e.fields.get("cDNA_position")
    if position is not None:
        start, _, end = position.partition("/")[0].partition("-")
        end = end or start
    else:
        start = e.fields.get("cdna_start", ABSENT).partition("/")[0]
        end = e.fields.get(CDNA_END_FIELD, ABSENT).partition("/")[0]
    if not (start.isdigit() and end.isdigit()):
        return None
    return min(int(start), int(end)), max(int(start), int(end))


def gained_start_position(e: Entry) -> int | None:
    """The ``CapDistanceToStart`` of the start codon the entry reports as gained: its
    ``5UTR_consequence`` is the gained-start term and its ``5UTR_annotation`` names one
    position. ``None`` otherwise."""
    if field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS) != UTRANNOTATOR_GAINED_START_TERM:
        return None
    positions = set(CAP_DISTANCE_RE.findall(field_of(e, UTRANNOTATOR_ANNOTATION_FIELDS)))
    return int(positions.pop()) if len(positions) == 1 else None


def lost_start_position(e: Entry) -> int | None:
    """The ``CapDistanceToStart`` of the start codon the entry reports as lost: its
    ``5UTR_consequence`` is the lost-start term and its ``5UTR_annotation`` names one
    position. ``None`` otherwise."""
    if field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS) != UTRANNOTATOR_LOST_START_TERM:
        return None
    positions = set(CAP_DISTANCE_RE.findall(field_of(e, UTRANNOTATOR_ANNOTATION_FIELDS)))
    return int(positions.pop()) if len(positions) == 1 else None


def allele_on_transcript_strand(e: Entry) -> str | None:
    """The entry's sequence allele as its transcript reads it: as written where the row's
    ``STRAND`` is ``1``, reverse-complemented where it is ``-1``. ``None`` for an allele that
    is no sequence and for a row that states no strand."""
    strand = field_of(e, STRAND_FIELDS)
    if re.fullmatch(r"[ACGTN]+", e.allele) is None or strand not in ("1", "-1"):
        return None
    return reverse_complement(e.allele) if strand == "-1" else e.allele


def gained_start_window_row(e: Entry, a: Entry, record_rows: Collection[Entry] | None) -> bool:
    """Whether the two entries of a row differ as the two gained-start windows make them
    (``utrannotator_gained_start_window``). Both windows open two bases before the allele's
    first base in the mutated UTR and close two bases after its last: the reference's
    (UTRAnnotator.pm:345) counts the reference allele's bases, the candidate's the alternate
    allele's, and each engine reports the first ``ATG`` of its window by the 0-based index of
    its ``A``. So the row is the class's when one engine alone reports the gained start, its
    codon begins inside that engine's allele length and at or past the other's, and the
    reporting engine is the one whose allele is the longer; or when the reference reports it
    at a negative index, which only a window opened before the UTR's first base yields. An
    insertion's row is refused where the reference reports a start codon lost and the
    candidate's gained start stands the allele's length past that codon's index: that codon
    is the reference's own, moved on by the inserted bases."""
    reference_gains = field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS) == UTRANNOTATOR_GAINED_START_TERM
    candidate_gains = field_of(a, UTRANNOTATOR_CONSEQUENCE_FIELDS) == UTRANNOTATOR_GAINED_START_TERM
    at = gained_start_position(e if reference_gains else a)
    if at is None:
        return False
    if reference_gains and at < 0:
        return True
    lengths, span = record_allele_lengths(e, record_rows), cdna_span(e)
    if reference_gains == candidate_gains or lengths is None or span is None:
        return False
    reference_length, allele_length = lengths
    first, last = span
    # The cDNA positions are the reference allele's own, an insertion's its two flanking bases.
    if last - first + 1 != (reference_length or 2):
        return False
    # A codon the reference reports lost that stands the allele's length further on was moved,
    # not made, by the inserted bases.
    lost = lost_start_position(e)
    if candidate_gains and not reference_length and lost is not None and at == lost + allele_length:
        return False
    start = first - 1 if reference_length else first
    reporting, other = (reference_length, allele_length) if reference_gains else (allele_length, reference_length)
    return start + other <= at < start + reporting


def gwas_annotated(e: Entry) -> bool:
    """Whether the entry carries a ``GWAS_`` field with a value: the plugin annotated its row."""
    return any(not blank(value) for name, value in e.fields.items() if field_named(name, (GWAS_PREFIX,)))


def class_field_values(e: Entry, patterns: Collection[str]) -> dict[str, str]:
    """The entry's filled fields among those ``patterns`` name, by name."""
    return {k: v for k, v in e.fields.items() if field_named(k, patterns) and not blank(v)}


def sibling_allele_value(
    e: Entry, patterns: Collection[str], record_rows: Collection[Entry], same_amino_acids: bool = True,
    candidate_rows: Collection[Entry] | None = None,
) -> bool:
    """Whether the reference's row ``e`` carries, in every field ``patterns`` name, the values
    its row of another allele of the record carries on the same feature and, under
    ``same_amino_acids``, beside the same ``Amino_acids``: what a plugin writes that matches
    its table against every alternate allele of the record and keeps the first row that
    matches. The row's fields are read together: a row that shares one field with its sibling
    and carries another of its own is no copy of the sibling's table row. The sibling is the
    allele the table names, whose row carries the values on both engines: it counts only where
    the candidate's row of that allele (``candidate_rows``, the candidate's entries of the
    same record) carries them too. So the named allele's own row, where the candidate lacks
    its value, has no such sibling and stays a difference, and without the candidate's rows
    no row is the shape."""
    values = class_field_values(e, patterns)
    amino_acids = field_of(e, AMINO_ACIDS_FIELDS)
    if not values or (same_amino_acids and blank(amino_acids)) or candidate_rows is None:
        return False
    named = {c.allele for c in sibling_rows(e, candidate_rows) if class_field_values(c, patterns) == values}
    return any(
        class_field_values(s, patterns) == values
        and (not same_amino_acids or field_of(s, AMINO_ACIDS_FIELDS) == amino_acids)
        and s.allele in named
        for s in sibling_rows(e, record_rows)
    )


class RevelRow(NamedTuple):
    """One row of a REVEL table as REVEL.pm:207-241 read it: the cell of the run's position
    column as written, the two allele cells (``None`` for a cell at or past the row's last
    filled one, the undefined value ``split /\\t/`` leaves), the two amino-acid cells, the
    score cell and the transcript ids cell (``None`` likewise, and for the seven-column layout,
    which has no ids)."""

    position: str
    reference: str | None
    alternate: str | None
    reference_amino_acid: str
    alternate_amino_acid: str
    score: str | None
    transcript_ids: str | None


class RevelTable:
    """The REVEL table a scored run's plugin read (``--revel-table``), read through its tabix
    index as the plugin reads it: the rows at a span of positions, in file order. The header,
    the last of the file's leading ``#`` lines, fixes the layout by its column count
    (REVEL.pm:127-130, :209-241): seven columns are ``chr, hg19_pos, ref, alt, aaref, aaalt,
    REVEL``; eight and nine are ``chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL,
    Ensembl_transcriptid``. The position column is the one the index was built on, which is
    the run's assembly's. ``no_match`` says the run set the plugin's ``no_match``
    (REVEL.pm:107-110), under which no row's transcript ids are read."""

    def __init__(self, path: Path, no_match: bool = False) -> None:
        self.path, self.no_match = path, no_match
        index = Path(f"{path}.tbi")
        try:
            raw = gzip.decompress(index.read_bytes())
            if raw[:4] != b"TBI\x01":
                raise ValueError("not a tabix index")
            references, _, sequence_column, start_column, _, meta, _, names_length = struct.unpack_from("<8i", raw, 4)
            names = [name.decode() for name in raw[36 : 36 + names_length].split(b"\0")[:references]]
            offset = 36 + names_length
            self._chunks: dict[str, dict[int, list[tuple[int, int]]]] = {}
            self._linear: dict[str, list[int]] = {}
            for name in names:
                (bins,) = struct.unpack_from("<i", raw, offset)
                offset += 4
                chunks: dict[int, list[tuple[int, int]]] = {}
                for _ in range(bins):
                    number, count = struct.unpack_from("<Ii", raw, offset)
                    bounds = struct.unpack_from(f"<{2 * count}Q", raw, offset + 8)
                    chunks[number] = list(zip(bounds[::2], bounds[1::2]))
                    offset += 8 + 16 * count
                (intervals,) = struct.unpack_from("<i", raw, offset)
                self._linear[name] = list(struct.unpack_from(f"<{intervals}Q", raw, offset + 4))
                offset += 4 + 8 * intervals
                self._chunks[name] = chunks
        except (OSError, EOFError, ValueError, struct.error, zlib.error) as exc:
            raise ComparatorError(f"--revel-table: cannot read the tabix index {index}: {getattr(exc, 'strerror', None) or exc}") from exc
        self._sequence_column, self._position_column, self._meta = sequence_column - 1, start_column - 1, chr(meta)
        header: list[str] = []
        with open(path, "rb") as handle:
            for line in self._lines(handle, 0, 1 << 62):
                if not line.startswith(self._meta):
                    break
                header = line[1:].split()
        if len(header) == 7:
            self._columns = (2, 3, 4, 5, 6, None)
        elif len(header) in (8, 9):
            self._columns = (3, 4, 5, 6, 7, 8)
        else:
            raise ComparatorError(f"--revel-table: {path} is not a REVEL table (a header of 7, 8 or 9 columns)")
        self._read: dict[tuple[str, int, int], list[RevelRow]] = {}

    def _block(self, handle, offset: int) -> tuple[bytes, int]:
        """The inflated bytes of the BGZF block at the file offset ``offset`` and the offset
        of the block after it; no bytes at the end of the file."""
        handle.seek(offset)
        head = handle.read(12)
        if len(head) < 12:
            return b"", offset
        extra = handle.read(struct.unpack_from("<H", head, 10)[0])
        size, at = None, 0
        while at + 4 <= len(extra):
            length = struct.unpack_from("<H", extra, at + 2)[0]
            if extra[at : at + 2] == b"BC" and length == 2:
                size = struct.unpack_from("<H", extra, at + 4)[0] + 1
            at += 4 + length
        if head[:4] != b"\x1f\x8b\x08\x04" or size is None:
            raise ComparatorError(f"--revel-table: {self.path} is not a BGZF file")
        try:
            return zlib.decompress(handle.read(size - len(extra) - 20), -15), offset + size
        except zlib.error as exc:
            raise ComparatorError(f"--revel-table: cannot read {self.path}: {exc}") from exc

    def _lines(self, handle, start: int, end: int) -> Iterator[str]:
        """The lines that begin at the virtual offset ``start`` or after it and before
        ``end``, without their terminators."""
        block, within = start >> 16, start & 0xFFFF
        data, following = self._block(handle, block)
        carried, begins = b"", start
        while data and begins < end:
            while True:
                newline = data.find(b"\n", within)
                if newline < 0:
                    carried += data[within:]
                    break
                yield (carried + data[within:newline]).decode("utf-8", "replace")
                carried, within = b"", newline + 1
                begins = (block << 16) | within
                if begins >= end:
                    return
            block, within = following, 0
            if not carried:
                begins = block << 16
            data, following = self._block(handle, block)
        if carried and begins < end:
            yield carried.decode("utf-8", "replace")

    def rows(self, chromosome: str, start: int, end: int) -> list[RevelRow]:
        """The table's rows at the positions ``start`` to ``end`` of ``chromosome``, in file
        order: those REVEL.pm:168 fetches for a variant of that span. The chromosome is
        looked up as given, then without and with a ``chr`` prefix."""
        key = (chromosome, start, end)
        if key in self._read:
            return self._read[key]
        bare = chromosome.removeprefix("chr")
        name = next((n for n in (chromosome, bare, f"chr{bare}") if n in self._chunks), None)
        found: list[RevelRow] = []
        if name is not None and start <= end:
            low, high = start - 1, end - 1
            bins = [0]
            for shift, first in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
                bins.extend(range(first + (low >> shift), first + (high >> shift) + 1))
            linear = self._linear[name]
            floor = linear[min(low >> 14, len(linear) - 1)] if linear else 0
            chunks = sorted((max(a, floor), b) for number in bins for a, b in self._chunks[name].get(number, ()) if b > floor)
            merged: list[list[int]] = []
            for a, b in chunks:
                if merged and a <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], b)
                else:
                    merged.append([a, b])
            with open(self.path, "rb") as handle:
                for a, b in merged:
                    for line in self._lines(handle, a, b):
                        cells = line.split("\t")
                        if line.startswith(self._meta) or len(cells) <= max(self._sequence_column, self._position_column):
                            continue
                        leading = re.match(r"\d+", cells[self._position_column])
                        if cells[self._sequence_column] != name or leading is None:
                            continue
                        if int(leading.group()) > end:
                            break
                        if int(leading.group()) >= start:
                            found.append(self._row(cells))
        self._read[key] = found
        return found

    def _row(self, cells: list[str]) -> RevelRow:
        """The cells of one line as the plugin's ``split /\\t/`` leaves them: the trailing
        empty cells dropped, a cell past the last filled one undefined."""
        while cells and cells[-1] == "":
            cells.pop()
        reference, alternate, aaref, aaalt, score, ids = (
            cells[column] if column is not None and column < len(cells) else None for column in self._columns
        )
        return RevelRow(cells[self._position_column], reference, alternate, aaref or "", aaalt or "", score, ids)


def trimmed_allele_keys(reference: str, alternate: str, position: int) -> set[tuple[str, str, int]]:
    """The (reference, alternate, start) of an allele pair after each trimming order
    ``get_matched_variant_alleles`` tries (ensembl-variation Utils/Sequence.pm:1190-1209 and
    :1250-1258, ``trim_sequences`` :965-1038 on the forward strand): the shared bases cut
    from the left and then from the right, and from the right first as well where either
    allele is longer than one base; an allele left empty is written ``-``."""
    keys: set[tuple[str, str, int]] = set()
    for end_first in (False, True) if len(reference) > 1 or len(alternate) > 1 else (False,):
        ref, alt, start = reference, alternate, position
        for left in (False, True) if end_first else (True, False):
            while ref and alt and (ref[0] == alt[0] if left else ref[-1] == alt[-1]):
                ref, alt, start = (ref[1:], alt[1:], start + 1) if left else (ref[:-1], alt[:-1], start)
        keys.add((ref or "-", alt or "-", start))
    return keys


def revel_record(
    e: Entry, record_rows: Collection[Entry] | None
) -> tuple[str, int, int, list[set[tuple[str, str, int]]], int] | None:
    """The record of the row ``e`` as the plugin's match reads it (REVEL.pm:168-188): its
    chromosome, the start and end of its span, the allele keys of each of its alternates
    (``trimmed_allele_keys``) and the place of the row's own alternate among them. The
    alleles are read as Ensembl VEP spells them, an indel record's shared first base dropped
    and its start one base on (``vcf_alleles_as_spelt``): from an identifier-less
    ``Uploaded_variation`` in the default and tab formats, else from the record entry among
    ``record_rows`` under ``e``'s own label (``record_label``), a VCF line's ``CHROM``,
    ``POS``, ``REF`` and ``ALT`` or the JSON ``seq_region_name``, ``start`` and
    ``allele_string``. ``None`` where neither states the record, and where no alternate, or
    more than one, spells the row's allele."""
    name = e.fields.get("Uploaded_variation", "")
    spelt = uploaded_alleles(name)
    stated: tuple[str, int | None, str, list[str]] | None = None
    if spelt is not None and e.location is not None:
        stated = (e.location.rpartition(":")[0], uploaded_position(name), *spelt)
    else:
        label = record_label(e)
        record = next((r for r in record_rows or () if entry_kind(r) == "record" and record_label(r) == label), None)
        fields = record.fields if record is not None else {}
        if all(key in fields for key in ("CHROM", "POS", "REF", "ALT")) and fields["POS"].isdigit():
            reference, alternates = vcf_alleles_as_spelt(fields["REF"], fields["ALT"].split(","))
            stated = (fields["CHROM"], int(fields["POS"]) + (reference != fields["REF"]), reference, alternates)
        elif all(key in fields for key in ("seq_region_name", "start", "allele_string")) and fields["start"].isdigit():
            reference, *alternates = fields["allele_string"].split("/")
            stated = (fields["seq_region_name"], int(fields["start"]), reference, alternates)
    if stated is None:
        return None
    chromosome, start, reference, alternates = stated
    if not chromosome or start is None or alternates.count(e.allele) != 1:
        return None
    keys = [trimmed_allele_keys(reference, alternate, start) for alternate in alternates]
    return chromosome, start, start + len(reference.strip("-")) - 1, keys, alternates.index(e.allele)


def revel_row_keys(row: RevelRow) -> set[tuple[str, str, str]] | None:
    """The keys ``get_matched_variant_alleles`` gives a table row (Utils/Sequence.pm:1219-1232):
    its alleles trimmed as ``trimmed_allele_keys`` trims them, with its position cell as text
    where no base is cut from the left. ``None`` for a row with no reference or alternate
    allele cell, which the match dies on (Utils/Sequence.pm:968); no key for a row whose
    position cell is false to Perl."""
    if row.reference is None or row.alternate is None:
        return None
    if row.position in ("", "0"):
        return set()
    if not row.position.isdigit():
        return {(row.reference or "-", row.alternate or "-", row.position)}
    start = int(row.position)
    return {
        (ref, alt, row.position if at == start else str(at))
        for ref, alt, at in trimmed_allele_keys(row.reference, row.alternate, start)
    }


def revel_row_passes(table: RevelTable, row: RevelRow, feature: str, alternate_residue: str) -> bool:
    """The tests REVEL.pm:191-197 make of a row whose alleles matched: where its transcript
    ids cell is true to Perl and the run did not set ``no_match``, one of its ``;``-separated
    ids is the row's ``Feature``; and its alternate amino acid is the row's."""
    ids = row.transcript_ids
    if ids not in (None, "", "0") and not table.no_match and feature not in ids.split(";"):
        return False
    return row.alternate_amino_acid == alternate_residue


def revel_walk(
    table: RevelTable, chromosome: str, start: int, end: int, keys: set[tuple[str, str, int]], feature: str,
    alternate_residue: str, reference_residue: str | None = None,
) -> RevelRow | None:
    """The first table row, in file order over the record's span ``start`` to ``end``, that
    carries one of the allele keys ``keys`` and passes ``revel_row_passes``: the row whose
    score REVEL.pm:170-199 return for an output row. With ``reference_residue`` the row's
    reference amino acid is tested as well, the candidate's walk. ``None`` where the walk
    ends without such a row, and where it meets a row with no allele cell first, on which
    the plugin's match dies."""
    wanted = {(ref, alt, str(at)) for ref, alt, at in keys}
    for row in table.rows(chromosome, start, end):
        row_keys = revel_row_keys(row)
        if row_keys is None:
            return None
        if not (row_keys & wanted and revel_row_passes(table, row, feature, alternate_residue)):
            continue
        if reference_residue is None or row.reference_amino_acid == reference_residue:
            return row
    return None


def revel_one_base_site(e: Entry) -> tuple[str, int] | None:
    """The chromosome and position of a row whose Location is one position and whose allele
    is one base, read from a default or tab row's Location or from the Location a VCF or JSON
    entry's record gives its allele (``Entry.allele_location``); ``None`` for any other row."""
    location = e.location if e.location is not None else e.allele_location
    if location is None or re.fullmatch(r"[ACGT]", e.allele) is None:
        return None
    chromosome, _, position = location.rpartition(":")
    return (chromosome, int(position)) if chromosome and position.isdigit() else None


def revel_other_reference_amino_acid_row(
    e: Entry, a: Entry | None, expected_value: str, actual_value: str,
    record_rows: Collection[Entry] | None, table: RevelTable | None,
) -> bool:
    """Whether the reference's ``REVEL`` on the row ``e`` is the score of a table row whose
    reference amino acid is not the row's own, and the candidate's is what the corrected walk
    gives (``revel_other_reference_amino_acid_rows``). Both outputs spell the row's
    ``Amino_acids`` alike, as one pair of residues.

    Where the record is stated (``revel_record``), the reference's value must be the score of
    the row its walk reaches over every alternate of the record (``revel_walk``), and that
    row's reference amino acid must be another than the first residue of ``Amino_acids``.
    Where it is not, a one-base row at one position (``revel_one_base_site``) is read
    without it: the rows the reference's walk can reach, whatever other alternates the
    record holds, are the first passing row of the row's own allele and, for each other
    allele, its first passing row before that one; every one of them that carries the
    reference's value must name another reference amino acid, and one must carry it.

    The candidate's value must be the score of the first row of the row's own allele that
    passes and names both residues of ``Amino_acids``, or be empty where there is none.
    Nothing is accepted without the table, without the candidate's row, on a row whose
    walk meets a table row with no allele cell, or on a value no such table row explains."""
    if table is None or a is None or blank(expected_value):
        return False
    amino_acids = field_of(e, AMINO_ACIDS_FIELDS)
    reference_residue, slash, alternate_residue = amino_acids.partition("/")
    if not (slash and reference_residue and alternate_residue) or field_of(a, AMINO_ACIDS_FIELDS) != amino_acids:
        return False
    record = revel_record(e, record_rows)
    reached: list[RevelRow | None]
    own: RevelRow | None
    if record is not None:
        chromosome, start, end, keys, place = record
        reached = [revel_walk(table, chromosome, start, end, set().union(*keys), e.feature, alternate_residue)]
        own = revel_walk(table, chromosome, start, end, keys[place], e.feature, alternate_residue, reference_residue)
    else:
        site = revel_one_base_site(e)
        if site is None:
            return False
        rows = table.rows(*site, site[1])
        if any(revel_row_keys(row) is None for row in rows):
            return False
        passing = [
            row for row in rows
            if row.position == str(site[1]) and len(row.alternate or "") == 1
            and revel_row_passes(table, row, e.feature, alternate_residue)
        ]
        first_own = next((n for n, row in enumerate(passing) if row.alternate == e.allele), len(passing))
        reached = []
        for row in passing[: first_own + 1]:
            if all(row.alternate != earlier.alternate for earlier in reached):
                reached.append(row)
        own = next(
            (row for row in passing if row.alternate == e.allele and row.reference_amino_acid == reference_residue), None
        )
    explaining = [row for row in reached if row is not None and row.score == expected_value]
    if not explaining or any(row.reference_amino_acid == reference_residue for row in explaining):
        return False
    wanted = own.score if own is not None else None
    return blank(actual_value) if blank(wanted or "") else actual_value == wanted


def vcf_alleles_as_spelt(reference: str, alternates: list[str]) -> tuple[str, list[str]]:
    """A VCF line's REF and ALTs as its entries spell them: Parser/VCF.pm:300-323 drop the
    first base of a record that holds an ALT of another length than its REF, a ``*`` among
    them, and whose alleles not spelt with ``*`` all share it, an empty remainder written ``-``."""
    plain = [a for a in alternates if "*" not in a]
    holds_indel = any(len(a) != len(reference) for a in alternates)
    if not (holds_indel and reference and all(a[:1] == reference[:1] for a in plain)):
        return reference, alternates
    return dash(reference[1:]), [a if "*" in a else dash(a[1:]) for a in alternates]


def vcf_first_alternate(reference: str, alternates: list[str]) -> str:
    """The first ALT of a VCF line as its entries spell the allele (``vcf_alleles_as_spelt``)."""
    return vcf_alleles_as_spelt(reference, alternates)[1][0]


def first_alternate(e: Entry, record_rows: Collection[Entry]) -> str | None:
    """The first alternate allele of ``e``'s record where an entry states the record's alleles
    in order, spelt as its rows spell the allele: a text row's identifier-less
    ``Uploaded_variation``, which names a VCF record by the line's own REF and ALT, else the
    ``ALT`` of the VCF format's record entry or the ``allele_string`` of the JSON format's
    among ``record_rows``; ``None`` where none does."""
    spelt = uploaded_alleles(e.fields.get("Uploaded_variation", ""))
    if spelt is not None:
        return vcf_first_alternate(*spelt)
    label = record_label(e)
    for r in record_rows:
        if entry_kind(r) == "record" and record_label(r) == label:
            if "ALT" in r.fields:
                return vcf_first_alternate(r.fields.get("REF", ""), r.fields["ALT"].split(","))
            alleles = r.fields.get("allele_string", "").split("/")
            return alleles[1] if len(alleles) > 1 else None
    return None


def alternates_in_order(e: Entry, record_rows: Collection[Entry] | None) -> list[str] | None:
    """The alternate alleles of ``e``'s record in the record's own order, spelt as its rows
    spell the allele, where an entry states them (``first_alternate`` names the entries);
    ``None`` where none does."""
    spelt = uploaded_alleles(e.fields.get("Uploaded_variation", ""))
    if spelt is not None:
        return vcf_alleles_as_spelt(*spelt)[1]
    label = record_label(e)
    for r in record_rows or ():
        if entry_kind(r) == "record" and record_label(r) == label:
            if "ALT" in r.fields:
                return vcf_alleles_as_spelt(r.fields.get("REF", ""), r.fields["ALT"].split(","))[1]
            return r.fields.get("allele_string", "").split("/")[1:] or None
    return None


def later_alternate(e: Entry, record_rows: Collection[Entry] | None) -> bool:
    """Whether ``e`` is the row of a second or later alternate allele of its record: an entry
    states the record's alternates in order (``alternates_in_order``), one of them alone
    spells the row's allele, and it is not the first."""
    alternates = alternates_in_order(e, record_rows)
    return alternates is not None and alternates.count(e.allele) == 1 and alternates[0] != e.allele


def first_allele_shift(e: Entry, record_rows: Collection[Entry] | None) -> int | None:
    """The HGVS_OFFSET of the reference's row of the first alternate allele of ``e``'s record
    on ``e``'s feature, for the row ``e`` of a later alternate (``later_alternate``); ``None``
    for any other row, without that one row among ``record_rows`` and where it carries no
    non-zero offset."""
    if not later_alternate(e, record_rows):
        return None
    first = alternates_in_order(e, record_rows)[0]
    rows = [r for r in sibling_rows(e, record_rows or ()) if r.allele == first]
    offset = field_of(rows[0], HGVS_OFFSET_FIELDS) if len(rows) == 1 else ABSENT
    return int(offset) if NON_ZERO_INTEGER_RE.fullmatch(offset) else None


def earlier_allele_shift(e: Entry, record_rows: Collection[Entry] | None) -> int | None:
    """The non-zero HGVS_OFFSET of the reference's row, on ``e``'s feature, of the last
    alternate allele before ``e``'s own in its record's order whose row carries one, for the
    row ``e`` of a later alternate (``later_alternate``): the shift whose CDS pair that
    allele's protein notation left on the TranscriptVariation the record's alleles share,
    where each allele's own notation moves it afresh (TranscriptVariationAllele.pm:844-858).
    ``None`` for any other row, where no earlier allele's row carries an offset, and where an
    earlier allele has two rows on the feature among ``record_rows``."""
    if not later_alternate(e, record_rows):
        return None
    alternates = alternates_in_order(e, record_rows)
    siblings = sibling_rows(e, record_rows or ())
    for allele in reversed(alternates[: alternates.index(e.allele)]):
        rows = [r for r in siblings if r.allele == allele]
        if len(rows) > 1:
            return None
        offset = field_of(rows[0], HGVS_OFFSET_FIELDS) if rows else ABSENT
        if NON_ZERO_INTEGER_RE.fullmatch(offset):
            return int(offset)
    return None


def nearest_exon_items(value: str) -> list[tuple[str, str, str, str]] | None:
    """The (exon id, distance, side, length) items of a ``NearestExonJB`` value in any
    format's spelling; ``None`` for a value that is not a list of such items."""
    if NEAREST_EXON_LIST_RE.fullmatch(value) is None:
        return None
    return NEAREST_EXON_ITEM_RE.findall(value)


def cdna_span_closed(e: Entry) -> bool:
    """Whether the row's cDNA position names two positions, both numbers: the ``a-b`` of the
    text formats' ``cDNA_position``, read up to a ``/total`` suffix, or the JSON ``cdna_start``
    and ``cdna_end``. A row that begins or ends outside an exon prints ``?`` or nothing there."""
    position = e.fields.get("cDNA_position")
    if position is not None:
        start, _, end = position.partition("/")[0].partition("-")
    else:
        start = e.fields.get("cdna_start", ABSENT).partition("/")[0]
        end = e.fields.get(CDNA_END_FIELD, ABSENT).partition("/")[0]
    return start.isdigit() and end.isdigit() and start != end


def notation_indel_length(hgvsc: str) -> int | None:
    """The number of bases an HGVSc inserts, duplicates or deletes: the ``ins`` bases, one for
    a single position, else the span between two positions without an intron offset or a
    ``*`` (the numbering has no position 0); ``None`` for any other notation."""
    m = INDEL_NOTATION_RE.fullmatch(hgvsc.partition(":")[2])
    if m is None:
        return None
    start, end, kind = m.groups()
    if kind.startswith("ins"):
        return len(kind) - len("ins")
    if end is None:
        return 1
    if INTEGER_RE.fullmatch(start) is None or INTEGER_RE.fullmatch(end) is None:
        return None
    first, last = int(start), int(end)
    length = last - first + 1
    if first < 0 < last:
        length -= 1
    return length


def lof_info_items(value: str) -> list[tuple[str, str]] | None:
    """The ``KEY:VALUE`` items of a ``LoF_info`` value in order, joined by ``,`` in the default,
    tab and JSON formats and by ``&`` in the VCF format; ``None`` for an empty value."""
    if blank(value):
        return None
    return [(key, item) for key, _, item in (part.partition(":") for part in re.split(r"[,&]", value))]


def gerp_totals(e: Entry, a: Entry | None) -> tuple[float, float, dict[str, str]] | None:
    """The reference's and the candidate's ``GERP_DIST`` on a row whose two ``LoF_info`` values
    carry the same items in the same order and differ in that item alone, two different
    numbers, with the row's other items by key; ``None`` for any other row."""
    if a is None:
        return None
    expected, actual = lof_info_items(field_of(e, LOF_INFO_FIELDS)), lof_info_items(field_of(a, LOF_INFO_FIELDS))
    if expected is None or actual is None or [key for key, _ in expected] != [key for key, _ in actual]:
        return None
    if [key for (key, value), (_, other) in zip(expected, actual) if value != other] != ["GERP_DIST"]:
        return None
    items = dict(expected)
    reference_total, candidate_total = as_number(items["GERP_DIST"]), as_number(dict(actual)["GERP_DIST"])
    if reference_total is None or candidate_total is None or reference_total == candidate_total:
        return None
    return reference_total, candidate_total, items


def lof_filter_items(value: str) -> list[str]:
    """The items of a ``LoF_filter`` value, joined by ``,`` in the default, tab and JSON
    formats and by ``&`` in the VCF format; none for an empty value."""
    return [] if blank(value) else re.split(r"[,&]", value)


def lof_filter_item_side(e: Entry, a: Entry | None, item: str) -> str | None:
    """The side, ``reference`` or ``candidate``, whose ``LoF_filter`` carries ``item`` on a row
    whose two ``LoF_filter`` values differ by that item alone and whose ``LoF`` reads ``LC``
    beside the item and, on the other side, ``LC`` beside another item and ``HC`` beside
    none; ``None`` for any other row."""
    if a is None:
        return None
    expected_items = lof_filter_items(field_of(e, LOF_FILTER_FIELDS))
    actual_items = lof_filter_items(field_of(a, LOF_FILTER_FIELDS))
    if (item in expected_items) == (item in actual_items):
        return None
    expected_call, actual_call = field_of(e, LOF_FIELDS), field_of(a, LOF_FIELDS)
    if item in actual_items:
        side = "candidate"
        carrying_items, carrying_call = actual_items, actual_call
        other_items, other_call = expected_items, expected_call
    else:
        side = "reference"
        carrying_items, carrying_call = expected_items, expected_call
        other_items, other_call = actual_items, actual_call
    if carrying_items.count(item) != 1 or [i for i in carrying_items if i != item] != other_items:
        return None
    if carrying_call != "LC" or other_call != ("LC" if other_items else "HC"):
        return None
    return side


def annotation_groups(value: str) -> tuple[tuple[str, ...], ...]:
    """A text-format ``5UTR_annotation`` as its groups of ``key=value`` pairs, the order of
    the groups and of each group's pairs set aside (UTRAnnotator.pm:282-283 walk both as
    hashes)."""
    return tuple(sorted(tuple(sorted(group.split(":"))) for group in value.split(",")))


class UtrannotatorSeedSweep:
    """What Ensembl VEP's UTRAnnotator plugin writes on each row of an input across runs
    under different Perl hash seeds (``PERL_HASH_SEED`` with ``PERL_PERTURB_KEYS=0``), each
    of the five flags of UTRAnnotator.pm:262-267 walked last under one of them: per row of
    the default-format outputs (the Location with the ``Uploaded_variation``, the Allele and
    the Feature), every ``5UTR_consequence`` term with the ``5UTR_annotation`` values written
    beside it. The plugin keeps the last flag that fired, so the terms of a row are the
    analyses that fire on it."""

    def __init__(self) -> None:
        self.written: dict[tuple[str, str, str], dict[str, set[tuple[tuple[str, ...], ...]]]] = defaultdict(dict)

    @classmethod
    def from_outputs(cls, paths: Collection[Path]) -> UtrannotatorSeedSweep:
        """The rows of the default-format outputs ``paths``, one run each. A file whose header
        line does not name the columns a row is keyed by is refused: read as an output, it
        would write no term on any row."""
        sweep = cls()
        for path in paths:
            text = read_text(path)
            headers = (line for line in text_lines(text) if line.startswith("#") and not line.startswith("##"))
            if not {"Location", "Allele", "Feature"} <= set(next(headers, "")[1:].split("\t")):
                raise ComparatorError(f"{path}: not a default-format output (columns Location, Allele, Feature)")
            sweep.add_output(parse_default(text))
        return sweep

    @staticmethod
    def row(e: Entry) -> tuple[str, str, str]:
        return record_label(e), e.allele, e.feature

    def add_output(self, entries: list[Entry]) -> None:
        """One run's entries: each row's term, where it is one of the plugin's five."""
        for e in entries:
            term = field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS)
            if term in UTRANNOTATOR_TERMS:
                annotation = annotation_groups(field_of(e, UTRANNOTATOR_ANNOTATION_FIELDS))
                self.written[self.row(e)].setdefault(term, set()).add(annotation)

    def gained_start(self, e: Entry) -> int | None:
        """The one ``CapDistanceToStart`` the runs write beside the gained-start term on
        ``e``'s row; ``None`` where no run writes that term there or they write more than one
        position. A VCF or JSON entry names no Location and matches no row."""
        written = self.written.get(self.row(e)) if e.location is not None else None
        annotations = written.get(UTRANNOTATOR_GAINED_START_TERM) if written is not None else None
        if not annotations:
            return None
        positions = {pair for groups in annotations for group in groups for pair in group if CAP_DISTANCE_RE.fullmatch(pair)}
        return int(CAP_DISTANCE_RE.fullmatch(positions.pop()).group(1)) if len(positions) == 1 else None

    def writes(self, e: Entry, a: Entry, with_annotation: bool) -> bool:
        """Whether the runs write, on ``e``'s row, the term of the reference's entry ``e``
        and the term of the candidate's entry ``a``, and under ``with_annotation`` beside
        each the entry's own ``5UTR_annotation``. A VCF or JSON entry names no Location and
        matches no row."""
        written = self.written.get(self.row(e)) if e.location is not None else None
        if written is None:
            return False
        for entry in (e, a):
            annotations = written.get(field_of(entry, UTRANNOTATOR_CONSEQUENCE_FIELDS))
            if annotations is None:
                return False
            if with_annotation and annotation_groups(field_of(entry, UTRANNOTATOR_ANNOTATION_FIELDS)) not in annotations:
                return False
        return True


def utrannotator_precedence_row(e: Entry, a: Entry) -> bool:
    """Whether the two entries' ``5UTR_consequence`` values are both among the plugin's five
    terms and differ (``utrannotator_consequence_precedence``)."""
    ev = next((e.fields[f] for f in UTRANNOTATOR_CONSEQUENCE_FIELDS if f in e.fields), ABSENT)
    av = next((a.fields[f] for f in UTRANNOTATOR_CONSEQUENCE_FIELDS if f in a.fields), ABSENT)
    return ev != av and ev in UTRANNOTATOR_TERMS and av in UTRANNOTATOR_TERMS


def structural_allele_text_row(e: Entry, a: Entry) -> bool:
    """Whether the reference's entry ``e`` carries the effect of a structural allele's ALT
    text read as deleted bases and the candidate's entry ``a`` none
    (``utrannotator_structural_allele_text``): the row's Allele is a structural class word on
    both sides, the reference's ``5UTR_consequence`` is one of the plugin's five terms, and
    the candidate fills neither ``5UTR_consequence`` nor ``5UTR_annotation``."""
    if e.allele not in STRUCTURAL_CLASS_WORDS or a.allele != e.allele:
        return False
    if field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS) not in UTRANNOTATOR_TERMS:
        return False
    return all(blank(field_of(a, names)) for names in (UTRANNOTATOR_CONSEQUENCE_FIELDS, UTRANNOTATOR_ANNOTATION_FIELDS))


def structural_deletion_span_row(e: Entry, a: Entry) -> bool:
    """Whether the candidate's entry ``a`` carries the effect of a symbolic deletion's stated
    span on a row where that span lies inside one exon of the 5' UTR
    (``utrannotator_structural_allele_text``). The row tells it by its own columns: its
    Allele is a deletion's class word on both sides; its Location's span is the distance
    between its two cDNA positions, so both ends lie in one exon; its Consequence names
    ``5_prime_UTR_variant`` beside no term but ``feature_truncation`` and
    ``NMD_transcript_variant``, so neither end leaves the 5' UTR. The candidate's
    ``5UTR_consequence`` is a term a deletion can fire, the frameshift only where the span is
    no multiple of three (UTRAnnotator.pm:923), beside a ``5UTR_annotation``; the reference's
    ``5UTR_consequence`` is one of the plugin's five terms, or both of its fields are empty.
    A VCF or JSON entry names no Location and is not the class's."""
    if e.allele not in STRUCTURAL_DELETION_WORDS or a.allele != e.allele or e.location is None:
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    positions = cdna_span(e)
    if not (start.isdigit() and (end.isdigit() or not end)) or positions is None:
        return False
    span = int(end or start) - int(start) + 1
    if span < 1 or positions[1] - positions[0] + 1 != span:
        return False
    if set(e.consequence.split(",")) - STRUCTURAL_DELETION_ROW_TERMS != {"5_prime_UTR_variant"}:
        return False
    reference_term = field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS)
    if blank(reference_term):
        if not blank(field_of(e, UTRANNOTATOR_ANNOTATION_FIELDS)):
            return False
    elif reference_term not in UTRANNOTATOR_TERMS:
        return False
    term = field_of(a, UTRANNOTATOR_CONSEQUENCE_FIELDS)
    if term not in UTRANNOTATOR_DELETION_TERMS or blank(field_of(a, UTRANNOTATOR_ANNOTATION_FIELDS)):
        return False
    return span % 3 != 0 or term != UTRANNOTATOR_FRAMESHIFT_TERM


def structural_insertion_length_row(e: Entry, a: Entry) -> bool:
    """Whether the candidate's entry ``a`` carries the frameshift a plain symbolic insertion's
    stated length decides (``utrannotator_structural_allele_text``). The row's Allele is
    ``insertion`` on both sides, its Consequence names ``5_prime_UTR_variant`` and its
    Location's span, the stated length, is no multiple of three (UTRAnnotator.pm:889, :923).
    The candidate's ``5UTR_consequence`` is the frameshift, and each effect of its
    ``5UTR_annotation`` holds exactly the six keys read from the reference UTR, or those six
    and ``alt_type=uORF``; the reference's ``5UTR_consequence`` is one of the plugin's five
    terms, or both of its fields are empty. The row's columns do not show where the candidate
    fills it, so the shape is exact by the value it admits. A VCF or JSON entry names no
    Location and is not the class's."""
    if e.allele != STRUCTURAL_INSERTION_WORD or a.allele != e.allele or e.location is None:
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    if not (start.isdigit() and (end.isdigit() or not end)):
        return False
    span = int(end or start) - int(start) + 1
    if span < 1 or span % 3 == 0 or "5_prime_UTR_variant" not in e.consequence.split(","):
        return False
    reference_term = field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS)
    if blank(reference_term):
        if not blank(field_of(e, UTRANNOTATOR_ANNOTATION_FIELDS)):
            return False
    elif reference_term not in UTRANNOTATOR_TERMS:
        return False
    annotation = field_of(a, UTRANNOTATOR_ANNOTATION_FIELDS)
    if field_of(a, UTRANNOTATOR_CONSEQUENCE_FIELDS) != UTRANNOTATOR_FRAMESHIFT_TERM or blank(annotation):
        return False
    for effect in annotation.split(","):
        pairs = [tuple(pair.partition("=")[::2]) for pair in effect.split(":")]
        keys = {key for key, _ in pairs}
        if len(keys) != len(pairs) or keys - {UTRANNOTATOR_ALTERNATE_TYPE[0]} != UTRANNOTATOR_REFERENCE_UTR_KEYS:
            return False
        if len(keys) > len(UTRANNOTATOR_REFERENCE_UTR_KEYS) and UTRANNOTATOR_ALTERNATE_TYPE not in pairs:
            return False
    return True


def lost_start_after_insertion_row(e: Entry, a: Entry, record_rows: Collection[Entry] | None) -> bool:
    """Whether the reference's entry reports as lost a start codon that an insertion only moves
    on and the candidate's entry, which carries the plugin's counts, reports no lost start
    (``utrannotator_lost_start_after_insertion``). An insertion's index in the 5' UTR is that
    of the base after it (UTRAnnotator.pm:1285-1313), which is the row's lower cDNA position
    where the two cDNA positions are the insertion's flanking bases. A codon reported at that
    very index begins directly after the insertion (``:769``) and stands whole whatever the
    allele. A codon reported one index below it holds the insertion after its first base and
    stands whole further on when the allele, read on the transcript's strand, ends with that
    base, ``A``; one reported two below holds it after its second base and stands whole when
    the allele so read ends with ``AT`` (``allele_on_transcript_strand``; a row that states no
    strand is refused there). Where the row's record states its alleles (``Entry.insertion``,
    else ``record_allele_lengths``) the reference allele is empty; a text row under an
    identifier with no record handed is read by its positions and its allele alone."""
    if field_of(a, UTRANNOTATOR_CONSEQUENCE_FIELDS) == UTRANNOTATOR_LOST_START_TERM or blank(field_of(a, UTRANNOTATOR_COUNT_FIELDS)):
        return False
    lost = lost_start_position(e)
    span = cdna_span(e)
    if lost is None or span is None or span[1] - span[0] != 1:
        return False
    if re.fullmatch(r"[ACGTN]+", e.allele) is None:
        return False
    # How many of the codon's bases stand before the insertion.
    before = span[0] - lost
    if before:
        read = allele_on_transcript_strand(e)
        if before not in (1, 2) or read is None or not read.endswith(START_CODON[:before]):
            return False
    if e.insertion is not None:
        return e.insertion
    lengths = record_allele_lengths(e, record_rows)
    return lengths is None or lengths[0] == 0


def gained_start_remade_in_place_row(
    e: Entry, a: Entry, record_rows: Collection[Entry] | None, utrannotator_sweep: UtrannotatorSeedSweep | None
) -> bool:
    """Whether the reference reports as gained, on an insertion's row, a start codon the
    inserted bases leave at its own index, and the candidate's entry, which carries the
    plugin's counts, reports no gained start there
    (``utrannotator_gained_start_remade_in_place``). The row's lower cDNA position is the 5'
    UTR index of the base after the insertion (UTRAnnotator.pm:1285-1313). A start codon
    reported one index below it holds the insertion after its first base, and the mutated UTR
    reads a start codon at that index when the allele, read on the transcript's strand
    (``allele_on_transcript_strand``; a row that states no strand is refused), opens with the
    codon's last two bases, ``TG``; one reported two below holds it after its second base
    and the allele so read opens with ``G``. The index is the one the reference's own
    ``5UTR_annotation`` names beside the gained-start term; where the reference writes another
    of the plugin's terms, the one index the runs of ``utrannotator_sweep`` write beside the
    gained start on the row. A gained start at a negative index, which the reference's
    window yields by wrapping before the UTR's first base, is no codon of the UTR and is
    refused. Where the row's record states its alleles (``Entry.insertion``, else
    ``record_allele_lengths``) the reference allele is empty."""
    if blank(field_of(a, UTRANNOTATOR_COUNT_FIELDS)):
        return False
    span = cdna_span(e)
    read = allele_on_transcript_strand(e)
    if span is None or span[1] - span[0] != 1 or read is None:
        return False
    term = field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS)
    if term == UTRANNOTATOR_GAINED_START_TERM:
        at = gained_start_position(e)
    elif term in UTRANNOTATOR_TERMS and utrannotator_sweep is not None:
        at = utrannotator_sweep.gained_start(e)
    else:
        at = None
    # A negative index is the wrapped window's (``utrannotator_gained_start_window``), no codon's.
    if at is None or at < 0:
        return False
    # How many of the codon's bases stand before the insertion.
    before = span[0] - at
    if before not in (1, 2) or not read.startswith(START_CODON[before:]):
        return False
    if field_of(a, UTRANNOTATOR_CONSEQUENCE_FIELDS) == UTRANNOTATOR_GAINED_START_TERM and gained_start_position(a) in (at, None):
        return False
    if e.insertion is not None:
        return e.insertion
    lengths = record_allele_lengths(e, record_rows)
    return lengths is None or lengths[0] == 0


def annotation_frames(value: str) -> list[dict[str, str]] | None:
    """The groups of a ``5UTR_annotation`` as maps of key to value, one group per upstream
    ORF: the text formats' ``key=value`` pairs joined by ``:``, with ``,`` between groups
    (``&`` in a VCF value), or the numbered objects of the JSON format's object. ``None`` for
    a value that is neither, an empty one included."""
    if value.startswith("{"):
        try:
            parsed = json.loads(value)
        except ValueError:
            return None
        if not isinstance(parsed, dict) or not all(isinstance(group, dict) for group in parsed.values()):
            return None
        return [{str(key): str(item) for key, item in group.items()} for group in parsed.values()]
    if blank(value) or ("," in value and "&" in value):
        return None
    groups = [group.split(":") for group in value.split("&" if "&" in value else ",")]
    if not all("=" in pair for group in groups for pair in group):
        return None
    return [dict(pair.partition("=")[::2] for pair in group) for group in groups]


def utr_base_after_insertion(e: Entry) -> int | None:
    """The 5' UTR base directly after the insertion of the entry's row as its record writes
    it, in bases before the coding start, read from the row's HGVSc and ``HGVS_OFFSET``. An
    ``ins`` between two neighbouring 5' UTR positions names that base as its second
    coordinate; a ``dup`` stands for the copy written directly after the duplicated bases, so
    the base is the one after its last coordinate. The notation stands 3' of the record's own
    insertion point on the transcript by the size of ``HGVS_OFFSET``, which is printed
    negative on a reverse-strand transcript (OutputFactory.pm:1755), so the record's point
    lies that many bases further from the coding start. ``None`` for any other notation and
    for a coordinate outside the 5' UTR."""
    notation = field_of(e, HGVSC_FIELDS).partition(":")[2]
    offset = field_of(e, HGVS_OFFSET_FIELDS)
    shift = abs(int(offset)) if re.fullmatch(r"-?\d+", offset) else 0
    inserted = UTR_HGVSC_INSERTION_RE.fullmatch(notation)
    if inserted is not None:
        first, last = int(inserted[1]), int(inserted[2])
        return -last + shift if last == first + 1 else None
    duplicated = UTR_HGVSC_DUPLICATION_RE.fullmatch(notation)
    if duplicated is not None:
        return -int(duplicated[2] or duplicated[1]) - 1 + shift
    return None


def frameshift_at_stop_boundary_row(e: Entry, a: Entry) -> bool:
    """Whether both entries report the frameshift and the candidate's ``5UTR_annotation``
    lists, beyond every group of the reference's, only upstream ORFs whose stop codon begins
    directly after the insertion (``utrannotator_frameshift_at_stop_boundary``).
    UTRAnnotator.pm:921 holds the record's own index against the check point one index
    before an ORF's stop codon (``:911-917``), so an insertion written directly before the
    stop codon's first base is past it and the ORF is left out, though it runs on past its
    stop codon's index and its length changes. Such a group is a uORF whose two printed
    lengths differ and whose stop codon's first base, ``ref_StartDistanceToCDS`` less
    ``ref_type_length`` plus three bases before the coding start, is the base
    ``utr_base_after_insertion`` reads from the reference's row; a row that states no
    insertion point is refused. No group both entries list may be such an ORF: there the
    reference's own test admits what it leaves out here."""
    terms = {field_of(entry, UTRANNOTATOR_CONSEQUENCE_FIELDS) for entry in (e, a)}
    if terms != {UTRANNOTATOR_FRAMESHIFT_TERM}:
        return False
    expected = annotation_frames(field_of(e, UTRANNOTATOR_ANNOTATION_FIELDS))
    actual = annotation_frames(field_of(a, UTRANNOTATOR_ANNOTATION_FIELDS))
    if expected is None or actual is None:
        return False
    listed = Counter(tuple(sorted(group.items())) for group in expected)
    written = Counter(tuple(sorted(group.items())) for group in actual)
    added = written - listed
    if listed - written or not added:
        return False
    after = utr_base_after_insertion(e)
    if after is None:
        return False

    def ends_at_the_insertion(pairs: tuple[tuple[str, str], ...]) -> bool:
        group = dict(pairs)
        start, before = group.get("ref_StartDistanceToCDS", ABSENT), group.get("ref_type_length", ABSENT)
        if group.get("ref_type") != "uORF" or not (start.isdigit() and before.isdigit()):
            return False
        return group.get("alt_type_length") != before and int(start) - int(before) + 3 == after

    return all(ends_at_the_insertion(pairs) for pairs in added) and not any(
        ends_at_the_insertion(pairs) for pairs in written & listed
    )


def lost_stop_written_before_stop_codon_row(e: Entry, a: Entry, record_rows: Collection[Entry] | None) -> bool:
    """Whether the reference reports the frameshift and the candidate, on an insertion's row,
    the lost stop of upstream ORFs the reference lists among the shifted ones
    (``utrannotator_lost_stop_written_before_stop_codon``). UTRAnnotator.pm:628-633 read a
    stop codon's index only for a record written inside the stop codon, so for the same
    mutated sequence written before the start codon that overlaps it the ORF is reported
    shifted (``:921``), with its new length. Each of the candidate's groups carries the
    ``KozakContext`` and ``KozakStrength`` of exactly one uORF of the reference's list, a
    different one for each group, and the alternative stop codon where that ORF's printed
    ``alt_type_length`` ends it: at ``ref_StartDistanceToCDS`` plus the allele's length less
    ``alt_type_length`` plus three bases before the coding start where the ORF still closes
    in the 5' UTR (``alt_type=uORF``), none otherwise. The row's cDNA positions are an
    insertion's two flanking bases, the candidate's entry carries the plugin's counts, and
    where the row's record states its alleles (``Entry.insertion``, else
    ``record_allele_lengths``) the reference allele is empty."""
    if field_of(e, UTRANNOTATOR_CONSEQUENCE_FIELDS) != UTRANNOTATOR_FRAMESHIFT_TERM:
        return False
    if field_of(a, UTRANNOTATOR_CONSEQUENCE_FIELDS) != UTRANNOTATOR_LOST_STOP_TERM or blank(field_of(a, UTRANNOTATOR_COUNT_FIELDS)):
        return False
    span = cdna_span(e)
    if span is None or span[1] - span[0] != 1 or re.fullmatch(r"[ACGTN]+", e.allele) is None:
        return False
    shifted = annotation_frames(field_of(e, UTRANNOTATOR_ANNOTATION_FIELDS))
    lost = annotation_frames(field_of(a, UTRANNOTATOR_ANNOTATION_FIELDS))
    if shifted is None or lost is None:
        return False
    taken: set[int] = set()
    for group in lost:
        kozak = (group.get("KozakContext"), group.get("KozakStrength"))
        frames = [
            at for at, frame in enumerate(shifted)
            if frame.get("ref_type") == "uORF" and (frame.get("KozakContext"), frame.get("KozakStrength")) == kozak
        ]
        if None in kozak or len(frames) != 1 or frames[0] in taken:
            return False
        taken.add(frames[0])
        frame = shifted[frames[0]]
        start, length = frame.get("ref_StartDistanceToCDS", ABSENT), frame.get("alt_type_length", ABSENT)
        if frame.get("alt_type") == "uORF":
            if not (start.isdigit() and length.isdigit()):
                return False
            stop = (group.get("AltStop"), group.get("AltStopDistanceToCDS"))
            if stop != ("True", str(int(start) + len(e.allele) - int(length) + 3)):
                return False
        elif (group.get("AltStop"), group.get("AltStopDistanceToCDS")) != ("False", "NA"):
            return False
    if e.insertion is not None:
        return e.insertion
    lengths = record_allele_lengths(e, record_rows)
    return lengths is None or lengths[0] == 0


def cdna_pair_spans(body: str) -> list[int | None]:
    """The distance between the two coordinates of each ``start_end`` pair of an HGVSc body, in
    positions, where the pair states it: two plain positions on one side of the stop codon
    (coding numbering has no position 0, so 1 follows -1), or two offsets from one exon
    boundary. ``None`` for a pair that mixes a plain position with an offset, names two
    boundaries or runs over the stop codon."""
    spans: list[int | None] = []
    for start_star, start, start_offset, end_star, end, end_offset in CDNA_COORDINATE_PAIR_RE.findall(body):
        if start_star != end_star:
            spans.append(None)
        elif start_offset and end_offset:
            spans.append(int(end_offset) - int(start_offset) if start == end else None)
        elif start_offset or end_offset:
            spans.append(None)
        else:
            first, last = (int(token) + (int(token) < 0) for token in (start, end))
            spans.append(last - first)
    return spans


class TranscriptSymbols:
    """The gene symbol of each transcript of a candidate JSON cache (``--json-cache``), by
    stable id: the symbol SpliceAI.pm:314 looks a transcript's score record up under. A
    chromosome's shards are read on the first lookup of a row whose Location names its
    directory, every shard on the first lookup of any other entry."""

    def __init__(self, cache: Path) -> None:
        self.root = cache / "transcripts"
        if not self.root.is_dir():
            raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
        self._symbols: dict[str | None, dict[str, str]] = {}

    def _read(self, chromosome: str | None) -> dict[str, str]:
        directories = [self.root / chromosome] if chromosome is not None else sorted(self.root.iterdir())
        symbols: dict[str, str] = {}
        for directory in directories:
            for shard in sorted(directory.glob("*.json")) + sorted(directory.glob("*.json.gz")):
                try:
                    parsed = json.loads(read_text(shard))
                except ValueError as exc:
                    raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
                for t in NearestSites.shard_transcripts(shard, parsed):
                    if t.get("stable_id") and t.get("gene_symbol"):
                        symbols[t["stable_id"]] = t["gene_symbol"]
        return symbols

    def symbol(self, e: Entry) -> str | None:
        """The gene symbol of the entry's transcript, or ``None`` for a transcript the cache
        does not hold or gives no symbol; a Feature carrying a version the cache's stable id
        lacks (``--transcript_version``) is looked up without it."""
        location = at_its_location(e).location
        chromosome = location.split(":", 1)[0] if location else None
        if chromosome is not None and not (self.root / chromosome).is_dir():
            chromosome = None
        if chromosome not in self._symbols:
            self._symbols[chromosome] = self._read(chromosome)
        symbols = self._symbols[chromosome]
        return symbols.get(e.feature) or symbols.get(e.feature.split(".", 1)[0])


class SpliceAiTable:
    """The score files a scored run's SpliceAI plugin read (``--spliceai-table``, once per
    file: the ``snv=`` file and the ``indel=`` file), each a VCF with its tabix index beside
    it (SpliceAI.pm:114-122), read through that index at the positions a row names and never
    whole. ``RevelTable``
    reads the rows of a tabix-indexed table at a span of positions, and a VCF's eight columns
    read as its eight-column layout: the position, REF and ALT cells stand in their own
    places and the INFO column in the score's. ``json_cache`` is the run's candidate JSON
    cache, whose transcripts name the gene of a row that states no ``SYMBOL``
    (``TranscriptSymbols``). ``refusals`` counts the rows ``spliceai_multi_entry_cell`` is
    asked about under the files and those it refuses (``SpliceAiRefusals``)."""

    def __init__(self, paths: Collection[Path], json_cache: Path | None = None) -> None:
        self.paths = list(paths)
        self.symbols = TranscriptSymbols(json_cache) if json_cache is not None else None
        self.refusals = SpliceAiRefusals()
        self._files: list[RevelTable] = []
        for path in self.paths:
            try:
                table = RevelTable(path)
            except ComparatorError as exc:
                raise self._refusal(exc) from exc
            if table._columns[0] != 3:
                raise ComparatorError(f"--spliceai-table: {path} is not a VCF (a header line of eight columns)")
            self._files.append(table)

    @staticmethod
    def _refusal(exc: ComparatorError) -> ComparatorError:
        """A refusal of the table reader under this option's name."""
        reason = str(exc).partition(": ")[2].replace(
            "is not a REVEL table (a header of 7, 8 or 9 columns)", "is not a VCF (a header line of eight columns)"
        )
        return ComparatorError(f"--spliceai-table: {reason}")

    def records(self, chromosome: str, low: int, high: int) -> list[tuple[int, str, str, str]]:
        """The POS, REF, ALT and INFO cells of every record at the positions ``low`` to
        ``high`` of ``chromosome``, the files in the order given and each in file order; a
        line that lacks one of the four cells is left out."""
        found: list[tuple[int, str, str, str]] = []
        for table in self._files:
            try:
                rows = table.rows(chromosome, low, high)
            except ComparatorError as exc:
                raise self._refusal(exc) from exc
            found += [
                (int(row.position), row.reference, row.alternate, row.score)
                for row in rows
                if row.position.isdigit() and row.reference and row.alternate and row.score
            ]
        return found


class SpliceAiRefusals:
    """The rows ``spliceai_multi_entry_cell`` is asked about under the score files and, of
    them, the rows it refuses, each under the first ground met for it (``SPLICEAI_GROUNDS``).
    A row is its entry's key, Location, allele and feature, and counts once however many of
    its fields the class is asked about."""

    def __init__(self) -> None:
        self._asked: set[tuple[str, str | None, str, str]] = set()
        self._refused: dict[tuple[str, str | None, str, str], str] = {}

    def note(self, e: Entry, ground: str | None) -> None:
        """Counts the row of ``e`` as asked about and, under a ``ground``, as refused."""
        row = (e.key, e.location, e.allele, e.feature)
        self._asked.add(row)
        if ground is not None:
            self._refused.setdefault(row, ground)

    def counts(self) -> dict:
        """The rows ``asked`` about, the rows ``refused`` and, under ``grounds``, the refused
        rows of each ground that refused any, in the order of ``SPLICEAI_GROUNDS``."""
        refused = Counter(self._refused.values())
        return {
            "asked": len(self._asked), "refused": len(self._refused),
            "grounds": {ground: refused[ground] for ground in SPLICEAI_GROUNDS if refused[ground]},
        }

    def clear(self) -> None:
        self._asked.clear()
        self._refused.clear()


def spliceai_refusal_lines(counts: dict, whose: str = "") -> list[str]:
    """The lines a summary prints of the counts of a ``SpliceAiRefusals``: the rows the class
    was asked about and refused, then the refused rows of each ground. ``whose`` names the
    rows' side of a summary of several."""
    head = f"{SPLICEAI_MULTI_ENTRY_CELL} under --spliceai-table{whose}: rows asked about {counts['asked']}, rows refused {counts['refused']}"
    return [head, *(f"{SPLICEAI_MULTI_ENTRY_CELL}{whose}: rows refused {rows}: {SPLICEAI_GROUNDS[ground]}" for ground, rows in counts["grounds"].items())]


class SpliceAiValue(NamedTuple):
    """A row's SpliceAI value: its ``SpliceAI_pred`` text where its format writes the one
    field, else ``None``; its nine members (``SPLICEAI_MEMBERS``); its cutoff flag, empty
    for none; and whether the members are JSON values, whose numbers compare as numbers."""

    whole: str | None
    members: tuple[str, ...]
    cutoff: str
    numbers: bool

    def writes(self, value: str) -> bool:
        """Whether the row writes ``value``, a ``|``-joined value as the default format
        prints it: as its whole text where the row carries the one field, else as its nine
        members, the parts of ``value`` past the ninth being those SpliceAI.pm:276-287 drop."""
        if self.whole is not None:
            return self.whole == value
        parts = value.split("|")[: len(SPLICEAI_MEMBERS)]
        return len(parts) == len(self.members) and all(
            member == part or (self.numbers and as_number(member) is not None and as_number(member) == as_number(part))
            for member, part in zip(self.members, parts)
        )


def spliceai_row_value(e: Entry | None) -> SpliceAiValue | None:
    """The SpliceAI value of a row, read from its ``SpliceAI_pred`` field, else from its nine
    ``SpliceAI_pred_`` fields (a CSQ field's ``&`` read as the comma it stands for), else from
    the members of its JSON ``spliceai`` object; ``None`` for a row that writes none."""
    if e is None:
        return None

    def written(name: str) -> str:
        value = field_of(e, (name,))
        return "" if blank(value) else value

    whole = written(SPLICEAI_VALUE_FIELD)
    if whole:
        return SpliceAiValue(whole, tuple(whole.split("|")[: len(SPLICEAI_MEMBERS)]), written(SPLICEAI_CUTOFF_FIELD), False)
    split = tuple(written(f"{SPLICEAI_VALUE_FIELD}_{member}").replace("&", ",") for member in SPLICEAI_MEMBERS)
    if any(split):
        return SpliceAiValue(None, split, written(SPLICEAI_CUTOFF_FIELD), False)
    members = tuple(json_object_member(e, "spliceai", member) or "" for member in SPLICEAI_MEMBERS)
    if any(members):
        return SpliceAiValue(None, members, json_object_member(e, "spliceai", SPLICEAI_CUTOFF_FIELD) or "", True)
    return None


def spliceai_entries_at(e: Entry, table: SpliceAiTable) -> tuple[list[list[list[str]]], str | None]:
    """The entries of the score records at a row's variant, one list per record in file order,
    each entry its ten ``|``-joined members: the allele, the symbol, four scores and four
    positions (SpliceAI.pm:276-287, 328-332). A record stands at the row's variant when its
    REF and ALT, spelt as Ensembl VEP
    spells a VCF line's (``vcf_alleles_as_spelt``), give the row's Location and Allele; a
    record whose ALT lists several alleles matches no allele (SpliceAI.pm:272). Beside them
    stands the ground on which they are not read (a key of ``SPLICEAI_GROUNDS``), ``None``
    where they are: the Location names no span, one of the records at the variant carries
    another INFO column than the one key or an entry of another length, the records there
    spell more than one REF, or no record stands there."""
    located = at_its_location(e)
    span = re.fullmatch(r"(.+):(\d+)(?:-(\d+))?", located.location or "")
    if span is None:
        return [], "no_span"
    low, high = sorted((int(span.group(2)), int(span.group(3) or span.group(2))))
    records: list[list[list[str]]] = []
    references: set[str] = set()
    for position, reference, alternate, info in table.records(span.group(1), low - 1, high):
        if "," in alternate:
            continue
        spelt, (allele,) = vcf_alleles_as_spelt(reference, [alternate])
        start = position + (spelt != reference)
        end = start + len(spelt.strip("-")) - 1
        if (min(start, end), max(start, end)) != (low, high) or allele != e.allele:
            continue
        entries = [entry.split("|") for entry in info.removeprefix(SPLICEAI_KEY).split(",")]
        if not info.startswith(SPLICEAI_KEY) or ";" in info or any(len(entry) != 1 + len(SPLICEAI_MEMBERS) for entry in entries):
            return [], "other_info"
        references.add(reference)
        records.append(entries)
    if len(references) > 1:
        return [], "several_references"
    return records, None if records else "no_record"


def spliceai_multi_entry_cell_ground(e: Entry, a: Entry | None, name: str, ev: str, av: str, table: SpliceAiTable) -> str | None:
    """The ground on which ``spliceai_multi_entry_cell`` refuses the differing SpliceAI field
    ``name`` of a row under the score files (a key of ``SPLICEAI_GROUNDS``), ``None`` where
    it accepts it. Every record line of ``table``
    at the row's variant is read, in file order (``spliceai_entries_at``). A line that holds
    an entry naming the row's gene holds several entries; the candidate's row writes the
    gene's entry of the last line that names it, the last such entry of that line, whether
    or not the lines give the gene entries that differ in text; and the reference's row writes
    the whole value, behind its first allele, of the last line whose first entry names the
    gene (the two cutoff flags then equal), SpliceAI.pm:305 keeping one line a first symbol,
    or no SpliceAI value and no flag where no line's first entry names it. The row's gene is
    the ``SYMBOL`` either row writes, else its transcript's symbol in the table's candidate
    cache. Nothing is accepted on a row whose gene neither names or on a row whose gene no
    entry names."""
    ours = spliceai_row_value(a)
    if ours is None:
        return "no_candidate_value"
    records, ground = spliceai_entries_at(e, table)
    if ground is not None:
        return ground
    stated = next((value for value in (field_of(e, ("SYMBOL", "gene_symbol")), field_of(a, ("SYMBOL", "gene_symbol"))) if not blank(value)), None)
    symbol = stated if stated is not None or table.symbols is None else table.symbols.symbol(e)
    if symbol is None:
        return "no_gene"
    holding = [record for record in records if any(entry[1] == symbol for entry in record)]
    own = [entry for record in holding for entry in record if entry[1] == symbol]
    if not own:
        return "gene_not_named"
    if all(len(record) < 2 for record in holding):
        return "lone_entry"
    if not ours.writes("|".join(own[-1][1:])):
        return "candidate_value"
    theirs = spliceai_row_value(e)
    flag = field_named(name, (SPLICEAI_CUTOFF_FIELD,))
    led = [record for record in holding if record[0][1] == symbol]
    if not led:
        silent = theirs is None and blank(field_of(e, (SPLICEAI_CUTOFF_FIELD,)))
        taken = silent and blank(ev) and (av in SPLICEAI_CUTOFF_FLAGS if flag else not blank(av))
        return None if taken else "reference_not_silent"
    whole = ",".join("|".join(entry) for entry in led[-1]).partition("|")[2]
    if theirs is None or not theirs.writes(whole):
        return "reference_value"
    if flag or theirs.cutoff != ours.cutoff:
        return "reference_flag"
    return None


def spliceai_multi_entry_cell_row(e: Entry, a: Entry | None, name: str, ev: str, av: str, table: SpliceAiTable | None) -> bool:
    """Whether the differing SpliceAI field ``name`` of a row is the difference a score record
    of several entries makes (``spliceai_multi_entry_cell``): nothing without the table, and
    under it the field ``spliceai_multi_entry_cell_ground`` refuses on no ground. The row
    is counted among the table's ``refusals``, under the ground where one is met."""
    if table is None:
        return False
    ground = spliceai_multi_entry_cell_ground(e, a, name, ev, av, table)
    table.refusals.note(e, ground)
    return ground is None


def field_class_accepts(
    klass: str, e: Entry, name: str, expected_value: str, actual_value: str, a: Entry | None = None,
    record_rows: Collection[Entry] | None = None,
    gwas_pairs: GwasCatalogPairs | None = None,
    utrannotator_sweep: UtrannotatorSeedSweep | None = None,
    codon_tables: OrganelleCodonTables | None = None,
    downstream_sequences: DownstreamSequences | None = None,
    transcript_edges: TranscriptEdges | None = None,
    candidate_rows: Collection[Entry] | None = None,
    declared: Collection[str] = (),
    *,
    intron_ends: IntronEnds | None = None,
    revel_table: RevelTable | None = None,
    spliceai_table: SpliceAiTable | None = None,
    assembly: str | None = None,
    splice_windows: SpliceWindows | None = None,
) -> bool:
    """Whether the built-in field class ``klass`` accepts the differing field ``name`` of the
    entry ``e``: the field is one the class reads (``FIELD_CLASSES``) and the two values have
    the class's shape; ``a`` is the candidate's entry, read by a class whose shape spans two
    fields of the row, and ``record_rows`` the reference's entries of ``e``'s own input record
    (every allele's rows, at least those on its feature), which the caller vouches are one
    record's: a class that reads a record's alleles accepts nothing without them.
    ``gwas_pairs`` is the catalog the two GWAS classes read and accept nothing without.
    ``revel_table`` is the table ``revel_other_reference_amino_acid_rows`` reads and accepts
    nothing without.
    ``spliceai_table`` is the score files ``spliceai_multi_entry_cell`` reads and accepts
    nothing without; it counts the rows the class is asked about and refuses under them.
    ``assembly`` is the run's assembly as its manifest names it, which
    ``plugin_contig_mitochondrial_spelling`` reads: without one the class takes the rows it
    takes on another assembly than GRCh37.
    ``codon_tables`` is the candidate cache's codon tables, without which
    ``organelle_peptide_paths_translation_table`` accepts nothing.
    ``downstream_sequences`` is the candidate cache's transcripts, without which the classes
    of ``TRANSCRIPT_READING_CLASSES`` accept nothing.
    ``transcript_edges`` is the candidate cache's exon spans, without which
    ``hgvs_clipped_reference_against_whole_alternate`` is its open-side test alone.
    ``candidate_rows`` is the candidate's entries of the same input record, without which the
    three sibling-allele classes accept nothing.
    ``utrannotator_sweep`` is the seeded runs ``utrannotator_consequence_precedence`` reads and
    accepts nothing without, and ``utrannotator_gained_start_remade_in_place`` reads on a row
    whose reference writes another term than the gained start. ``intron_ends`` is the candidate cache's intron ends, which
    ``hgvsp_intron_end_insertion`` and ``hgvsp_intron_start_insertion`` read and accept nothing
    without.
    ``splice_windows`` is the candidate cache's introns, which
    ``maxentscan_overlapping_windows`` reads and accepts nothing without.
    ``declared`` is the classes named with ``klass``: a Downstream value two causes change is
    accepted only where it names the second cause's class too (``downstream_class``).
    A class of ``TERM_PAIR_FOLLOWING_CLASSES`` accepts a cell of a row its own predicate
    takes, where ``term_pair_following_cell`` reads the two values from the two sets of terms;
    it reads ``a`` and accepts nothing without it."""
    if not field_named(name, FIELD_CLASSES[klass]):
        return False
    ev = "-" if blank(expected_value) else expected_value
    av = "-" if blank(actual_value) else actual_value
    if klass in TERM_PAIR_FOLLOWING_CLASSES:
        return a is not None and TERM_PAIR_CLASSES[klass](e, a) and term_pair_following_cell(e, a, name, ev, av)
    if klass == GXF_REVERSE_STRAND_UTR_SEQUENCE:
        if blank(ev) or blank(av) or not any(e.fields.get(f) == "-1" for f in STRAND_FIELDS):
            return False
        expected_match, actual_match = TER_COUNT_RE.match(ev), TER_COUNT_RE.match(av)
        return expected_match is not None and actual_match is not None and expected_match.group(1) == actual_match.group(1)
    if klass == ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE:
        return (
            codon_tables is not None and a is not None and organelle_peptide_path(e, a, ev, av)
            and codon_tables.table(e) is not None
        )
    if klass == CONSERVATION_MAX_TEXT_ORDER:
        reference, candidate = as_number(ev), as_number(av)
        return reference is not None and candidate is not None and reference < candidate < 0
    if klass in (GWAS_HAPLOTYPE_CROSS_PRODUCT, GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE):
        site = gwas_site(e, record_rows) if gwas_pairs is not None else None
        if site is None:
            return False
    if klass == GWAS_INTERACTION_FIRST_SNP_ALLELE_WHITESPACE:
        if not blank(ev) or gwas_annotated(e):
            return False
        # The plugin writes a field the catalog row leaves empty as a defined empty string
        # (``GWAS_odds_ratio=``), which the reference's unannotated row has no key for: the
        # class's beside the candidate's filled fields, as the cross-product class reads it.
        if blank(av) and not (actual_value == "" and a is not None and gwas_annotated(a)):
            return False
        return gwas_pair_named(e, field_of(a, GWAS_STUDY_FIELDS), site, gwas_pairs.candidate_only)
    if klass == GWAS_HAPLOTYPE_CROSS_PRODUCT:
        if not gwas_pair_named(e, field_of(e, GWAS_STUDY_FIELDS), site, gwas_pairs.reference_only):
            return False
        # The plugin writes a field the catalog row leaves empty as a defined empty string
        # (``GWAS_beta_coef=``), which a row the candidate never annotated has no key for; such
        # a pair is the class's when the reference row carries any ``GWAS_`` value. The empty
        # string is the reference's: a key the candidate alone writes stays a difference.
        if not blank(ev) or blank(av):
            return not blank(ev) or (expected_value == "" and blank(actual_value) and gwas_annotated(e))
        # A value the candidate alone writes is its own catalog row's where it reports another
        # study than the reference, whose row leaves the field empty.
        own_study = field_of(a, GWAS_STUDY_FIELDS)
        return gwas_annotated(e) and not blank(own_study) and own_study != field_of(e, GWAS_STUDY_FIELDS)
    if klass == SPLICEVAULT_SIBLING_ALLELE_ROWS:
        if blank(ev) or record_rows is None:
            return False
        return sibling_allele_value(e, FIELD_CLASSES[klass], record_rows, same_amino_acids=False, candidate_rows=candidate_rows)
    if klass in (ALPHAMISSENSE_TWO_DIGIT_TRANSCRIPT_VERSION, SPLICEVAULT_JSON_LATER_ROWS):
        return blank(ev) and not blank(av)
    if klass == CADD_STRUCTURAL_ALLELE_MATCH:
        return blank(ev) and not blank(av) and e.allele != "-" and not set(e.allele.upper()) <= set("ACGTN")
    if klass == UTRANNOTATOR_ANNOTATION_PAIR_ORDER:
        return reordered_pair_groups(ev, av)
    if klass == NEAREST_EXON_JB_ITEM_ORDER:
        return reordered_items(ev, av)
    if klass == UTRANNOTATOR_GAINED_START_WINDOW:
        return a is not None and gained_start_window_row(e, a, record_rows)
    if klass == UTRANNOTATOR_CONSEQUENCE_PRECEDENCE:
        if a is None or utrannotator_sweep is None or not utrannotator_precedence_row(e, a):
            return False
        # An annotation the two entries share is no differing field of its own, so the
        # consequence is accepted only where the runs write it beside both terms.
        annotations = [field_of(entry, UTRANNOTATOR_ANNOTATION_FIELDS) for entry in (e, a)]
        shared = ABSENT not in annotations and annotation_groups(annotations[0]) == annotation_groups(annotations[1])
        return utrannotator_sweep.writes(e, a, shared or field_named(name, UTRANNOTATOR_ANNOTATION_FIELDS))
    if klass == UTRANNOTATOR_STRUCTURAL_ALLELE_TEXT:
        if a is None:
            return False
        if blank(av):
            return not blank(ev) and structural_allele_text_row(e, a)
        if e.allele == STRUCTURAL_INSERTION_WORD:
            return structural_insertion_length_row(e, a)
        return structural_deletion_span_row(e, a)
    if klass == UTRANNOTATOR_LOST_START_AFTER_INSERTION:
        return a is not None and lost_start_after_insertion_row(e, a, record_rows)
    if klass == UTRANNOTATOR_GAINED_START_REMADE_IN_PLACE:
        return a is not None and gained_start_remade_in_place_row(e, a, record_rows, utrannotator_sweep)
    if klass == UTRANNOTATOR_FRAMESHIFT_AT_STOP_BOUNDARY:
        return a is not None and frameshift_at_stop_boundary_row(e, a)
    if klass == UTRANNOTATOR_LOST_STOP_WRITTEN_BEFORE_STOP_CODON:
        return a is not None and lost_stop_written_before_stop_codon_row(e, a, record_rows)
    if klass == LOEUF_ZERO_MINIMUM:
        return blank(ev) and av == "0"
    if klass == MAXENTSCAN_OVERLAPPING_WINDOWS:
        if splice_windows is None or not splice_windows.holds(e):
            return False
        if field_named(name, MAXENTSCAN_WINDOW_TEXT_FIELDS):
            return MAXENTSCAN_DONOR_TEXT_RE.fullmatch(ev) is not None and MAXENTSCAN_ACCEPTOR_TEXT_RE.fullmatch(av) is not None
        return as_number(ev) is not None and as_number(av) is not None
    if klass in (ALPHAMISSENSE_SIBLING_ALLELE_ROWS, REVEL_SIBLING_ALLELE_ROWS):
        return not blank(ev) and record_rows is not None and sibling_allele_value(
            e, FIELD_CLASSES[klass], record_rows, candidate_rows=candidate_rows
        )
    if klass == REVEL_OTHER_REFERENCE_AMINO_ACID_ROWS:
        return revel_other_reference_amino_acid_row(e, a, ev, av, record_rows, revel_table)
    if klass == GNOMADMT_INFO_VALUE_SECOND_EQUALS:
        kept, sign, _ = decode_extra(av).partition("=")
        return not blank(ev) and bool(sign) and decode_extra(ev) == kept
    if klass == LOFTEE_ANC_ALLELE_MULTIALLELIC_RECORD:
        if record_rows is None or lof_filter_item_side(e, a, ANC_ALLELE_FILTER) != "candidate":
            return False
        if not LOF_SUBSTITUTION_TERMS & set(consequence_set(e.consequence).split(",")):
            return False
        if re.fullmatch(r"[ACGT]", e.allele) is None or record_indel_row(e, record_rows):
            return False
        return bool(sibling_rows(e, record_rows))
    if klass == LOFTEE_FIRST_ALT_DONOR_FILTER:
        if record_rows is None or lof_filter_item_side(e, a, GC_TO_GT_DONOR_FILTER) is None:
            return False
        if LOF_DONOR_TERM not in consequence_set(e.consequence).split(",") or not sibling_rows(e, record_rows):
            return False
        # The first alternate's own row reads its own base in both engines, and rows that state
        # no allele order do not tell it from a later allele's.
        first = first_alternate(e, record_rows)
        return first is not None and first != e.allele
    if klass == LOFTEE_GERP_ZOOM_LEVEL_SUMMARY:
        totals = gerp_totals(e, a)
        if totals is None or not LOF_TRUNCATING_TERMS & set(consequence_set(e.consequence).split(",")):
            return False
        reference_total, candidate_total, items = totals
        walked = as_number(items.get("BP_DIST", ""))
        # ``BP_DIST`` is the interval's bases less one: 79 walks the 80-base interval.
        if walked is None or walked + 1 < GERP_EXACT_BELOW_BASES:
            return False
        if field_named(name, LOF_INFO_FIELDS):
            return True
        side = lof_filter_item_side(e, a, END_TRUNC_FILTER)
        reach = as_number(items.get("DIST_FROM_LAST_EXON", ""))
        if side is None or reach is None or reach > LAST_EXON_REACH:
            return False
        lower, higher = (reference_total, candidate_total) if side == "reference" else (candidate_total, reference_total)
        return lower <= GERP_END_TRUNC_CUTOFF < higher
    if klass == SPLICE_FAMILY_SWAP:
        if a is None or blank(ev) or not blank(av):
            return False
        # Downstream.pm:95-99 write on a frameshift row none of whose terms names a splice site.
        terms = consequence_set(e.consequence).split(",")
        if FRAMESHIFT_TERM not in terms or any("splice" in term for term in terms):
            return False
        pair = registry_class(consequence_set(e.consequence), consequence_set(a.consequence), e.allele, (SPLICE_FAMILY_SWAP,))
        return pair == SPLICE_FAMILY_SWAP
    if klass == NEAREST_EXON_JB_STRETCH_PRIMING:
        distance = field_of(e, DISTANCE_FIELDS)
        return blank(ev) and not blank(av) and distance.isdigit() and int(distance) <= FRAMESHIFT_INTRON_STRETCH
    if klass == NEAREST_EXON_JB_SHIFTED_SPAN:
        if not positive_shift(e) or feature_kind(entry_kind(e)) != "transcript":
            return False
        if blank(ev) or blank(av):
            return blank(ev) != blank(av)
        expected_items, actual_items = nearest_exon_items(ev), nearest_exon_items(av)
        if expected_items is None or actual_items is None:
            return False
        return {exon_id for exon_id, *_ in expected_items} != {exon_id for exon_id, *_ in actual_items}
    if klass == NEAREST_EXON_JB_ID_LESS_EXON_FOLD:
        expected_items, actual_items = nearest_exon_items(ev), nearest_exon_items(av)
        if expected_items is None or actual_items is None or len(expected_items) != 1:
            return False
        if any(exon_id for exon_id, *_ in (*expected_items, *actual_items)):
            return False
        _, folded_distance, folded_side, _ = expected_items[0]
        if (folded_distance, folded_side) not in {(distance, side) for _, distance, side, _ in actual_items}:
            return False
        # One candidate item beside another length is the fold of two overlapped exons: the row
        # begins in one exon and ends in another, across both splice sites between them.
        terms = set(consequence_set(e.consequence).split(","))
        return len(actual_items) > 1 or (SPLICE_SITE_TERMS <= terms and cdna_span_closed(e))
    if klass == DOWNSTREAM_REVERSE_STRAND_ROTATION:
        if blank(ev) or blank(av) or FRAMESHIFT_TERM not in consequence_set(e.consequence).split(","):
            return False
        if not any(e.fields.get(f) == "-1" for f in STRAND_FIELDS):
            return False
        offset = field_of(e, HGVS_OFFSET_FIELDS)
        length = notation_indel_length(field_of(e, HGVSC_FIELDS))
        if INTEGER_RE.fullmatch(offset) is None or length is None or length <= 0:
            return False
        bases_shifted = abs(int(offset))
        return bases_shifted > length and bases_shifted % length != 0
    if klass == SPLICEAI_CHOPPED_MULTIALLELIC:
        if record_rows is None or not blank(ev) or blank(av):
            return False
        if re.fullmatch(r"[ACGTN]+", e.allele) is None or record_indel_row(e, record_rows):
            return False
        return "-" in record_alleles(e, record_rows)
    if klass == SPLICEAI_MULTI_ENTRY_CELL:
        return spliceai_multi_entry_cell_row(e, a, name, ev, av, spliceai_table)
    if klass == PLUGIN_CONTIG_MITOCHONDRIAL_SPELLING:
        if not blank(ev) or blank(av):
            return False
        if field_named(name, (GNOMAD_PLUGIN_PREFIX,)) and not name.lower().endswith(GNOMADC_FIELD_SUFFIX):
            return False
        on_grch37 = assembly is not None and assembly.strip().lower() in GRCH37_NAMES
        if on_grch37 and not field_named(name, ALPHAMISSENSE_PREFIXES):
            return False
        return record_chromosome(e, record_rows) in MITOCHONDRIAL_SPELLINGS
    if klass == REFSEQ_EDITED_MODEL_HGVSC_OFFSET:
        if blank(ev) or blank(av) or not any(e.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS):
            return False
        expected_name, _, expected_body = ev.partition(":")
        actual_name, _, actual_body = av.partition(":")
        if expected_name != actual_name or CDNA_COORDINATE_RE.sub("#", expected_body) != CDNA_COORDINATE_RE.sub("#", actual_body):
            return False
        spans = zip(cdna_pair_spans(expected_body), cdna_pair_spans(actual_body))
        return all(expected is None or actual is None or expected == actual for expected, actual in spans)
    if klass == REFSEQ_EDITED_MODEL_SHIFTED_INSERTION_BASES:
        return edited_model_shifted_insertion_bases(e, a, ev, av)
    if klass == REFSEQ_EDITED_MODEL_SEAM_INSERTION_HGVSC:
        return edited_model_seam_insertion_replacement(e, a, ev, av)
    if klass == REFSEQ_EDITED_MODEL_UNMINIMISED_ALLELE_HGVSC:
        return edited_model_clipped_alleles(e, a, ev, av)
    if klass == REFSEQ_EDITED_MODEL_HGVS_OFFSET:
        if not blank(ev) or not any(e.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS):
            return False
        hgvsc = next((e.fields[f] for f in HGVSC_FIELDS if f in e.fields), ABSENT)
        coding = any(not blank(e.fields.get(f, ABSENT)) for f in CDS_POSITION_FIELDS)
        return coding and hgvsc.endswith("dup") and NON_ZERO_INTEGER_RE.fullmatch(av) is not None
    if klass == REFSEQ_EDITED_MODEL_DELETED_SPAN_CODONS:
        return a is not None and deleted_span_codons_pair(e, a, ev, av)
    if klass == HGVS_INTRON_OFFSET_ABUTTING_INSERTION:
        return intron_offset_abutting_insertion(e, a, name, ev, av)
    if klass == HGVS_SHIFTED_INSERTION_CLAMPED_END:
        offset = next((e.fields[f] for f in HGVS_OFFSET_FIELDS if f in e.fields), ABSENT)
        if blank(ev) or blank(av) or blank(offset) or not insertion_row(e) or not positive_shift(e):
            return False
        # The shift carried the insertion from before the first base to after base k, so its
        # cDNA position starts at the shift's length.
        shift = next((e.fields[f] for f in SHIFT_LENGTH_FIELDS if f in e.fields), ABSENT)
        cdna_start = next((e.fields[f] for f in CDNA_START_FIELDS if f in e.fields), ABSENT).partition("-")[0]
        if cdna_start != shift:
            return False
        expected_name, _, expected_body = ev.partition(":")
        actual_name, _, actual_body = av.partition(":")
        clamped = clamped_insertion_position(expected_body)
        shifted = shifted_insertion_end(actual_body)
        if expected_name != actual_name or clamped is None or shifted is None:
            return False
        clamped_numbering, clamped_position, clamped_bases, stripped_count = clamped
        shifted_numbering, shifted_position, shifted_bases, shifted_count = shifted
        if (clamped_numbering, clamped_position) != (shifted_numbering, shifted_position):
            return False
        if shifted_count != len(e.allele) or len(clamped_bases) + stripped_count != len(e.allele):
            return False
        return shifted_bases is None or shifted_bases.startswith(clamped_bases)
    if klass == HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE:
        return shifted_indel_written_span_gate(e, a, name, ev, av)
    if klass == HGVSP_SHIFTED_INDEL_WRITTEN_SPAN_GATE:
        return shifted_indel_written_span_protein_gate(e, a, ev, av)
    if klass == HGVS_SHIFTED_DELETION_CLAMPED_SPAN:
        if field_named(name, HGVS_INTRON_START_OFFSET_FIELDS):
            return a is not None and shifted_deletion_start_offsets(e, a, ev, av)
        return shifted_deletion_clamped_span(e, ev, av)
    if klass == HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE:
        if field_named(name, HGVS_INTRON_OFFSET_FIELDS):
            return clipped_notation_intron_offsets(e, a, name, ev, av, transcript_edges)
        return clipped_reference_against_whole_alternate(e, ev, av, transcript_edges)
    if klass == HGVS_ADDED_BASES_LEFT_AT_THE_SPAN_END:
        return added_bases_left_at_the_span_end(e, a, name)
    if klass == HGVS_INSERTION_LEFT_WHERE_THE_SHIFT_STOPS:
        return insertion_left_where_the_shift_stops(e, a, name)
    if klass == HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT:
        return added_bases_written_outside_the_transcript(e, a, name, ev, av)
    if klass == HGVSP_STOP_CODON_PROTEIN_COMPARISON:
        return stop_codon_protein_comparison(e.consequence, field_of(e, FLAGS_FIELDS), ev, av)
    if klass == HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE:
        return stop_named_on_incomplete_coding_sequence(e.consequence, field_of(e, FLAGS_FIELDS), ev, av)
    if klass == HGVSP_INTRON_END_INSERTION:
        return hgvsp_intron_end_insertion(e, a, ev, av, intron_ends)
    if klass == HGVSP_INTRON_START_INSERTION:
        return hgvsp_intron_start_insertion(e, a, ev, av, intron_ends)
    if klass == SHIFT_INTRON_TERM_UNSHIFTED_SPAN:
        if blank(ev) or blank(av) or not positive_shift(e):
            return False
        expected_terms = set(consequence_set(ev).split(","))
        actual_terms = set(consequence_set(av).split(","))
        if INTRON_TERM in actual_terms or expected_terms != actual_terms | {INTRON_TERM}:
            return False
        hgvsc = next((e.fields[f] for f in HGVSC_FIELDS if f in e.fields), ABSENT)
        if blank(hgvsc):
            return insertion_row(e) and SPLICE_ACCEPTOR_TERM in actual_terms
        offsets = acceptor_site_offsets(hgvsc.partition(":")[2])
        if offsets is None:
            return False
        if offsets:
            return True
        # A notation without an intron offset is numbered over the exons alone, so a deletion
        # over a whole intron reads as exonic too; such a span covers both splice sites, and the
        # row carries intron_variant from its shifted span on both engines.
        return not actual_terms & SPLICE_SITE_TERMS
    if klass == FREQUENCY_FILTER_ZERO_VALUE:
        return freqs_differ_by_zero_read_as_na(ev, av)
    if klass == HGNC_ID_COPY:
        return hgnc_id_copy_pair(e, a, ev, av)
    if klass == HGVS_SHIFT_SHARED_ACROSS_ALLELES:
        if a is None or downstream_sequences is None:
            return False
        return later_allele_shared_shift(e, a, name, record_rows, downstream_sequences)
    if klass == BLOSUM62_HGVS_SHIFT_SHARED_ACROSS_ALLELES:
        if a is None or downstream_sequences is None:
            return False
        return blosum62_later_allele_shared_shift(e, a, ev, av, record_rows, downstream_sequences)
    if klass in TRANSCRIPT_READING_CLASSES:
        if a is None or downstream_sequences is None:
            return False
        if field_named(name, NMD_FIELDS):
            return nmd_reads_shifted_pair(e, a, expected_value, actual_value, downstream_sequences, record_rows)
        # The values as the entries carry them: an empty one is the plugin's, a ``-`` or an
        # absent one is no value.
        decided = downstream_class(e, a, name, expected_value, actual_value, downstream_sequences, record_rows)
        return decided is not None and decided[0] == klass and (decided[1] is None or decided[1] in declared)
    if klass == PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN:
        if a is None or not splice_term_past_the_record_span_pair(e, a):
            return False
        return term_pair_following_cell(e, a, name, ev, av)
    if klass == PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR:
        if a is None or not any(pair(e, a) for pair in STOP_CODON_PLUGIN_VALUE_PAIRS):
            return False
        return term_pair_following_cell(e, a, name, ev, av)
    return False


def peptide_description_shape(description: str) -> str:
    """``dup``, ``ins`` or ``other`` for the text after ``p.``; a ``delins`` is ``other``."""
    if description.endswith("dup"):
        return "dup"
    if re.search(r"(?<!del)ins", description):
        return "ins"
    return "other"


def organelle_peptide_path(e: Entry, a: Entry, expected_value: str, actual_value: str) -> bool:
    """Whether a differing HGVSp pair names one protein and differs in a path of
    TranscriptVariationAllele.pm that translates a sequence bare, on a row the bare translation
    alone can explain: ``_get_fs_peptides`` and ``_stop_loss_extra_AA`` write the whole
    description of a frameshift (its residues, an immediate stop, the stop distance and the
    lost-start shapes), ``_stop_loss_extra_AA`` the count of an extension after a lost stop, and
    ``_check_for_peptide_duplication`` the choice between a duplication and an insertion. A
    substitution, a ``delins`` and a plain deletion read their residues through ``peptide``,
    which follows the transcript's table, and a value one side leaves empty names no protein,
    so each stays a difference. The class changes the HGVSp alone: the reference entry ``e``
    and the candidate's ``a`` state one coding change, the same consequence set and the same
    ``CODING_CHANGE_FIELDS``, so a row whose terms, positions, codons, coding notation or shift
    differ (a later allele annotated at another allele's shifted span, an indel shifted short at
    a sequence edge) is left to the class of that difference.
    An inserted allele of two or more bases on a reverse-strand transcript, shifted past its
    own length by a length that is no multiple of it, stays a difference too:
    TranscriptVariationAllele.pm:1337-1339 rotate a reverse-strand transcript's allele by its
    length less the shift, nothing once the shift is the longer, so the HGVSp is read from
    bases the record does not produce. A row that does not state a forward strand is read as
    a reverse-strand one."""
    expected_name, _, expected_description = expected_value.rpartition(":p.")
    actual_name, _, actual_description = actual_value.rpartition(":p.")
    if not expected_name or expected_name != actual_name:
        return False
    def coding_change(entry: Entry) -> dict[str, str]:
        return {k.lower(): v for k, v in entry.fields.items() if field_named(k, CODING_CHANGE_FIELDS) and not blank(v)}

    if consequence_set(e.consequence) != consequence_set(a.consequence) or coding_change(e) != coding_change(a):
        return False
    offset = field_of(e, HGVS_OFFSET_FIELDS)
    if len(e.allele) >= 2 and re.fullmatch(r"[ACGTN]+", e.allele) and INTEGER_RE.fullmatch(offset):
        shift = abs(int(offset))
        if shift > len(e.allele) and shift % len(e.allele) and field_of(e, STRAND_FIELDS) != "1":
            return False
    terms = set(consequence_set(e.consequence).split(","))
    if "frameshift_variant" in terms:
        return True
    if "stop_lost" in terms:
        return "ext" in expected_description and "ext" in actual_description
    shapes = {peptide_description_shape(expected_description), peptide_description_shape(actual_description)}
    return bool(terms & {"inframe_insertion", "protein_altering_variant"}) and shapes == {"dup", "ins"}


def hgnc_id_copy_pair(e: Entry, a: Entry | None, expected_value: str, actual_value: str) -> bool:
    """Whether the HGNC_ID ``expected_value`` of the reference entry ``e`` and ``actual_value``
    of the candidate entry ``a`` have the shape of ``hgnc_id_copy``: a Transcript entry, one
    value empty and the other an HGNC id, and the SYMBOL_SOURCE either side prints a source
    other than HGNC. A gene HGNC names carries its own id on both engines, so a differing id
    there, and two differing ids anywhere, stay differences."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript":
        return False
    if blank(expected_value) == blank(actual_value):
        return False
    filled = actual_value if blank(expected_value) else expected_value
    if HGNC_ID_RE.fullmatch(filled) is None:
        return False
    sources = {field_of(entry, SYMBOL_SOURCE_FIELDS) for entry in (e, a) if entry is not None}
    sources.discard(ABSENT)
    return bool(sources) and all(not blank(source) and source != HGNC_SOURCE for source in sources)


def codon_table_id(value: object) -> int | None:
    """The genetic code a cache transcript's ``codon_table`` names: a number in a native cache,
    the dump's string in a converted one, table 1 for 0 or no key (Downstream.pm:132); ``None``
    for a table ``GENETIC_CODES`` lacks and for any other value."""
    if value is None:
        return 1
    text = str(value)
    if not text.isdecimal():
        return None
    table = int(text) or 1
    return table if table in GENETIC_CODES else None


def cache_integer(value: object) -> int | None:
    """An integer of a cache transcript, a number in a native cache and the dump's string in a
    converted one; ``None`` for any other value."""
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and INTEGER_RE.fullmatch(value):
        return int(value)
    return None


@dataclass(frozen=True)
class TranscriptModel:
    """A cache transcript's coordinates: its ``strand``, its ``exons`` in transcript order as
    (start, end), its ``cdna_coding_start`` and ``cdna_coding_end``, the count of ``N`` its
    ``translateable_seq`` opens with (the start exon's phase, which the CDS coordinates count
    in) and its ``mapper`` pairs as (cDNA start, cDNA end, genomic start, genomic end,
    orientation)."""

    strand: int
    exons: tuple[tuple[int, int], ...]
    coding_start: int
    coding_end: int
    padding: int
    pairs: tuple[tuple[int, int, int, int, int], ...]

    @classmethod
    def from_cache(cls, transcript: dict, padding: int) -> TranscriptModel | None:
        """The model of a cache transcript object, ``None`` where the cache states a part of
        it in no form this reads."""
        mapper = transcript["variation_effect_feature_cache"].get("mapper")
        pairs = mapper.get("pairs") if isinstance(mapper, dict) else None
        exons = transcript.get("exons")
        strand = cache_integer(transcript.get("strand"))
        coding_start = cache_integer(transcript.get("cdna_coding_start"))
        coding_end = cache_integer(transcript.get("cdna_coding_end"))
        if strand not in (1, -1) or coding_start is None or coding_end is None:
            return None
        if not isinstance(pairs, list) or not isinstance(exons, list) or not pairs or not exons:
            return None
        if not all(isinstance(item, dict) for item in (*pairs, *exons)):
            return None
        read_pairs = [tuple(cache_integer(pair.get(key)) for key in MAPPER_PAIR_KEYS) for pair in pairs]
        read_exons = [(cache_integer(exon.get("start")), cache_integer(exon.get("end"))) for exon in exons]
        if any(None in item for item in (*read_pairs, *read_exons)):
            return None
        return cls(strand, tuple(sorted(read_exons, reverse=strand == -1)), coding_start, coding_end, padding, tuple(read_pairs))

    def cdna(self, genomic: int) -> int | None:
        """The cDNA coordinate of a genomic base, ``None`` for a base outside the mapper pairs."""
        for cdna_start, _, genomic_start, genomic_end, orientation in self.pairs:
            if genomic_start <= genomic <= genomic_end:
                return cdna_start + (genomic_end - genomic if orientation == -1 else genomic - genomic_start)
        return None

    def genomic(self, cds: int) -> int | None:
        """The genomic base of a CDS coordinate, ``None`` for a coordinate outside the coding
        bounds or the mapper pairs."""
        cdna = cds - self.padding + self.coding_start - 1
        if not self.coding_start <= cdna <= self.coding_end:
            return None
        for cdna_start, cdna_end, genomic_start, genomic_end, orientation in self.pairs:
            if cdna_start <= cdna <= cdna_end:
                return genomic_end - (cdna - cdna_start) if orientation == -1 else genomic_start + cdna - cdna_start
        return None

    def cds_pair(self, start: int, end: int) -> tuple[int | None, int | None]:
        """The ``cds_start`` and ``cds_end`` of the genomic span ``start`` to ``end``
        (BaseTranscriptVariation.pm:252-264 over ``TranscriptMapper::genomic2cds``), ``None``
        for an end that is no coding base: the CDS coordinates of the span's 5' and 3' bases;
        for an insertion (``start`` is ``end`` + 1) ``map_insert``'s result
        (Mapper.pm:504-558), the base after the 5' flank's coordinate and the base before the
        3' flank's, a flank outside the mapper pairs read as the base beside the other."""
        to_cds = self.padding - self.coding_start + 1
        if start == end + 1:
            five, three = (end, start) if self.strand == 1 else (start, end)
            before, after = self.cdna(five), self.cdna(three)
            if before is None and after is None:
                return None, None
            before = after - 1 if before is None else before
            after = before + 1 if after is None else after
            return (
                before + 1 + to_cds if self.coding_start <= before < self.coding_end else None,
                after - 1 + to_cds if self.coding_start < after <= self.coding_end else None,
            )
        five, three = (start, end) if self.strand == 1 else (end, start)
        first, last = (self.cdna(base) for base in (five, three))
        return (
            first + to_cds if first is not None and self.coding_start <= first <= self.coding_end else None,
            last + to_cds if last is not None and self.coding_start <= last <= self.coding_end else None,
        )

    def genomic_span(self, cds_start: int, cds_end: int) -> tuple[int, int] | None:
        """The one genomic span (start, end; an insertion's start one past its end) whose
        ``cds_pair`` is ``cds_start`` and ``cds_end``, or ``None``: an end the mapper pairs do
        not hold, or an insertion's pair (``cds_start`` is ``cds_end`` + 1) between two bases
        the genome does not keep adjacent, which the point after the one and the point before
        the other both give."""
        first, last = self.genomic(min(cds_start, cds_end)), self.genomic(max(cds_start, cds_end))
        if first is None or last is None or (cds_start > cds_end and abs(first - last) != 1):
            return None
        return (max(first, last), min(first, last)) if cds_start > cds_end else (min(first, last), max(first, last))

    def ends_where_nmd_escapes(self, genomic_end: int) -> bool:
        """``variant_exon_check`` (NMD.pm:132-176) on a record's genomic end: it lies in the
        last exon, or in the second to last within 51 bases of that exon's 3' end, its end on
        the forward strand and its start on the reverse."""
        start, end = self.exons[-1]
        if start <= genomic_end <= end:
            return True
        if len(self.exons) < 2:
            return False
        start, end = self.exons[-2]
        return start <= genomic_end <= start + 51 if self.strand == -1 else end - 51 <= genomic_end <= end


@dataclass(frozen=True)
class DownstreamTranscript:
    """What the Downstream and NMD plugins read of one cache transcript: the
    ``translateable_seq``, the ``three_prime_utr`` (empty where the cache spells none), the
    codon table (``codon_table_id``, ``None`` for a table ``GENETIC_CODES`` lacks) and the
    length of the ``peptide`` (``None`` without one) of its
    ``variation_effect_feature_cache``, and its ``TranscriptModel`` (``None`` where the cache
    does not state one)."""

    coding: str
    utr: str
    table: int | None
    peptide_length: int | None
    model: TranscriptModel | None

    def downstream_value(self, start: int, end: int, allele: str, length_change: bool, ensembl: bool) -> str | None:
        """The Downstream plugin's value on the CDS pair ``start`` and ``end`` (an insertion's
        start one past its end) with the bases ``allele`` on the transcript's strand: the
        DownstreamProtein, or with ``length_change`` the ProteinLengthChange as printed. Both
        rules write the allele over the pair (Downstream.pm:108-110) and translate the changed
        coding sequence with the 3' UTR up to the first stop; the ProteinLengthChange is the
        whole translation's length less the peptide's (:141-148; the translated coding
        sequence's for a transcript the cache gives no peptide). With ``ensembl`` the
        DownstreamProtein starts at the pair's lower coordinate rounded down to a multiple of
        three (:112-115) and there is no value for a start past the coding sequence's end,
        where the assignment of :110 fails; without it, the candidate's rule, at the codon
        holding the first changed base."""
        if ensembl and start - 1 > len(self.coding):
            return None
        changed = self.coding[: start - 1] + allele + self.coding[max(end, start - 1) :]
        if length_change:
            reference = len(translate_until_stop(self.coding, self.table)) if self.peptide_length is None else self.peptide_length
            return str(len(translate_until_stop(changed + self.utr, self.table)) - reference)
        first = min(start, end) if ensembl else min(start - 1, len(self.coding))
        return translate_until_stop(changed[first - first % 3 :] + self.utr, self.table)


class DownstreamSequences:
    """The transcripts of a candidate JSON cache as the Downstream and NMD plugins read them
    (``DownstreamTranscript``), by transcript stable id, kept for a transcript with a
    ``translateable_seq``. A chromosome's shards are read on the first lookup of an entry
    whose Location names its directory, every shard on the first lookup of an entry without
    one (a VCF or JSON entry). A stable id the shards read give one reading under two models
    (a transcript the cache files under two chromosomes, each with its own mapper pairs) is
    kept without a model; one they give two different readings otherwise is left out."""

    def __init__(self, cache: Path) -> None:
        self.root = cache / "transcripts"
        if not self.root.is_dir():
            raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
        self._chromosomes = {directory.name for directory in self.root.iterdir() if directory.is_dir()}
        self._sequences: dict[str | None, dict[str, DownstreamTranscript]] = {}

    def _read(self, chromosome: str | None) -> dict[str, DownstreamTranscript]:
        sequences: dict[str, DownstreamTranscript] = {}
        conflicting: set[str] = set()
        for name in sorted(self._chromosomes) if chromosome is None else [chromosome]:
            directory = self.root / name
            for shard in sorted(directory.glob("*.json")) + sorted(directory.glob("*.json.gz")):
                try:
                    parsed = json.loads(read_text(shard))
                except ValueError as exc:
                    raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
                for t in NearestSites.shard_transcripts(shard, parsed):
                    stable_id, vefc = t.get("stable_id"), t.get("variation_effect_feature_cache")
                    if not isinstance(stable_id, str) or not isinstance(vefc, dict):
                        continue
                    coding, utr = vefc.get("translateable_seq"), vefc.get("three_prime_utr")
                    if not isinstance(coding, str) or not coding or not isinstance(utr, (str, type(None))):
                        continue
                    coding, peptide = coding.upper(), vefc.get("peptide")
                    read = DownstreamTranscript(
                        coding, (utr or "").upper(), codon_table_id(vefc.get("codon_table")),
                        len(peptide) if isinstance(peptide, str) else None,
                        TranscriptModel.from_cache(t, len(coding) - len(coding.lstrip("N"))),
                    )
                    held = sequences.setdefault(stable_id, read)
                    if held != read:
                        if replace(held, model=None) == replace(read, model=None):
                            sequences[stable_id] = replace(read, model=None)
                        else:
                            conflicting.add(stable_id)
        for stable_id in conflicting:
            del sequences[stable_id]
        return sequences

    def transcript(self, e: Entry) -> DownstreamTranscript | None:
        """The entry's transcript, or ``None`` for one the cache does not hold with a
        ``translateable_seq``; a Feature carrying a version the cache's stable id lacks is
        looked up without it."""
        chromosome = e.location.rpartition(":")[0] if e.location else None
        if chromosome not in self._chromosomes:
            chromosome = None
        if chromosome not in self._sequences:
            self._sequences[chromosome] = self._read(chromosome)
        sequences = self._sequences[chromosome]
        return sequences.get(e.feature, sequences.get(re.sub(r"\.\d+$", "", e.feature)))


def translate_until_stop(sequence: str, table: int) -> str:
    """The translation of ``sequence`` under the genetic code ``table`` (one of
    ``GENETIC_CODES``) up to its first stop, as ``Bio::Tools::CodonTable::translate`` reads it
    for Downstream.pm:139-142: the whole codons in order, a trailing one or two bases dropped;
    a codon with an ambiguous letter by the amino acids its readings code under the table, one
    amino acid as itself, aspartate with asparagine ``B``, glutamate with glutamine ``Z`` and
    anything else ``X``, as is a codon with a letter that is no nucleotide's."""
    code = GENETIC_CODES[table]
    residues: list[str] = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i : i + 3]
        residue = code.get(codon)
        if residue is None:
            readings = [NUCLEOTIDE_READINGS.get(base, "") for base in codon]
            coded = {code[x + y + z] for x in readings[0] for y in readings[1] for z in readings[2]}
            residue = coded.pop() if len(coded) == 1 else AMBIGUOUS_RESIDUE_PAIRS.get(frozenset(coded), "X")
        if residue == "*":
            break
        residues.append(residue)
    return "".join(residues)


def cds_pair_readings(e: Entry) -> list[tuple[int | None, int | None]]:
    """The ``cds_start`` and ``cds_end`` an entry's coding span may spell, a gap (``?``, or a
    JSON bound the entry lacks) ``None``, and no pair for an entry without a coding
    coordinate: the text formats' ``CDS_position`` or the JSON ``cds_start`` and ``cds_end``
    as given and, for a sequence allele at two adjacent numbers, swapped as well, the pair of
    an insertion, whose start is one past its end and which every format spells lower first
    (``Utils::format_coords``; the JSON entry's two keys likewise). Each bound is read up to
    the ``/total`` suffix ``--total_length`` appends."""
    position = e.fields.get("CDS_position")
    if position is None:
        bounds = [e.fields.get(key, ABSENT).partition("/")[0] for key in ("cds_start", "cds_end")]
    else:
        first, _, last = position.partition("/")[0].partition("-")
        bounds = [first, last or first]
    start, end = (int(bound) if bound.isdigit() else None for bound in bounds)
    if start is None and end is None:
        return []
    readings = [(start, end)]
    if start is not None and end == start + 1 and re.fullmatch(r"[ACGTNacgtn]+", e.allele):
        readings.append((end, start))
    return readings


def record_span_readings(e: Entry) -> list[tuple[int, int]]:
    """The genomic start and end a text entry's Location may spell: the span as printed and,
    for a sequence allele at two adjacent bases, the point between them as well, an insertion's
    span, whose start is one past its end."""
    first, _, last = (e.location or "").rpartition(":")[2].partition("-")
    last = last or first
    if not (first.isdigit() and last.isdigit()):
        return []
    start, end = int(first), int(last)
    readings = [(start, end)]
    if end == start + 1 and re.fullmatch(r"[ACGTNacgtn]+", e.allele):
        readings.append((end, start))
    return readings


def shifted_pair_readings(
    e: Entry, model: TranscriptModel, shift: int, own_pairs: list[tuple[int | None, int | None]]
) -> list[tuple[tuple[int | None, int | None], tuple[int, int] | None, int]]:
    """The readings of a row whose reference entry ``e`` carries the HGVS_OFFSET ``shift``,
    each the own CDS pair (one of ``own_pairs``), the pair a plugin reads under --hgvs where
    the shift reaches it (the record's span moved by ``shift`` has a CDS coordinate at both
    ends), else ``None``, and the record's genomic end. The record's span is a reading of the
    Location that maps to an own pair or, for an entry without a Location, the one span an own
    pair maps back to (``TranscriptModel.genomic_span``)."""
    if e.location:
        spans = [(span, model.cds_pair(*span)) for span in record_span_readings(e)]
        spans = [(span, own) for span, own in spans if own in own_pairs]
    else:
        spans = [(model.genomic_span(*own), own) for own in own_pairs if None not in own]
        spans = [(span, own) for span, own in spans if span is not None]
    readings = []
    for (start, end), own in spans:
        shifted = model.cds_pair(start + shift, end + shift)
        readings.append((own, None if None in shifted else shifted, end))
    return readings


def hgvs_rotated_allele(allele: str, shift: int, strand: int) -> str:
    """The bases a plugin reads under --hgvs for the allele as written, ``allele``, on a row
    the notation moved by ``shift``: ``shift_feature_seqs``
    (TranscriptVariationAllele.pm:1330-1345) turns the allele left one base as many times as
    the shift's length on a forward-strand transcript and as the allele's length less it on a
    reverse-strand one, and ``feature_seq`` reads the result on the transcript's strand."""
    turns = abs(shift) if strand == 1 else len(allele) - abs(shift)
    if allele and turns > 0:
        turns %= len(allele)
        allele = allele[turns:] + allele[:turns]
    return reverse_complement(allele) if strand == -1 else allele


def later_allele_prints_shifted_pair(
    e: Entry, a: Entry, model: TranscriptModel | None, shift: int, record_rows: Collection[Entry] | None
) -> bool:
    """Whether the two CDS pairs of a row are those of a later allele under --hgvs: ``e`` is
    the reference's row of a second or later alternate allele of its record
    (``later_alternate``), the candidate entry ``a`` prints the record's own pair, and the
    reference prints the record's span moved by its HGVS_OFFSET ``shift`` or by the offset
    of the last earlier allele's row that carries one (``earlier_allele_shift``): the pair an
    earlier allele's protein notation left on the TranscriptVariation the record's alleles
    share (OutputFactory.pm:1737 over BaseTranscriptVariation.pm:252-267), which is the
    row's own offset only where the two alleles move alike."""
    if model is None or not later_alternate(e, record_rows):
        return False
    moves = [move for move in dict.fromkeys((shift, earlier_allele_shift(e, record_rows))) if move]
    return any(printed_pair_is_moved(e, a, model, move) for move in moves)


def printed_pair_is_moved(e: Entry, a: Entry, model: TranscriptModel, move: int, own: tuple[int | None, int | None] | None = None) -> bool:
    """Whether the reference entry ``e`` prints the record's span moved by ``move`` bases
    (``shifted_pair_readings``), the span that of a CDS pair the candidate entry ``a`` prints,
    or of its pair ``own`` alone."""
    printed = cds_pair_readings(e)
    own_pairs = cds_pair_readings(a) if own is None else [own]
    return any(moved in printed for _, moved, _ in shifted_pair_readings(e, model, move, own_pairs))


def downstream_class(
    e: Entry, a: Entry, name: str, expected_value: str, actual_value: str, sequences: DownstreamSequences,
    record_rows: Collection[Entry] | None = None,
) -> tuple[str, str | None] | None:
    """The one class that accepts the Downstream field ``name`` (a DownstreamProtein or a
    ProteinLengthChange) whose values are ``expected_value`` on the reference entry ``e`` and
    ``actual_value`` on the candidate entry ``a``, with the class that must be declared beside
    it (``None`` for a row one cause changes), or ``None``. Both are values, the entries
    share one consequence set, or two that differ as the registry's
    ``start_cooccurrence_swap`` admits, with ``frameshift_variant`` and no splice term and
    one CDS span without a gap, and the cache holds the transcript under a table of
    ``GENETIC_CODES``. Two CDS spans are one on the row of a later allele of ``record_rows``
    whose reference entry prints the record's span moved by its HGVS_OFFSET or by an earlier
    allele's (``later_allele_prints_shifted_pair``): the own pair is then the candidate's. The
    candidate's value is its rule on the own pair (``DownstreamTranscript.downstream_value``)
    and the reference's is Downstream.pm:104-149 on the pair the plugin read: the pair and
    the allele the notation shifted where the reference's non-zero HGVS_OFFSET reaches the
    plugin (``shifted_pair_readings``, ``hgvs_rotated_allele``), whatever earlier pair the
    row prints; on a later allele's row without an offset, the printed pair an earlier
    allele's notation left (``earlier_allele_shift``), with the allele as written; else the
    own pair. Where
    Downstream.pm gives another value on the own pair the class is
    ``plugin_hgvs_shifted_cds_coordinates``: alone where the candidate's value is that
    own-pair value, the value Ensembl VEP writes without the flag, and beside the own pair's
    class (``downstream_start_rule_class``) where it is not, since the start rule changes the
    value as well. Where Downstream.pm gives one value on both pairs, a DownstreamProtein on
    a row with one ProteinLengthChange is the own pair's class."""
    if any(value in ("-", ABSENT) for value in (expected_value, actual_value)):
        return None
    length_change = field_named(name, PROTEIN_LENGTH_CHANGE_FIELDS)
    if length_change and not (INTEGER_RE.fullmatch(expected_value) and INTEGER_RE.fullmatch(actual_value)):
        return None
    terms, candidate_terms = consequence_set(e.consequence), consequence_set(a.consequence)
    # One consequence set, or the pair the registry sets aside as ``start_cooccurrence_swap``.
    start_pair = (START_COOCCURRENCE_SWAP,)
    if terms != candidate_terms and registry_class(terms, candidate_terms, e.allele, start_pair) is None:
        return None
    if FRAMESHIFT_TERM not in terms.split(",") or "splice" in terms:
        return None
    own_pairs = cds_pair_readings(a)
    two_pairs = own_pairs != cds_pair_readings(e)
    if not own_pairs or None in own_pairs[0]:
        return None
    allele = "" if e.allele == "-" else e.allele.upper()
    if e.allele != "-" and re.fullmatch(r"[ACGTN]+", allele) is None:
        return None
    offset = field_of(e, HGVS_OFFSET_FIELDS)
    shift = int(offset) if NON_ZERO_INTEGER_RE.fullmatch(offset) else 0
    start, end = own_pairs[0]
    # A later allele's row that prints another pair reads the pair at its own offset, or with
    # none the pair the last earlier allele's notation left.
    earlier = (earlier_allele_shift(e, record_rows) or 0) if two_pairs else 0
    read_at = shift or earlier
    # The cache is read for a row the notation shifted or whose span starts on a codon's
    # third base.
    if not read_at and (two_pairs or length_change or start % 3 or start > end):
        return None
    transcript = sequences.transcript(e)
    if transcript is None or transcript.table is None:
        return None
    model = transcript.model
    if read_at and model is None:
        return None
    if two_pairs and not later_allele_prints_shifted_pair(e, a, model, shift, record_rows):
        return None
    strand = model.strand if model is not None else -1 if field_of(e, STRAND_FIELDS) == "-1" else 1
    own_allele = reverse_complement(allele) if strand == -1 else allele
    readings = shifted_pair_readings(e, model, read_at, own_pairs) if read_at else [(own_pairs[0], None, 0)]
    one_length_change = not blank(field_of(e, PROTEIN_LENGTH_CHANGE_FIELDS)) and field_of(
        e, PROTEIN_LENGTH_CHANGE_FIELDS
    ) == field_of(a, PROTEIN_LENGTH_CHANGE_FIELDS)
    for (start, end), shifted, _ in readings:
        # A later allele's reading is the one whose pair, moved by the row's offset or by the
        # earlier allele's, the reference prints.
        moves = [move for move in dict.fromkeys((shift, earlier)) if move]
        if two_pairs and not any(printed_pair_is_moved(e, a, model, move, (start, end)) for move in moves):
            continue
        if actual_value != transcript.downstream_value(start, end, own_allele, length_change, ensembl=False):
            continue
        reference = unshifted = transcript.downstream_value(start, end, own_allele, length_change, ensembl=True)
        if shifted is not None:
            reference = transcript.downstream_value(
                *shifted, hgvs_rotated_allele(allele, shift, strand), length_change, ensembl=True
            )
        if expected_value != reference:
            continue
        start_rule = None if length_change else downstream_start_rule_class(transcript, start, end, own_allele)
        if reference != unshifted:
            if actual_value == unshifted:
                return PLUGIN_HGVS_SHIFTED_CDS_COORDINATES, None
            if start_rule is not None:
                return PLUGIN_HGVS_SHIFTED_CDS_COORDINATES, start_rule
        elif one_length_change and start_rule is not None:
            return start_rule, None
    return None


def downstream_start_rule_class(transcript: DownstreamTranscript, start: int, end: int, allele: str) -> str | None:
    """The class of ``DOWNSTREAM_CLASSES`` whose shape the own pair ``start`` to ``end`` has
    with the bases ``allele`` in its place, or ``None``: a span that starts on a codon's third
    base, ``downstream_start_past_own_codon`` where the offset of Downstream.pm:113 lies
    inside the changed sequence or at its end, ``downstream_utr_frame_past_coding_end`` where
    it lies past that end and the transcript has a 3' UTR."""
    if start % 3 or start > end or end > len(transcript.coding):
        return None
    if start <= start - 1 + len(allele) + len(transcript.coding) - end:
        return DOWNSTREAM_START_PAST_OWN_CODON
    return DOWNSTREAM_UTR_FRAME_PAST_CODING_END if transcript.utr else None


def field_absent(e: Entry, value: str) -> bool:
    """Whether ``value`` is how the entry's format spells a field the row does not carry: no
    Extra pair of a default-format row and no key of a JSON entry (``ABSENT``; ``-`` where
    the Extra text is compared pair by pair), the ``-`` of a tab cell, and the empty cell of a
    CSQ entry (no Location and no ``variant_allele``), which holds every field of its header.
    An empty value anywhere else is a value."""
    return value in ("-", ABSENT) or (value == "" and e.location is None and "variant_allele" not in e.fields)


def nmd_reads_shifted_pair(
    e: Entry, a: Entry, expected_value: str, actual_value: str, sequences: DownstreamSequences,
    record_rows: Collection[Entry] | None = None,
) -> bool:
    """Whether the NMD ``expected_value`` of the reference entry ``e`` and ``actual_value`` of
    the candidate entry ``a`` have the shape of ``plugin_hgvs_shifted_cds_coordinates``: one
    is ``NMD_escaping_variant`` and the other no value (``field_absent``), the entries share one consequence set,
    or two that differ as ``stop_retained_without_stop_codon_test_pair`` admits, each
    with a term of ``NMD_GATE_TERMS``, and one CDS span, or the two of a later allele of
    ``record_rows`` (``later_allele_prints_shifted_pair``), whose own pair is the
    candidate's and whose printed pair is the span moved by the row's offset or an earlier
    allele's, the reference's non-zero HGVS_OFFSET
    reaches the plugin (``shifted_pair_readings``), the reference carries the value exactly
    when the rule of NMD.pm:98-109 holds on the shifted pair and the candidate exactly when it
    holds on the own pair: both ends of the pair are coordinates and the CDS end is at most
    101 or the record's genomic end passes ``variant_exon_check``
    (``TranscriptModel.ends_where_nmd_escapes``)."""
    reference, candidate = expected_value == NMD_ESCAPING, actual_value == NMD_ESCAPING
    if reference == candidate or not (field_absent(a, actual_value) if reference else field_absent(e, expected_value)):
        return False
    terms, candidate_terms = consequence_set(e.consequence), consequence_set(a.consequence)
    # One consequence set, or the pair ``stop_retained_without_stop_codon_test_pair`` admits.
    if terms != candidate_terms and not stop_retained_without_stop_codon_test_pair(e, a):
        return False
    if not all(NMD_GATE_TERMS & set(side.split(",")) for side in (terms, candidate_terms)):
        return False
    offset = field_of(e, HGVS_OFFSET_FIELDS)
    own_pairs, printed = cds_pair_readings(a), cds_pair_readings(e)
    if NON_ZERO_INTEGER_RE.fullmatch(offset) is None or not own_pairs:
        return False
    transcript = sequences.transcript(e)
    if transcript is None or transcript.model is None:
        return False
    two_pairs = own_pairs != printed
    if two_pairs and not later_allele_prints_shifted_pair(e, a, transcript.model, int(offset), record_rows):
        return False
    moves = [move for move in dict.fromkeys((int(offset), earlier_allele_shift(e, record_rows))) if move]
    for own, shifted, genomic_end in shifted_pair_readings(e, transcript.model, int(offset), own_pairs):
        if shifted is None or (two_pairs and not any(printed_pair_is_moved(e, a, transcript.model, move, own) for move in moves)):
            continue
        escapes = transcript.model.ends_where_nmd_escapes(genomic_end)
        if (shifted[1] <= 101 or escapes) == reference and (None not in own and (own[1] <= 101 or escapes)) == candidate:
            return True
    return False


def protein_span(e: Entry) -> tuple[int, int] | None:
    """The first and last codon an entry's protein position names, the lower first: the
    ``Protein_position`` read up to a ``/total`` suffix, or the JSON ``protein_start`` and
    ``protein_end``; ``None`` where either is no number."""
    position = e.fields.get("Protein_position")
    if position is None:
        bounds = [e.fields.get(key, ABSENT).partition("/")[0] for key in ("protein_start", "protein_end")]
    else:
        first, _, last = position.partition("/")[0].partition("-")
        bounds = [first, last or first]
    if not all(bound.isdigit() for bound in bounds):
        return None
    return min(map(int, bounds)), max(map(int, bounds))


def codon_span(pair: tuple[int, int]) -> tuple[int, int]:
    """The first and last codon of a CDS pair, the lower first; an insertion's pair names the
    codons of its two flanks."""
    return (min(pair) + 2) // 3, (max(pair) + 2) // 3


def codon_residues(transcript: DownstreamTranscript, codons: tuple[int, int]) -> str | None:
    """The coding sequence's translation over its codons ``codons`` (first, last) under the
    transcript's table, a stop ``*``; ``None`` for a codon past the sequence's end or one the
    table does not read."""
    code = GENETIC_CODES[transcript.table]
    residues = [code.get(transcript.coding[3 * n - 3 : 3 * n]) for n in range(codons[0], codons[1] + 1)]
    return None if None in residues else "".join(residues)


def later_allele_shared_shift(
    e: Entry, a: Entry, name: str, record_rows: Collection[Entry] | None, sequences: DownstreamSequences
) -> bool:
    """``hgvs_shift_shared_across_alleles`` on the differing field ``name`` of a later allele's
    row: the reference entry ``e`` is a second or later alternate's whose first allele's row
    carries a non-zero HGVS_OFFSET (``first_allele_shift``), the candidate entry ``a``
    prints the record's own CDS pair and the reference that span moved by the first allele's
    offset (``shifted_pair_readings``). A CDS position is the class's there. A protein
    position is where each entry's is the codons of its pair (``codon_span``) and the two
    differ. An ``Amino_acids`` is on such a row where the two values share their alternate
    part and each reference part is the translation of its entry's codons
    (``codon_residues``)."""
    shift = first_allele_shift(e, record_rows)
    own_pairs, printed = cds_pair_readings(a), cds_pair_readings(e)
    if shift is None or not own_pairs or None in own_pairs[0] or own_pairs == printed:
        return False
    transcript = sequences.transcript(e)
    if transcript is None or transcript.model is None:
        return False
    for own, moved, _ in shifted_pair_readings(e, transcript.model, shift, own_pairs):
        if moved is None or moved == own or moved not in printed:
            continue
        if field_named(name, CDS_POSITION_FIELDS):
            return True
        own_codons, moved_codons = codon_span(own), codon_span(moved)
        if own_codons == moved_codons or protein_span(a) != own_codons or protein_span(e) != moved_codons:
            continue
        if field_named(name, PROTEIN_POSITION_FIELDS):
            return True
        if transcript.table is None:
            continue
        reference, _, changed = field_of(e, AMINO_ACIDS_FIELDS).partition("/")
        own_reference, _, own_changed = field_of(a, AMINO_ACIDS_FIELDS).partition("/")
        if not changed or changed != own_changed:
            continue
        if reference == codon_residues(transcript, moved_codons) and own_reference == codon_residues(transcript, own_codons):
            return True
    return False


def blosum62_later_allele_shared_shift(
    e: Entry, a: Entry, ev: str, av: str, record_rows: Collection[Entry] | None, sequences: DownstreamSequences
) -> bool:
    """``blosum62_hgvs_shift_shared_across_alleles`` on a differing BLOSUM62: the reference
    entry ``e`` writes a value and the candidate entry ``a`` none, on a later allele's row
    whose CDS position and protein position ``hgvs_shift_shared_across_alleles`` takes
    (``later_allele_shared_shift``: the candidate prints the record's own pair and the
    reference that span moved by the first allele's offset, each protein position the codons
    of its pair). The reference's ``Amino_acids`` is one residue over one (Blosum62.pm:110),
    the first the translation of the moved pair's one codon (``codon_residues``) and the
    second the alternate part of the candidate's own ``Amino_acids``, which fails that test;
    and the value is the matrix cell of the reference's pair (``BLOSUM62_MATRIX``)."""
    if blank(ev) or not blank(av) or not blank(field_of(a, BLOSUM62_FIELDS)):
        return False
    if not all(later_allele_shared_shift(e, a, fields[0], record_rows, sequences) for fields in (CDS_POSITION_FIELDS, PROTEIN_POSITION_FIELDS)):
        return False
    scored = PEPTIDE_PAIR_RE.fullmatch(field_of(e, AMINO_ACIDS_FIELDS))
    own = field_of(a, AMINO_ACIDS_FIELDS)
    if scored is None or PEPTIDE_PAIR_RE.fullmatch(own) is not None or own.partition("/")[2] != scored.group(2):
        return False
    transcript, shift = sequences.transcript(e), first_allele_shift(e, record_rows)
    if transcript is None or transcript.model is None or transcript.table is None or shift is None:
        return False
    printed = cds_pair_readings(e)
    moved_codons = {codon_span(moved) for _, moved, _ in shifted_pair_readings(e, transcript.model, shift, cds_pair_readings(a))
                    if moved is not None and moved in printed}
    if not any(codon_residues(transcript, codons) == scored.group(1) for codons in moved_codons):
        return False
    return ev == BLOSUM62_MATRIX.get(scored.groups())


def positive_shift(e: Entry) -> bool:
    """Whether the entry's ``SHIFT_LENGTH`` (the JSON ``shift_length``) is a positive integer:
    a row ``--shift_3prime 1`` moved, printed under ``--shift_length``. Digits ``int`` does
    not read are a malformed input (``integer_of``)."""
    value = next((e.fields[f] for f in SHIFT_LENGTH_FIELDS if f in e.fields), ABSENT)
    return value.isdigit() and integer_of(value, "SHIFT_LENGTH") > 0


def shift_field_consistent(e: Entry) -> bool:
    """Whether the entry's ``SHIFT_LENGTH``, where the run prints one, is positive: a row
    ``--shift_3prime 1`` moved; a run without the field says nothing about the shift."""
    return all(f not in e.fields for f in SHIFT_LENGTH_FIELDS) or positive_shift(e)


def field_of(e: Entry | None, patterns: Collection[str]) -> str:
    """The value of the entry's first field named by one of ``patterns`` (``ABSENT`` for none),
    as ``FIELD_CLASSES`` name fields: a full name or a prefix, the case set aside."""
    if e is None:
        return ABSENT
    return next((v for k, v in e.fields.items() if field_named(k, patterns)), ABSENT)


def indel_notation_of(e: Entry, hgvsc: str) -> str | None:
    """The kind of an HGVSc that names the entry's transcript (its stable id, with or without a
    version) and describes an insertion or deletion: ``del``, ``dup`` or ``ins``; ``None`` for
    any other notation."""
    name, _, body = hgvsc.partition(":")
    if name.split(".")[0] != e.feature.split(".")[0]:
        return None
    m = INDEL_NOTATION_RE.fullmatch(body)
    return m.group(3)[:3] if m else None


def shifted_indel_written_span_gate(e: Entry, actual: Entry | None, name: str, ev: str, av: str) -> bool:
    """``hgvs_shifted_indel_written_span_gate`` on one differing field: the reference's HGVSc,
    HGVS_OFFSET and HGVSIntronOffset fields are all empty, the candidate's HGVSc is an insertion
    or deletion notation of the row's transcript beside a non-zero HGVS_OFFSET, the differing
    field is one of those (an intron offset an integer, and beside a ``del`` notation the one
    its coordinates give: ``deletion_offsets_follow_notation``), and the row is an indel's
    annotated inside the transcript whose SHIFT_LENGTH, where printed, is positive. An allele of two or more
    bases at a Location of as many, under an identifier that does not spell the record's alleles,
    reads the same for an insertion and a substitution in the text formats (``indel_row`` refuses
    it); there the candidate's ``dup`` or ``ins`` notation names the row an insertion's."""
    if not blank(ev) or blank(av) or actual is None or not shift_field_consistent(e):
        return False
    if any(not blank(field_of(e, (f,))) for f in SHIFTED_INDEL_GATE_FIELDS):
        return False
    # A row annotated outside the transcript has no notation on either engine: one the
    # candidate writes there is a difference, not this class.
    if OUTSIDE_TRANSCRIPT_TERMS & set(consequence_set(e.consequence).split(",")):
        return False
    hgvsc = av if field_named(name, HGVSC_FIELDS) else field_of(actual, HGVSC_FIELDS)
    offset = av if field_named(name, HGVS_OFFSET_FIELDS) else field_of(actual, HGVS_OFFSET_FIELDS)
    if field_named(name, HGVS_INTRON_OFFSET_FIELDS) and INTEGER_RE.fullmatch(av) is None:
        return False
    kind = indel_notation_of(e, hgvsc)
    if kind is None or NON_ZERO_INTEGER_RE.fullmatch(offset) is None:
        return False
    if field_named(name, HGVS_INTRON_OFFSET_FIELDS) and kind == "del" and not deletion_offsets_follow_notation(actual, hgvsc):
        return False
    return indel_row(e) or (kind in ("dup", "ins") and equal_span_allele_row(e))


def coordinate_intron_offset(coordinate: str) -> str:
    """The intron offset of an HGVS cDNA coordinate as the HGVSIntronOffset plugin writes it:
    ``9`` for ``29+9``, ``-944`` for ``30-944``, ``0`` for an exonic coordinate."""
    offset = INTRON_OFFSET_RE.search(coordinate)
    return str(int(offset.group())) if offset else "0"


def deletion_offsets_follow_notation(a: Entry, hgvsc: str) -> bool:
    """Whether the candidate row's two HGVSIntronOffset fields are the intron offsets of the
    two coordinates of its ``del`` notation ``hgvsc``, one coordinate standing for both, in
    either order: the plugin reads them before the coordinates are ordered for printing. A
    coordinate before the stop codon's last base (``*-5``) prints an intron distance without
    its exon base, so a notation carrying one says nothing here."""
    m = DELETION_NOTATION_RE.fullmatch(hgvsc.partition(":")[2])
    if m is None:
        return False
    first, last = m.group(2), m.group(3) or m.group(2)
    if first.startswith("*-") or last.startswith("*-"):
        return True
    written = sorted(field_of(a, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS)
    return written == sorted(coordinate_intron_offset(c) for c in (first, last))


def shifted_indel_written_span_protein_gate(e: Entry, actual: Entry | None, ev: str, av: str) -> bool:
    """``hgvsp_shifted_indel_written_span_gate`` on a differing HGVSp: the reference's is empty
    and the candidate's a protein notation that neither starts at the stop codon, in
    parentheses or not, nor names it as a flank of an insertion or as the end of a duplicated,
    deleted or replaced range, neither engine's row carries a ``cds_start_NF`` or ``cds_end_NF``
    flag, the row is an indel's annotated inside the transcript with no
    cDNA, CDS or protein position and no positive SHIFT_LENGTH, both engines write the same
    HGVSc and the same non-zero HGVS_OFFSET, or the reference's are an insertion's left where
    its shift stops against a candidate that describes it further 3' or not at all
    (``insertion_left_where_the_shift_stops``), and the reference's HGVSc is an insertion or deletion
    notation of the row's transcript that places the shifted allele in the coding sequence,
    in an exon whose first base is coding base 2 or later (``shifted_exon_start``), so that
    the intron the allele is written in lies inside the coding region. The two-base allele at
    a two-base Location is read as ``shifted_indel_written_span_gate`` reads it."""
    if not blank(ev) or blank(av) or actual is None or positive_shift(e):
        return False
    if PROTEIN_NOTATION_RE.fullmatch(av) is None:
        return False
    if any(stop.match(av) is not None for stop in (STOP_CODON_NOTATION_RE, STOP_CODON_FLANK_INSERTION_RE, STOP_CODON_RANGE_END_RE)):
        return False
    if any(INCOMPLETE_CDS_FLAG_RE.search(field_of(row, FLAGS_FIELDS)) for row in (e, actual)):
        return False
    if OUTSIDE_TRANSCRIPT_TERMS & set(consequence_set(e.consequence).split(",")):
        return False
    if any(not blank(field_of(e, (f,))) for f in (*CDNA_POSITION_FIELDS, *CDS_POSITION_FIELDS, *PROTEIN_POSITION_FIELDS)):
        return False
    hgvsc, offset = field_of(e, HGVSC_FIELDS), field_of(e, HGVS_OFFSET_FIELDS)
    alike = hgvsc == field_of(actual, HGVSC_FIELDS) and offset == field_of(actual, HGVS_OFFSET_FIELDS)
    if blank(hgvsc) or not (alike or insertion_left_where_the_shift_stops(e, actual, HGVSC_FIELDS[0])):
        return False
    kind = indel_notation_of(e, hgvsc)
    if kind is None or NON_ZERO_INTEGER_RE.fullmatch(offset) is None:
        return False
    exon_start = shifted_exon_start(hgvsc.partition(":")[2], abs(int(offset)))
    if exon_start is None or exon_start < 2:
        return False
    return indel_row(e) or (kind in ("dup", "ins") and equal_span_allele_row(e))


def shifted_exon_start(body: str, shift: int) -> int | None:
    """The coding coordinate of the first base of the exon into which an HGVSc body places an
    insertion or deletion shifted ``shift`` bases out of the intron before it, or a lower bound
    of that coordinate; ``None`` where the body does not place the allele in the coding
    sequence. A ``del`` lies there when its every coordinate is a coding base (no UTR sign, no
    intron offset); a ``dup`` when it ends on a coding base, its first coordinate a coding
    base or a base of the intron before; an ``ins`` between two coding bases. A ``dup`` that
    ends on an intron base and an ``ins`` beside one are refused, and with them every
    insertion the shift leaves between the intron's last base and the exon's first base,
    the copy of that one last base included, and an insertion before the stop codon's last
    base, whose intron base carries no coordinate (``*-1``). A ``dup``'s first coordinate in
    the intron names the exon's first base, the base it counts back from or the one after
    the coding base it counts forward from, and the stop codon's last base reads as the
    largest one there; a notation of coding bases alone gives the bound
    ``start - shift + 1``, ``start`` the first base of a ``del`` or the base after a ``dup``
    or an ``ins``."""
    m = INDEL_NOTATION_RE.fullmatch(body)
    if m is None or not body.startswith("c."):
        return None
    first, last, kind = m.groups()

    def coding(coordinate: str) -> bool:
        return CODING_COORDINATE_RE.fullmatch(coordinate) is not None

    def named_exon_start(coordinate: str) -> int | None:
        m = ACCEPTOR_INTRON_COORDINATE_RE.fullmatch(coordinate)
        if m is None:
            m = DONOR_INTRON_COORDINATE_RE.fullmatch(coordinate)
            return int(m.group(1)) + 1 if m is not None else None
        return sys.maxsize if m.group(1) == "*" else int(m.group(1))

    if kind == "del":
        return int(first) - shift + 1 if coding(first) and (last is None or coding(last)) else None
    if kind == "dup":
        end = first if last is None else last
        if not coding(end):
            return None
        if last is None or coding(first):
            return int(end) + 1 - shift + 1
        return named_exon_start(first)
    if last is None or not coding(first) or not coding(last):
        return None
    return int(last) - shift + 1


def equal_span_allele_row(e: Entry) -> bool:
    """Whether the entry's allele is two or more bases over a Location of as many, with no
    record alleles to read: the text an insertion and a substitution of that length both
    print, which ``indel_row`` reads as a substitution."""
    if e.indel is not None or e.location is None or not re.fullmatch(r"[ACGTN]{2,}", e.allele):
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    return start.isdigit() and end.isdigit() and int(end) - int(start) + 1 == len(e.allele)


def cdna_coordinate_order(coordinate: str) -> tuple[int, int, int]:
    """A key ordering HGVS cDNA coordinates along the transcript: the 5' UTR (``-N``), the
    coding sequence (``N``), the 3' UTR (``*N``), each by its position, then an intron offset
    (``+M`` after the exon base, ``-M`` before it). ``*-M``, an intron distance before the stop
    codon's last base printed without that base's coordinate
    (TranscriptVariationAllele.pm:2854-2859), orders after every coding coordinate and before
    ``*1``."""
    region, number = 0, coordinate
    if coordinate.startswith("*"):
        region, number = 1, coordinate[1:]
    elif coordinate.startswith("-"):
        region, number = -1, coordinate[1:]
    m = re.fullmatch(r"(\d*)([+-]\d+)?", number)
    if m is None:
        return (region, 0, 0)
    position, offset = integer_of(m.group(1) or "0", "HGVS coordinate"), integer_of(m.group(2) or "0", "HGVS coordinate offset")
    return (region, -position if region < 0 else position, offset)


def intron_offset_abutting_insertion(e: Entry, actual: Entry | None, name: str, ev: str, av: str) -> bool:
    """``hgvs_intron_offset_abutting_insertion`` on one differing HGVSIntronOffset field: the
    reference fills both of the plugin's fields and the candidate neither, on an insertion's row
    both engines annotate outside the transcript, its one term ``upstream_gene_variant`` or
    ``downstream_gene_variant`` beside a DISTANCE of 0, with no HGVSc and no HGVS_OFFSET on
    either side. ``HGVSIntronOffset.pm:62`` asks for the notation on every Transcript row, past
    the ``within_feature`` test that keeps HGVSc off such a row (``OutputFactory.pm:1721``); an
    insertion directly beside the transcript's first or last base has one slice coordinate
    inside (``TranscriptVariationAllele.pm:2729``) and is clamped onto that base (:2732-2733),
    so the offsets stored at :1543-1545 belong to a notation the row never prints. The pair is
    0 and 0; or 1 and 0 on a downstream row of a reverse-strand transcript, where :2802 wraps
    the position past the last base to the far exon and :2822-2826 build "last base +1"; or 0
    and -1 on an upstream row of a forward-strand transcript, where ``_clip_alleles`` strips
    the first base off the inserted allele's end and moves the notation's end to the position
    before it (:2231-2238), which :2802 and :2822-2826 read as "first base -1" (a row without
    STRAND admits either). Any other pair, a row with a term inside the transcript, a value on
    the candidate's side and the fields on the candidate's side alone stay differences."""
    if actual is None or not field_named(name, HGVS_INTRON_OFFSET_FIELDS):
        return False
    if blank(ev) or not blank(av) or field_of(e, (name,)) != ev:
        return False
    if any(not blank(field_of(actual, (f,))) for f in HGVS_INTRON_OFFSET_FIELDS):
        return False
    term = consequence_set(e.consequence)
    if term not in OUTSIDE_TRANSCRIPT_TERMS or consequence_set(actual.consequence) != term:
        return False
    for entry in (e, actual):
        if field_of(entry, ("DISTANCE",)) != "0":
            return False
        if any(not blank(field_of(entry, (f,))) for f in (*HGVSC_FIELDS, *HGVS_OFFSET_FIELDS)):
            return False
    if not insertion_row(e):
        return False
    pair = tuple(field_of(e, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS)
    if pair == ("0", "0"):
        return True
    strand = field_of(e, STRAND_FIELDS)
    if term == "downstream_gene_variant":
        return strand in ("-1", ABSENT) and pair == ("1", "0")
    return strand in ("1", ABSENT) and pair == ("0", "-1")


def edited_model_shifted_insertion_bases(e: Entry, a: Entry | None, ev: str, av: str) -> bool:
    """``refseq_edited_model_shifted_insertion_bases`` on a differing HGVSc: on an insertion's
    row carrying ``BAM_EDIT`` ``OK`` whose HGVS_OFFSET is one non-zero integer on both sides,
    two ``ins`` notations of the row's transcript, one numbering and one coordinate pair, the
    reference's bases the row's Allele read on the transcript's strand (reverse-complemented where
    the row's STRAND is -1; a row without STRAND admits the allele as written or
    reverse-complemented) and the candidate's those bases turned towards the transcript's 3' end,
    which on the transcript's own strand moves the first bases to the end, by the offset modulo
    their length."""
    if blank(ev) or blank(av) or not insertion_row(e) or not any(e.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS):
        return False
    offset = field_of(e, HGVS_OFFSET_FIELDS)
    if NON_ZERO_INTEGER_RE.fullmatch(offset) is None or (a is not None and field_of(a, HGVS_OFFSET_FIELDS) != offset):
        return False
    expected_name, _, expected_body = ev.partition(":")
    actual_name, _, actual_body = av.partition(":")
    reference = EDITED_INSERTION_NOTATION_RE.fullmatch(expected_body)
    candidate = EDITED_INSERTION_NOTATION_RE.fullmatch(actual_body)
    if expected_name != actual_name or reference is None or candidate is None or indel_notation_of(e, ev) != "ins":
        return False
    if reference.group(1, 2) != candidate.group(1, 2):
        return False
    reference_bases = reference.group(3)
    allele = e.allele.upper()
    strand = field_of(e, STRAND_FIELDS)
    if strand == "1":
        on_transcript_strand = (allele,)
    elif strand == "-1":
        on_transcript_strand = (reverse_complement(allele),)
    else:
        on_transcript_strand = (allele, reverse_complement(allele))
    if reference_bases not in on_transcript_strand:
        return False
    turn = abs(int(offset)) % len(reference_bases)
    return candidate.group(3) == reference_bases[turn:] + reference_bases[:turn]


def edited_model_seam_insertion_replacement(e: Entry, a: Entry | None, ev: str, av: str) -> bool:
    """Whether ``ev`` and ``av``, the HGVSc of the reference entry ``e`` and of the candidate
    entry ``a``, are the two notations of an insertion standing on an edited model's bases: on
    an insertion's row carrying ``BAM_EDIT`` ``OK`` or ``FAILED`` on both sides, the reference's
    an ``ins`` of the row's Allele read on the transcript's strand in ``c.`` numbering, and the
    candidate's, under the reference's sequence name, a notation that describes the replacement
    of the model's bases by that allele, empty where the allele is those bases. The bases and
    their span are read from the candidate's row: the CDS position (``CDS_position``, or the
    JSON ``cds_start`` and ``cds_end``) names the span, the reference codons, which start at
    the first base of the codon the protein position names, hold the bases, and the alternate
    codons hold the allele in their place. The candidate's notation describes the replacement
    when, applied to the reference codons at plain coding positions inside them, it gives the
    alternate codons. A value that is not the entry's own HGVSc is refused."""
    if blank(ev) or a is None or not insertion_row(e):
        return False
    own = field_of(a, HGVSC_FIELDS)
    if ev != field_of(e, HGVSC_FIELDS) or (av != own and not (blank(av) and blank(own))):
        return False
    for entry in (e, a):
        if not any(entry.fields.get(f) in ("OK", "FAILED") for f in BAM_EDIT_FIELDS):
            return False
    name, _, body = ev.partition(":")
    reference = EDITED_INSERTION_NOTATION_RE.fullmatch(body)
    if reference is None or reference.group(1) != "c" or indel_notation_of(e, ev) != "ins":
        return False
    inserted = reference.group(3)
    allele = e.allele.upper()
    strand = field_of(e, STRAND_FIELDS)
    if strand == "1":
        on_transcript_strand = (allele,)
    elif strand == "-1":
        on_transcript_strand = (reverse_complement(allele),)
    else:
        on_transcript_strand = (allele, reverse_complement(allele))
    if inserted not in on_transcript_strand:
        return False
    position = a.fields.get("CDS_position")
    if position is not None:
        # Under --total_length a position carries the sequence's length after a slash.
        start, separator, end = position.partition("/")[0].partition("-")
        end = end if separator else start
        codon = a.fields.get("Protein_position", ABSENT).partition("/")[0].partition("-")[0]
    else:
        start = a.fields.get("cds_start", ABSENT)
        end = a.fields.get("cds_end", ABSENT)
        codon = a.fields.get("protein_start", ABSENT)
    if not (start.isdecimal() and end.isdecimal() and codon.isdecimal()):
        return False
    first, last, window_start = int(start), int(end), 3 * int(codon) - 2
    reference_codons, _, alternate_codons = field_of(a, CODONS_FIELDS).partition("/")
    lead, tail = first - window_start, last - window_start + 1
    if not 0 <= lead < tail <= len(reference_codons):
        return False
    bases = reference_codons[lead:tail].upper()
    if re.fullmatch(r"[ACGTN]+", bases) is None:
        return False
    if alternate_codons.upper() != (reference_codons[:lead] + inserted + reference_codons[tail:]).upper():
        return False
    if inserted == bases:
        return blank(av)
    prefix = f"{name}:c."
    if not av.startswith(prefix):
        return False
    # A ``*``, a ``-`` or an intron offset names a base outside the coding sequence, which the
    # row's codons do not hold.
    notation = re.fullmatch(
        r"(?P<x>[1-9]\d*)(?:_(?P<y>[1-9]\d*))?"
        r"(?:(?P<old>[ACGTN])>(?P<new>[ACGTN])|(?P<plain>del|dup|inv)|(?P<kind>delins|ins)(?P<added>[ACGTN]+)|\[(?P<copies>[1-9]\d*)\])",
        av[len(prefix):],
    )
    if notation is None:
        return False
    window, described = reference_codons.upper(), alternate_codons.upper()
    paired = notation.group("y") is not None
    x = int(notation.group("x")) - window_start
    y = int(notation.group("y")) - window_start if paired else x
    kind = notation.group("plain") or notation.group("kind")
    if kind == "ins":
        # An insertion reads no base, so one of its two positions may lie one base outside the
        # codons: before their first base or after their last.
        if y != x + 1 or not 0 <= y <= len(window):
            return False
        return window[:y] + notation.group("added") + window[y:] == described
    if not 0 <= x <= y < len(window):
        return False
    span = window[x:y + 1]
    if notation.group("old") is not None:
        if paired or span != notation.group("old"):
            return False
        replacement = notation.group("new")
    elif notation.group("copies") is not None:
        copies = int(notation.group("copies"))
        if not paired or copies * len(span) > len(described):
            return False
        replacement = span * copies
    elif kind == "inv":
        if not paired:
            return False
        replacement = reverse_complement(span)
    elif kind == "dup":
        replacement = span * 2
    elif kind == "delins":
        replacement = notation.group("added")
    else:
        replacement = ""
    return window[:x] + replacement + window[y + 1:] == described


# The alleles an HGVSc spells on an edited model. On a RefSeq model carrying an ``_rna_edit``
# that the alignment applied (the row carries ``BAM_EDIT`` ``OK``) ``hgvs_transcript`` clips the
# notation's alleles of the bases they share (``_clip_alleles``,
# TranscriptVariationAllele.pm:1492, :2188), which moves the coordinates and re-reads the
# operation, and then writes both alleles over whole from ``feature_seq`` (:1507-1513). An
# allele that shares a base with its reference, as a record with a second alternate leaves
# it, so carries that base in the notation: an ``ins`` of the inserted bases and the shared
# base, a ``delins`` that names a base its span does not replace, a substitution that names
# two bases a side at one position. The candidate spells what the clip leaves. The shape: an
# HGVSc on a row carrying ``BAM_EDIT`` ``OK`` and one GIVEN_REF and one USED_REF on both
# sides, the two values two substitutions, two ``ins`` or two ``delins`` of one sequence
# name, numbering and coordinate text. The reference spells the row's own pair, its USED_REF
# and its Allele read on the transcript's strand (reverse-complemented where STRAND is -1; a
# row without STRAND admits either reading): both alleles of a substitution, the Allele of
# an ``ins`` or a ``delins``, which print no reference allele. The candidate spells that
# pair without the bases the clip takes, at least one. The clip reads the genome's bases,
# which GIVEN_REF holds: it takes the longest front GIVEN_REF and the Allele share and then
# the longest back of what is left (``shared_ends``), one split and no other, so a repeat
# loses its front and the same bases turned by one are not the shape. USED_REF shares the
# bases taken, and what they leave of GIVEN_REF is as long as the positions the coordinates
# cover: one base of a substitution, whose alternate allele keeps one base, none of an
# ``ins``, and of a ``delins`` the span its coordinates state (one position for a lone
# coordinate, the count ``cdna_pair_spans`` reads from a pair; a pair that states none is
# not measured). USED_REF is longer than GIVEN_REF where the allele spans bases the model
# holds and the genome lacks, and the reference allele of a substitution there keeps those
# bases on both sides. ``renumbered`` admits two values whose coordinates differ as
# ``refseq_edited_model_hgvsc_offset`` reads them, the row's notation renumbered on the
# edited model beside its alleles. A candidate that drops a base the clip does not take, or
# more or fewer bases than the coordinates leave, another operation, a reference that does
# not spell the row's pair, a ``del`` or a ``dup`` (which print no allele), a row without
# ``BAM_EDIT`` ``OK`` (a model whose edit failed carries ``FAILED`` and is read on the
# genome), GIVEN_REF or USED_REF and the two sides reversed stay differences.
GIVEN_REF_FIELDS = ("GIVEN_REF", "given_ref")
CLIPPED_ALLELES_NOTATION_RE = re.compile(
    r"(?P<numbering>[cn])\.(?P<coordinates>\*?-?\d+(?:[+-]\d+)?(?:_\*?-?\d+(?:[+-]\d+)?)?)"
    r"(?P<alleles>(?P<old>[ACGT]+)>(?P<new>[ACGT]+)|(?P<operation>delins|ins)(?P<added>[ACGT]+))"
)


def shared_ends(old: str, new: str) -> tuple[int, int]:
    """The bases ``_clip_alleles`` (TranscriptVariationAllele.pm:2188) takes from the alleles
    ``old`` and ``new``, as a count at the front and a count at the back: the longest front
    the two share, then the longest back what is left of them shares."""
    limit = min(len(old), len(new))
    front = 0
    while front < limit and old[front] == new[front]:
        front += 1
    back = 0
    while back < limit - front and old[-1 - back] == new[-1 - back]:
        back += 1
    return front, back


def edited_model_clipped_alleles(e: Entry, a: Entry | None, ev: str, av: str, renumbered: bool = False) -> bool:
    """Whether ``ev`` and ``av``, the HGVSc of the reference entry ``e`` and of the candidate
    entry ``a``, are the two spellings of one notation on an edited model: on a row carrying
    ``BAM_EDIT`` ``OK``, one GIVEN_REF and one USED_REF on both sides, two substitutions, two
    ``ins`` or two ``delins`` of one sequence name, numbering and coordinate text, the
    reference's spelling the row's USED_REF and Allele whole on the transcript's strand and
    the candidate's spelling them without the bases the clip takes from GIVEN_REF and the
    Allele (``shared_ends``), which USED_REF shares and which leave of GIVEN_REF's length the
    positions the coordinates cover. With ``renumbered`` the coordinates may also differ as
    ``refseq_edited_model_hgvsc_offset`` accepts. A value that is not the entry's own HGVSc is
    refused."""
    if blank(ev) or blank(av) or a is None:
        return False
    if ev != field_of(e, HGVSC_FIELDS) or av != field_of(a, HGVSC_FIELDS):
        return False
    if any(not any(entry.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS) for entry in (e, a)):
        return False
    given_ref, used_ref, allele = field_of(e, GIVEN_REF_FIELDS).upper(), field_of(e, USED_REF_FIELDS).upper(), e.allele.upper()
    if (field_of(a, GIVEN_REF_FIELDS).upper(), field_of(a, USED_REF_FIELDS).upper()) != (given_ref, used_ref):
        return False
    if any(re.fullmatch(r"[ACGT]+", bases) is None for bases in (given_ref, used_ref, allele)):
        return False
    expected_name, _, expected_body = ev.partition(":")
    actual_name, _, actual_body = av.partition(":")
    reference = CLIPPED_ALLELES_NOTATION_RE.fullmatch(expected_body)
    candidate = CLIPPED_ALLELES_NOTATION_RE.fullmatch(actual_body)
    if expected_name != actual_name or reference is None or candidate is None:
        return False
    if reference.group("numbering", "operation") != candidate.group("numbering", "operation"):
        return False
    operation, coordinates = reference.group("operation", "coordinates")
    if candidate.group("coordinates") != coordinates:
        # The candidate's coordinates under the reference's alleles: the reference's notation
        # renumbered, which the offset class reads.
        moved = f"{actual_name}:{candidate.group('numbering')}.{candidate.group('coordinates')}{reference.group('alleles')}"
        if not renumbered or not field_class_accepts(REFSEQ_EDITED_MODEL_HGVSC_OFFSET, e, HGVSC_FIELDS[0], ev, moved, a):
            return False
    paired = "_" in coordinates
    # The positions the coordinates cover: ``None`` where a ``delins`` pair states no count.
    if operation is None:
        if paired:
            return False
        covered: int | None = 1
    elif operation == "ins":
        covered = 0
    elif paired:
        (span,) = cdna_pair_spans(coordinates)
        covered = None if span is None else abs(span) + 1
    else:
        covered = 1
    strand = field_of(e, STRAND_FIELDS)
    readings = []
    if strand != "-1":
        readings.append((given_ref, used_ref, allele))
    if strand != "1":
        readings.append(tuple(reverse_complement(bases) for bases in (given_ref, used_ref, allele)))
    for given, old, new in readings:
        spelled = reference.group("old", "new") == (old, new) if operation is None else reference.group("added") == new
        if not spelled:
            continue
        front, back = shared_ends(given, new)
        left = len(given) - front - back
        if not front + back or (left != covered if covered is not None else left < 1):
            continue
        if len(old) < front + back or old[:front] != new[:front] or old[len(old) - back:] != new[len(new) - back:]:
            continue
        left_old, left_new = old[front:len(old) - back], new[front:len(new) - back]
        if operation is None:
            if len(left_new) == 1 and candidate.group("old", "new") == (left_old, left_new):
                return True
        elif left_new and candidate.group("added") == left_new:
            return True
    return False


def shifted_deletion_clamped_span(e: Entry, ev: str, av: str) -> bool:
    """``hgvs_shifted_deletion_clamped_span`` on a differing HGVSc: two ``del`` notations of
    the row's transcript and one numbering, ending at the same cDNA coordinate, the reference's
    starting later than the candidate's, on a deletion's row (allele ``-``) whose HGVS_OFFSET is
    a non-zero integer on the reference side and whose SHIFT_LENGTH, where printed, is
    positive; the candidate's coordinates agree with the row's cDNA position
    (``span_matches_cdna_position``)."""
    if blank(ev) or blank(av) or e.allele != "-" or not shift_field_consistent(e):
        return False
    if NON_ZERO_INTEGER_RE.fullmatch(field_of(e, HGVS_OFFSET_FIELDS)) is None:
        return False
    expected_name, _, expected_body = ev.partition(":")
    actual_name, _, actual_body = av.partition(":")
    reference, candidate = DELETION_NOTATION_RE.fullmatch(expected_body), DELETION_NOTATION_RE.fullmatch(actual_body)
    if expected_name != actual_name or reference is None or candidate is None or indel_notation_of(e, ev) is None:
        return False
    if reference.group(1) != candidate.group(1):
        return False
    reference_start, reference_end = reference.group(2), reference.group(3) or reference.group(2)
    candidate_start, candidate_end = candidate.group(2), candidate.group(3) or candidate.group(2)
    if reference_end != candidate_end or cdna_coordinate_order(candidate_start) >= cdna_coordinate_order(reference_start):
        return False
    return span_matches_cdna_position(e, candidate.group(1), candidate_start, candidate_end)


def shifted_deletion_start_offsets(e: Entry, a: Entry, ev: str, av: str) -> bool:
    """``hgvs_shifted_deletion_clamped_span`` on a differing ``HGVS_IntronStartOffset``: the
    row's HGVSc pair has the class's shape (``shifted_deletion_clamped_span``) and each
    engine's offset is the intron offset of its own notation's start coordinate, where the two
    notations differ."""
    notations = field_of(e, HGVSC_FIELDS), field_of(a, HGVSC_FIELDS)
    if not shifted_deletion_clamped_span(e, *notations):
        return False
    starts = [DELETION_NOTATION_RE.fullmatch(notation.partition(":")[2]).group(2) for notation in notations]
    return [ev, av] == [coordinate_intron_offset(start) for start in starts]


class ReferenceSequence:
    """The reference FASTA ``hgvs_shift_sequence_edge_window`` reads: an uncompressed FASTA
    beside its ``.fai`` index, a region's length and its bases by 1-based inclusive
    coordinates, upper-cased and clipped to the region."""

    def __init__(self, path: Path) -> None:
        if path.suffix == ".gz":
            raise ComparatorError(f"--reference-fasta: {path} is compressed; an uncompressed FASTA beside its .fai is read")
        index = Path(f"{path}.fai")
        try:
            lines = index.read_text(encoding="utf-8").splitlines()
            self.handle = path.open("rb")
        except OSError as exc:
            raise ComparatorError(f"--reference-fasta: cannot read {exc.filename}: {exc.strerror or exc}") from exc
        self.regions: dict[str, tuple[int, int, int, int]] = {}
        for line in lines:
            cols = line.split("\t")
            if len(cols) < 5 or not all(c.isdigit() for c in cols[1:5]) or int(cols[3]) == 0:
                raise ComparatorError(f"{index}: not a FASTA index line: {line}")
            self.regions[cols[0]] = (int(cols[1]), int(cols[2]), int(cols[3]), int(cols[4]))

    def length(self, region: str) -> int | None:
        entry = self.regions.get(region)
        return entry[0] if entry else None

    def bases(self, region: str, start: int, end: int) -> str:
        length, offset, per_line, width = self.regions[region]
        start, end = max(start, 1), min(end, length)
        if start > end:
            return ""
        span = end - start + 1
        self.handle.seek(offset + (start - 1) // per_line * width + (start - 1) % per_line)
        raw = self.handle.read(span + (span // per_line + 2) * (width - per_line))
        return raw.translate(None, b"\r\n")[:span].decode("ascii").upper()


def perform_shift_walk(head: str, tail: str, size: int, pre: str, post: str, forward: bool, whole_flanks: bool = False) -> int:
    """``perform_shift`` (TranscriptVariationAllele.pm:303-349) over the flanks given: how
    many bases an allele of ``size`` bases moves 3', forward along ``post`` from its first
    base or, for a reverse-strand transcript, back along ``pre`` from its last. The allele
    turns by one base per step, so the walk along ``post`` reads its first bases (``head``)
    and the walk along ``pre`` its last (``tail``); the loop bound is the flank's length less
    the allele's, one more along ``pre``, and ``post``'s length where that is negative
    (:303-304). ``whole_flanks`` bounds the walk as two flanks of 1,000 bases do and ends it
    where the flank given ends: a region's edge cuts the allele's own flank short without
    lowering the bound on what the sequence admits."""
    pre_length, post_length = (SEQUENCE_EDGE_WINDOW, SEQUENCE_EDGE_WINDOW) if whole_flanks else (len(pre), len(post))
    if forward:
        limit = post_length - size
        if limit < 0:
            limit = post_length
        n = 0
        while n <= limit and n < len(post) and post[n] == head[n % len(head)]:
            n += 1
        return n
    limit = pre_length - size + 1
    if limit < 0:
        limit = post_length
    n = 1
    while n <= limit and n <= len(pre) and pre[len(pre) - n] == tail[(len(tail) - n) % len(tail)]:
        n += 1
    return n - 1


def indel_span(e: Entry, record: Entry | None) -> tuple[str, int, int] | None:
    """The region, start and end of an insertion's or deletion's row as Ensembl VEP numbers
    them, an insertion's start one past its end: the Location of a text row; the record line's
    CHROM, POS, REF and ALT of a VCF entry, one alternate sharing the reference's first base;
    the record object's ``seq_region_name``, ``start`` and ``end`` of a JSON entry. ``None``
    for any other row."""
    deletion = e.allele == "-"
    if e.location is not None:
        region, _, span = e.location.rpartition(":")
        low, _, high = span.partition("-")
        if not region or not low.isdigit() or not (high or low).isdigit():
            return None
        first, last = int(low), int(high or low)
        if deletion:
            return region, first, last
        return (region, last, first) if last == first + 1 else None
    fields = record.fields if record is not None else {}
    if "POS" in fields:
        reference, alternate = fields.get("REF", ""), fields.get("ALT", "")
        if not fields["POS"].isdigit() or "," in alternate or not reference or reference[:1] != alternate[:1]:
            return None
        position = int(fields["POS"])
        if deletion:
            return (fields["CHROM"], position + 1, position + len(reference) - 1) if len(alternate) == 1 < len(reference) else None
        return (fields["CHROM"], position + 1, position) if len(reference) == 1 and alternate[1:] == e.allele else None
    start, end = fields.get("start", ""), fields.get("end", "")
    if "seq_region_name" not in fields or not start.isdigit() or not end.isdigit():
        return None
    if (int(start) <= int(end)) if deletion else (int(start) == int(end) + 1):
        return fields["seq_region_name"], int(start), int(end)
    return None


class EdgeTranscript(NamedTuple):
    """A transcript of the candidate JSON cache as ``hgvs_shift_sequence_edge_window`` reads
    it: the name a notation gives it, its span and strand, its exons in genomic order and, on
    a coding transcript, the cDNA positions of the coding sequence's first and last base. A
    coordinate is counted along the span from the transcript's 5' end, introns included, as
    the transcript slice counts it (TranscriptVariationAllele.pm:2704-2713)."""

    name: str
    start: int
    end: int
    strand: int
    exons: tuple[tuple[int, int], ...]
    coding: tuple[int, int] | None

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def coordinate(self, position: int) -> int:
        return position - self.start + 1 if self.strand == 1 else self.end - position + 1

    def notation_coordinate(self, coordinate: int) -> str | None:
        """The HGVS coordinate of a coordinate inside an exon: its cDNA position, on a coding
        transcript counted from the coding sequence's first base, negative before it and
        starred after its last base; ``None`` in an intron or outside the transcript."""
        position = self.start + coordinate - 1 if self.strand == 1 else self.end - coordinate + 1
        cdna = 0
        for start, end in self.exons if self.strand == 1 else reversed(self.exons):
            if start <= position <= end:
                cdna += position - start + 1 if self.strand == 1 else end - position + 1
                break
            cdna += end - start + 1
        else:
            return None
        if self.coding is None:
            return str(cdna)
        first, last = self.coding
        if cdna < first:
            return str(cdna - first)
        return f"*{cdna - last}" if cdna > last else str(cdna - first + 1)

    def notation(self, first: int, last: int, change: str) -> str | None:
        """The HGVSc body, numbering included, of ``change`` over the coordinates ``first`` to
        ``last``; ``None`` where either lies outside an exon."""
        ends = [self.notation_coordinate(coordinate) for coordinate in ((first,) if first == last else (first, last))]
        if None in ends:
            return None
        return f"{'n' if self.coding is None else 'c'}.{'_'.join(ends)}{change}"


def edge_transcript(record: dict) -> EdgeTranscript | None:
    """One transcript record of a JSON cache shard as an ``EdgeTranscript``; the coding
    bounds are the record's own or its mapper's. A notation names the transcript by its stable
    id, the version after a dot where the record has one and the id ends in none
    (TranscriptVariationAllele.pm:1498-1500)."""
    mapper = (record.get("variation_effect_feature_cache") or {}).get("mapper") or {}
    coding = [record.get(key, mapper.get(key)) for key in ("cdna_coding_start", "cdna_coding_end")]
    name, version = str(record.get("stable_id")), record.get("version")
    if version is not None and not re.search(r"\.\d+$", name) and "LRG" not in name:
        name = f"{name}.{version}"
    try:
        exons = tuple(sorted((int(exon["start"]), int(exon["end"])) for exon in record.get("exons") or []))
        bounds = None if None in coding else (int(coding[0]), int(coding[1]))
        return EdgeTranscript(name, int(record["start"]), int(record["end"]), int(record["strand"]), exons, bounds)
    except (KeyError, TypeError, ValueError):
        return None


class EdgeWalks(NamedTuple):
    """One row's two walks along its strand: over the allele's own flanks and over the
    clipped slice's two cuts."""

    own: int
    cut: int


class SequenceEdges:
    """What ``hgvs_shift_sequence_edge_window`` reads beside the entries: the reference
    sequence, for each row's two walks, and the candidate JSON cache's transcripts, for
    where each walk leaves the notation on the row's transcript (``SEQUENCE_EDGE_WINDOW``'s
    comment)."""

    def __init__(self, fasta: Path, cache: Path) -> None:
        self.reference = ReferenceSequence(fasta)
        if not (cache / "transcripts").is_dir():
            raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
        self.transcripts: dict[tuple[str, str], EdgeTranscript] = {}
        for region, record in structural_variant_comparator().iter_cache_transcripts(cache):
            transcript = edge_transcript(record)
            if transcript is not None:
                self.transcripts[(region, record["stable_id"])] = transcript
        self._walks: dict[tuple[str, int, int, str, bool], EdgeWalks | None] = {}

    def walks(self, e: Entry, record: Entry | None) -> EdgeWalks | None:
        """The row's two walks; ``None`` for a row that is no insertion's or deletion's on a
        known region, or whose strand walks towards an edge 1,000 bases away or more, where
        the cut flank is the allele's own."""
        strand = field_of(e, STRAND_FIELDS)
        span = indel_span(e, record)
        if span is None or strand not in ("1", "-1"):
            return None
        region, start, end = span
        length = self.reference.length(region)
        if length is None or not 0 < min(start, end) <= max(start, end) <= length:
            return None
        forward = strand == "1"
        key = (region, start, end, e.allele, forward)
        if key not in self._walks:
            self._walks[key] = self._walk(region, length, start, end, e.allele, forward)
        return self._walks[key]

    def _walk(self, region: str, length: int, start: int, end: int, allele: str, forward: bool) -> EdgeWalks | None:
        window = SEQUENCE_EDGE_WINDOW
        if (end + window <= length) if forward else (start - window >= 1):
            return None
        bases = self.reference.bases
        if allele == "-":
            head, tail, size = bases(region, start, min(end, start + window)), bases(region, max(start, end - window), end), end - start + 1
        elif re.fullmatch(r"[ACGTacgt]+", allele):
            head = tail = allele.upper()
            size = len(allele)
        else:
            return None
        low, high = max(start - window, 1), min(end + window, length)
        own = perform_shift_walk(
            head, tail, size, bases(region, start - window, start - 1), bases(region, end + 1, end + window), forward, whole_flanks=True
        )
        cut = perform_shift_walk(
            head, tail, size, bases(region, low, min(low + window - 1, high)), bases(region, max(high - window + 1, low), high), forward
        )
        return EdgeWalks(own, cut)

    def _deletion(self, transcript: EdgeTranscript, walks: EdgeWalks, start: int, end: int) -> tuple[str | None, str | None]:
        """The two notation bodies of a deletion's row, the empty string where none is
        written. Ensembl VEP clamps the written span onto the transcript and adds its walk,
        writing nothing for a span wholly outside or a clamped end the walk carries past the
        last base (TranscriptVariationAllele.pm:2729-2733, :1472, :1477-1478); the candidate
        moves the span by its walk and clamps what is left inside, writing nothing where
        nothing is."""
        length = transcript.length
        first, last = sorted((transcript.coordinate(start), transcript.coordinate(end)))
        if last < 1 or first > length or length < min(last, length) + walks.cut:
            reference: str | None = ""
        else:
            reference = transcript.notation(max(first, 1) + walks.cut, min(last, length) + walks.cut, "del")
        if last + walks.own < 1 or first + walks.own > length:
            return reference, ""
        return reference, transcript.notation(max(first + walks.own, 1), min(last + walks.own, length), "del")

    def _insertion(
        self, transcript: EdgeTranscript, walks: EdgeWalks, region: str, start: int, end: int, allele: str
    ) -> tuple[str | None, str | None]:
        """The two notation bodies of the row of an insertion between two bases of the
        transcript, the empty string where none is written. At its walk's place each side
        duplicates the bases before the insertion point where they repeat the allele as the
        walk turned it and inserts that allele where they do not; Ensembl VEP duplicates the
        bases after the point instead where those repeat it (Utils/Sequence.pm:573-596) and
        writes nothing where its walk carries the point past the last base
        (TranscriptVariationAllele.pm:1472), as the candidate writes nothing there. At the
        point after the last base a plain insertion has no second coordinate (:1524-1527)
        and neither side writes one."""
        forward, length, size = transcript.strand == 1, transcript.length, len(allele)
        point = transcript.coordinate(end if forward else start)
        if not 1 <= point < length:
            return None, None

        def repeats(steps: int, after: bool) -> bool:
            if forward == after:
                first = end + steps + 1 if forward else start - steps
            else:
                first = end + steps - size + 1 if forward else start - steps - size
            return self.reference.bases(region, first, first + size - 1) == turned_allele(allele, steps, forward)

        def placed(steps: int, after_first: bool) -> str | None:
            at = point + steps
            if at > length:
                return ""
            if after_first and at + size <= length and repeats(steps, True):
                return transcript.notation(at + 1, at + size, "dup")
            if at - size + 1 >= 1 and repeats(steps, False):
                return transcript.notation(at - size + 1, at, "dup")
            if at == length:
                return ""
            turned = turned_allele(allele, steps, forward)
            return transcript.notation(at, at + 1, f"ins{turned if forward else reverse_complement(turned)}")

        return placed(walks.cut, True), placed(walks.own, False)

    def accepts(self, e: Entry, a: Entry | None, record: Entry | None, name: str) -> bool:
        """``hgvs_shift_sequence_edge_window`` on the differing field ``name`` of an insertion's
        or deletion's Transcript row, the candidate's entry ``a`` beside the reference's ``e``:
        each side's HGVSc is the notation its walk leaves on the row's transcript, or empty
        where it leaves none, and its HGVS_OFFSET that walk's length beside a notation or an
        HGVSp (OutputFactory.pm:1755, :1767), so either field is accepted only where all four
        values are those and each notation names the transcript as the cache's record does. A
        text row whose allele is as long as its Location prints alike for an insertion and a
        substitution; its two notations name it."""
        if a is None or not field_named(name, SEQUENCE_EDGE_FIELDS) or e.indel is False:
            return False
        if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript":
            return False
        span, walks = indel_span(e, record), self.walks(e, record)
        if span is None or walks is None or walks.own == walks.cut:
            return False
        region, start, end = span
        # Under --transcript_version the Feature is the stable id with its version after a dot
        # (OutputFactory.pm:1465); the cache files the transcript by the stable id alone.
        transcript = self.transcripts.get((region, e.feature)) or self.transcripts.get((region, e.feature.rpartition(".")[0]))
        forward = field_of(e, STRAND_FIELDS) == "1"
        if transcript is None or (transcript.strand == 1) != forward:
            return False
        if e.allele == "-":
            bodies = self._deletion(transcript, walks, start, end)
        else:
            bodies = self._insertion(transcript, walks, region, start, end, e.allele.upper())
        for entry, body, walked in ((e, bodies[0], walks.cut), (a, bodies[1], walks.own)):
            notation = field_of(entry, HGVSC_FIELDS)
            if body is None or ("" if blank(notation) else notation) != (f"{transcript.name}:{body}" if body else ""):
                return False
            shown = walked if body or not blank(field_of(entry, HGVSP_FIELDS)) else 0
            offset = field_of(entry, HGVS_OFFSET_FIELDS)
            if (0 if blank(offset) else int(offset) if INTEGER_RE.fullmatch(offset) else None) != (shown if forward else -shown):
                return False
        return True


def turned_allele(allele: str, steps: int, forward: bool) -> str:
    """An inserted allele after ``steps`` steps of ``perform_shift``: its first base moved to
    its end at each step of a forward walk (TranscriptVariationAllele.pm:338-341), its last to
    its front at each step back (:331-334)."""
    turn = steps % len(allele) if forward else -steps % len(allele)
    return allele[turn:] + allele[:turn]


def cdna_open_side(e: Entry) -> str | None:
    """The side of the row's cDNA position that lies outside the transcript: ``start`` for a
    ``?`` before the dash of the text formats' column (``?-6``), ``end`` for one after it
    (``8485-?``), in the JSON format the side of ``cdna_start`` and ``cdna_end`` left out while
    the other is present (JSON.pm:340-341); ``None`` for a closed position."""
    text_position_field, json_start_field = CDNA_START_FIELDS
    if text_position_field in e.fields:
        start, _, end = e.fields[text_position_field].partition("/")[0].partition("-")
        if start == "?":
            return "start"
        return "end" if end == "?" else None
    has_start = not blank(e.fields.get(json_start_field, ABSENT))
    has_end = not blank(e.fields.get(CDNA_END_FIELD, ABSENT))
    if has_start == has_end:
        return None
    return "end" if has_start else "start"


def clipped_reference_against_whole_alternate(e: Entry, ev: str, av: str, edges: TranscriptEdges | None = None) -> bool:
    """``hgvs_clipped_reference_against_whole_alternate`` on a differing HGVSc: on a row of a
    sequence allele whose cDNA position has an open side (``cdna_open_side``), the reference's
    notation names the row's transcript, or is empty where the whole alternate equals the bases
    inside (Utils/Sequence.pm:529-531), and the candidate's names the same transcript and
    numbering, or is empty where the aligned change lies outside. Both engines clip their
    alleles against the bases the clamp left inside the transcript, so both notations lie within
    the allele's span over the transcript, Ensembl VEP's at most one base before it where the
    outside bases print as an insertion before exon 1 and an ``ins`` of either at most one base
    after it: every coordinate of the candidate's lies within the span's length plus one of every
    coordinate of the reference's (``edge_row_reach``; a row whose reach cannot be read, the end
    open and no Location, passes on the transcript and numbering alone). Against an empty
    reference the candidate is a ``del``, the bases the 3'-aligned change deletes inside; against
    an empty candidate the reference is an ``ins`` or a ``dup``, what the outside bases read as
    once the clip has stripped the unchanged bases inside, or, on a row whose allele is shorter
    than its span, a ``del``, a ``delins`` or a substitution (``deletion_over_the_edge``): the
    candidate describes no deletion that runs over the transcript's edge. The kinds and bases of the two
    notations are otherwise unconstrained: Ensembl VEP's alignment of the whole alternate against
    the bases inside gives any kind at any span inside, and the golden harness pins vep-rs's exact
    notation on every corpus row of the class. A substitution's row whose two notations show
    the clamp against the transcript's exon spans (``whole_alternate_over_fewer_bases``, which
    reads ``edges``, the candidate cache) is the class without the coordinate test; a row whose
    cDNA position is empty, the allele's other end lying in an intron or past the other edge,
    is the class only so, and without a cache not at all."""
    side = cdna_open_side(e)
    if e.allele == "-":
        return False
    reference = clipped_notation(ev, e.feature)
    candidate = clipped_notation(av, e.feature)
    clamp_shown = whole_alternate_over_fewer_bases(e, reference, candidate, edges)
    if side is None and not clamp_shown:
        return False
    if blank(ev):
        return candidate is not None and candidate.group("plain") == "del"
    if reference is None:
        return False
    if blank(av):
        return reference.group("kind") == "ins" or reference.group("plain") == "dup" or deletion_over_the_edge(e, side, reference)
    if candidate is None or ev.partition(":")[0] != av.partition(":")[0]:
        return False
    if candidate.group("numbering") != reference.group("numbering"):
        return False
    if clamp_shown:
        return True
    reach = edge_row_reach(e, side)
    if reach is None:
        return True
    return all(
        abs(c - r) <= reach
        for c in notation_coordinate_indexes(candidate)
        for r in notation_coordinate_indexes(reference)
    )


def insertion_without_reference_notation(e: Entry, candidate: re.Match) -> bool:
    """Whether the row is one on which the candidate prints an ``ins`` where the reference
    prints no HGVSc: a sequence allele's row on a reverse-strand transcript (STRAND ``-1``)
    whose cDNA position is open at its start (``cdna_open_side``), the transcript's first
    base lying at its genomic end, the candidate's notation an ``ins`` of the bases a longer
    allele adds inside the span (``insertion_inside_the_span``). There the clip leaves
    Ensembl VEP an insertion at slice position 0, past the last exon in genomic order, for
    which ``_get_cDNA_position`` returns no coordinate (TranscriptVariationAllele.pm:2834),
    so its row carries no notation and the HGVSIntronOffset plugin writes no field beside
    it. The test places the candidate's two fields; the row's HGVSc is no pair of the class.
    A forward-strand transcript's row, a position open at its end and a row without STRAND
    are refused."""
    if e.allele == "-" or cdna_open_side(e) != "start" or field_of(e, STRAND_FIELDS) != "-1":
        return False
    return insertion_inside_the_span(e, candidate)


def insertion_inside_the_span(e: Entry, candidate: re.Match) -> bool:
    """Whether the candidate's ``ins`` notation (a ``CLIPPED_NOTATION_RE`` match) writes the
    bases a longer allele adds inside the span of a row whose cDNA position is open at its
    start (``?-N``): the row's allele is a sequence longer than its Location, the notation
    inserts as many bases as it is longer, between two exon coordinates next to each other,
    and under ``n.`` numbering, whose coordinates are cDNA positions, at or beside a base
    inside: its first coordinate at most the position's closed end. Under ``c.`` numbering
    the coordinates' place is not read, the coding start not being among the row's fields; a
    row without a Location, whose span is not stated, is refused."""
    first, last = candidate.group("start"), candidate.group("end")
    if candidate.group("kind") != "ins" or last is None or intron_coordinate(first) or intron_coordinate(last):
        return False
    if cdna_coordinate_index(last) - cdna_coordinate_index(first) != 1:
        return False
    if e.location is None or re.fullmatch(r"[ACGTN]+", e.allele) is None:
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    if not start.isdigit() or (end and not end.isdigit()):
        return False
    span_bases = integer_of(end, "Location position") - integer_of(start, "Location position") + 1 if end else 1
    if len(candidate.group("bases")) != len(e.allele) - span_bases:
        return False
    if candidate.group("numbering") != "n":
        return True
    text_position_field, _json_start_field = CDNA_START_FIELDS
    _, _, high = e.fields.get(text_position_field, "").partition("/")[0].partition("-")
    closed = high or e.fields.get(CDNA_END_FIELD, "")
    if not closed.isdigit() or not first.isdigit() or not last.isdigit():
        return False
    return integer_of(first, "HGVS coordinate") <= integer_of(closed, "cDNA position")


def deletion_over_the_edge(e: Entry, side: str | None, reference: re.Match) -> bool:
    """Whether the reference's notation is one the class takes against a candidate that prints
    none on the row of an allele that deletes bases over the transcript's edge: the row's cDNA
    position has an open side, its allele holds fewer bases than its Location spans, and the
    notation is a ``del``, a ``delins`` or a substitution, what ``hgvs_variant_notation`` makes
    of the whole alternate against the bases the clamp left inside
    (TranscriptVariationAllele.pm:1474-1483, 2732-2733). The candidate cuts the alternate as
    the clamp cut the reference, describes a deletion lying inside the transcript there and one
    lying outside not at all, so its empty notation beside such a reference notation is that of
    a deletion with no window wholly inside the transcript. A row without a Location, an allele as
    long as its span or longer and a closed position stay differences."""
    if side is None or e.location is None:
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    if not (start.isdigit() and end.isdigit()):
        return False
    span_bases = integer_of(end, "Location position") - integer_of(start, "Location position") + 1
    if not 0 < len(e.allele) < span_bases:
        return False
    return reference.group("plain") == "del" or reference.group("kind") == "delins" or reference.group("alternate_base") is not None


def exonic_run(first: str, last: str) -> int | None:
    """The count of bases from the cDNA coordinate ``first`` to ``last``, both exonic and in
    one region of the transcript (5' UTR, coding sequence or 3' UTR), ``None`` where either
    carries an intron offset or the two lie in different regions: the scale
    ``cdna_coordinate_index`` reads does not carry a distance across the stop codon or an
    intron."""
    (first_region, _, first_offset), (last_region, _, last_offset) = cdna_coordinate_order(first), cdna_coordinate_order(last)
    if first_region != last_region or first_offset or last_offset:
        return None
    return cdna_coordinate_index(last) - cdna_coordinate_index(first) + 1


def added_bases_left_at_the_span_end(e: Entry, actual: Entry | None, name: str) -> bool:
    """``hgvs_added_bases_left_at_the_span_end`` on one differing field of a row whose allele
    holds more bases than its Location spans, the bases an allele adds to a reference it
    keeps whole, where the reference prints an HGVSc of the row's transcript and no
    HGVS_OFFSET: its record does not shift (TranscriptVariationAllele.pm:120-123). Three
    shapes, read on the two rows together.

    Past a reverse-strand transcript's last base: the reference's notation is an ``ins`` of
    the added bases between a coordinate and the same coordinate with an intron offset of 1
    (:2802, :2820-2827), beside which it writes the plugin's start offset 1 and end offset 0.
    The candidate's is a duplication of as many bases ending at that coordinate, with no
    HGVS_OFFSET and its own two intron offsets, where the last added base repeats the kept
    base beside it: the HGVSc and the start offset are the class. Or the candidate writes no
    HGVSc, no HGVS_OFFSET and neither field, where the last added base does not repeat that
    kept base, so the added bases have no placement inside the transcript: the HGVSc and
    both fields are the class.

    Carried along a run: the reference's notation is an ``ins`` between two neighbouring
    exonic bases, or a ``dup``, and the candidate's a ``dup`` of as many bases or an ``ins``
    of the same bases turned, a number of bases further 3' that its HGVS_OFFSET states,
    negative on a reverse strand, the two rows writing the same intron offsets: the HGVSc
    and the HGVS_OFFSET are the class. Where nothing moves, the candidate's is the ``dup``
    of the bases that end at the reference's ``ins``, beside no HGVS_OFFSET: the HGVSc is
    the class.

    A row whose allele holds no more bases than its span, a reference HGVS_OFFSET, a notation of
    another transcript or numbering, a count of bases that is not the allele's surplus, and
    every other pair of values stay differences."""
    if actual is None or e.location is None or re.fullmatch(r"[ACGT]+", e.allele) is None:
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    end = end or start
    if not (start.isdigit() and end.isdigit()):
        return False
    added = len(e.allele) - (integer_of(end, "Location position") - integer_of(start, "Location position") + 1)
    reference = clipped_notation(field_of(e, HGVSC_FIELDS), e.feature)
    if added <= 0 or reference is None or not blank(field_of(e, HGVS_OFFSET_FIELDS)):
        return False
    candidate_hgvsc, candidate_offset = field_of(actual, HGVSC_FIELDS), field_of(actual, HGVS_OFFSET_FIELDS)
    candidate = clipped_notation(candidate_hgvsc, e.feature)
    if not blank(candidate_hgvsc) and (candidate is None or candidate.group("numbering") != reference.group("numbering")):
        return False
    reference_offsets = tuple(field_of(e, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS)
    candidate_offsets = tuple(field_of(actual, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS)
    strand = field_of(e, STRAND_FIELDS)
    first, last = reference.group("start"), reference.group("end")
    past = last is not None and next((m.group(1) for c, o in ((first, last), (last, first))
                                      if (m := ONE_PAST_COORDINATE_RE.fullmatch(c)) and m.group(1) == o), None)
    if reference.group("kind") == "ins" and past:
        bases = reference.group("bases")
        if len(bases) != added or reverse_complement(e.allele[:added]) != bases:
            return False
        # A run without the plugin writes the two fields on neither row.
        without_plugin = all(blank(v) for v in (*reference_offsets, *candidate_offsets))
        if not (reference_offsets == ("1", "0") or without_plugin) or strand not in ("-1", ABSENT):
            return False
        repeats_the_last_base = e.allele[0] == e.allele[added]
        if blank(candidate_hgvsc):
            if repeats_the_last_base or not blank(candidate_offset) or any(not blank(v) for v in candidate_offsets):
                return False
            return field_named(name, HGVSC_FIELDS) or field_named(name, HGVS_INTRON_OFFSET_FIELDS)
        duplicated = (candidate.group("start"), candidate.group("end") or candidate.group("start"))
        if candidate.group("plain") != "dup" or duplicated[1] != past or exonic_run(*duplicated) != added:
            return False
        if not repeats_the_last_base or not blank(candidate_offset) or not (candidate_offsets == ("0", "0") or without_plugin):
            return False
        return field_named(name, HGVSC_FIELDS) or field_named(name, HGVS_INTRON_START_OFFSET_FIELDS)
    if candidate is None or not (blank(candidate_offset) or NON_ZERO_INTEGER_RE.fullmatch(candidate_offset)):
        return False
    if reference_offsets != candidate_offsets or blank(reference_offsets[0]) != blank(reference_offsets[1]):
        return False
    # The point the reference describes the added bases at, as the exonic base before it,
    # and the point the candidate does. An ``ins`` past the stop codon prints its 3' flank
    # first (TranscriptVariationAllele.pm:1535-1538).
    if reference.group("kind") == "ins":
        if last is None or len(reference.group("bases")) != added:
            return False
        if exonic_run(first, last) == 2:
            described = first
        elif exonic_run(last, first) == 2:
            described = last
        else:
            return False
    elif reference.group("plain") == "dup":
        if exonic_run(first, last or first) != added:
            return False
        described = last or first
    else:
        return False
    placed_first, placed_last = candidate.group("start"), candidate.group("end")
    if candidate.group("plain") == "dup":
        if exonic_run(placed_first, placed_last or placed_first) != added:
            return False
        placed = placed_last or placed_first
    elif candidate.group("kind") == "ins":
        if placed_last is None or exonic_run(placed_first, placed_last) != 2 or len(candidate.group("bases")) != added:
            return False
        placed = placed_first
    else:
        return False
    carried = exonic_run(described, placed)
    if blank(candidate_offset):
        typed = reference.group("kind") == "ins" and candidate.group("plain") == "dup" and carried == 1
        return typed and field_named(name, HGVSC_FIELDS)
    if carried is None or carried < 2:
        return False
    moved = integer_of(candidate_offset.lstrip("-"), "HGVS_OFFSET")
    if moved != carried - 1 or candidate_offset.startswith("-") != (strand == "-1") and strand != ABSENT:
        return False
    if candidate.group("kind") == "ins":
        inserted = reference.group("bases") if reference.group("kind") == "ins" else None
        if inserted is None or candidate.group("bases") != turned_allele(inserted, moved, True):
            return False
    return field_named(name, HGVSC_FIELDS) or field_named(name, HGVS_OFFSET_FIELDS)


def added_bases_written_outside_the_transcript(e: Entry, actual: Entry | None, name: str, ev: str, av: str) -> bool:
    """``hgvs_added_bases_written_outside_the_transcript`` on one differing field of a row
    both engines annotate outside the transcript, its one term ``upstream_gene_variant`` or
    ``downstream_gene_variant``, whose allele is an insertion's or holds more bases than its
    Location spans: the reference prints no HGVSc and no HGVS_OFFSET (OutputFactory.pm:1721)
    and writes the HGVSIntronOffset plugin's two fields or neither (HGVSIntronOffset.pm:62);
    the candidate prints a ``dup`` or an ``ins`` of the row's transcript, of as many bases as
    the allele holds or adds to its span where its coordinates count them, and both plugin
    fields, the intron offsets of that notation's coordinates in their printed order, or
    neither where the reference writes neither, as in a run without the plugin. On a
    downstream row nothing moved: the candidate's is a ``dup`` beside no
    HGVS_OFFSET, the duplication of the transcript's last bases. On an upstream row the 3'
    rule carried the bases inside: its HGVS_OFFSET is a non-zero integer, negative on a
    reverse strand, larger in size than the row's DISTANCE on an insertion's row and no
    smaller than it on the row of an allele that holds more bases than its span, whose
    DISTANCE counts from that span. The differing HGVSc, HGVS_OFFSET and
    plugin fields are the class. A row with a term inside the transcript, a reference HGVSc
    or HGVS_OFFSET, one plugin field on the reference's side without the other, a candidate
    without a notation, one without both fields beside a reference that writes them, fields
    that are not its notation's, an
    HGVS_OFFSET on a downstream row or none on an upstream one, and another kind of notation
    stay differences. A VCF or JSON entry is read at the Location its record gives its allele
    (``at_its_location``) and stays a difference where its record gives it none."""
    e = at_its_location(e)
    if actual is None or e.location is None or re.fullmatch(r"[ACGT]+", e.allele) is None:
        return False
    span_start, _, span_end = e.location.split(":", 1)[-1].partition("-")
    span_end = span_end or span_start
    if not (span_start.isdigit() and span_end.isdigit()):
        return False
    positions = integer_of(span_end, "Location position") - integer_of(span_start, "Location position") + 1
    surplus = len(e.allele) - positions
    counts = ({len(e.allele)} if insertion_row(e) else set()) | ({surplus} if surplus > 0 else set())
    if not counts:
        return False
    term = consequence_set(e.consequence)
    if term not in OUTSIDE_TRANSCRIPT_TERMS or consequence_set(actual.consequence) != term:
        return False
    if any(not blank(field_of(e, (f,))) for f in (*HGVSC_FIELDS, *HGVS_OFFSET_FIELDS)):
        return False
    reference_offsets = [field_of(e, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS]
    if blank(reference_offsets[0]) != blank(reference_offsets[1]):
        return False
    candidate = clipped_notation(field_of(actual, HGVSC_FIELDS), e.feature)
    if candidate is None or not (candidate.group("plain") == "dup" or candidate.group("kind") == "ins"):
        return False
    first, last = candidate.group("start"), candidate.group("end")
    if candidate.group("kind") == "ins":
        if last is None:
            return False
        counted = len(candidate.group("bases"))
    else:
        counted = exonic_run(first, last or first)
    if counted is not None and counted not in counts:
        return False
    if any(c.startswith("*-") for c in (first, last or first)):
        return False
    written = tuple(field_of(actual, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS)
    # A run without the plugin writes the two fields on neither row.
    without_plugin = all(blank(v) for v in (*reference_offsets, *written))
    if not without_plugin and written != (coordinate_intron_offset(first), coordinate_intron_offset(last or first)):
        return False
    offset, strand = field_of(actual, HGVS_OFFSET_FIELDS), field_of(e, STRAND_FIELDS)
    if term == "downstream_gene_variant":
        if not blank(offset) or candidate.group("plain") != "dup":
            return False
    else:
        distance = field_of(e, DISTANCE_FIELDS)
        if NON_ZERO_INTEGER_RE.fullmatch(offset) is None or not distance.isdigit():
            return False
        # An insertion's DISTANCE counts from the nearer base beside its point, a longer allele's from its span.
        added_to_a_span = surplus > 0 and (positions != 2 or counted == surplus)
        if integer_of(offset.lstrip("-"), "HGVS_OFFSET") < integer_of(distance, "DISTANCE") + (0 if added_to_a_span else 1):
            return False
        if strand != ABSENT and offset.startswith("-") != (strand == "-1"):
            return False
    if field_named(name, HGVS_INTRON_OFFSET_FIELDS):
        return not blank(av) and field_of(actual, (name,)) == av and ev != av
    return blank(ev) and not blank(av) and (field_named(name, HGVSC_FIELDS) or field_named(name, HGVS_OFFSET_FIELDS))


def printed_intron_offsets(notation: re.Match) -> tuple[str, str] | None:
    """The intron offsets of the two coordinates of a notation (a ``CLIPPED_NOTATION_RE``
    match) in their printed order, one coordinate standing for both: the pair the
    HGVSIntronOffset plugin writes beside an insertion's ``ins`` or ``dup``, whose span
    ``hgvs_variant_notation`` holds from the 5' flank to the 3' flank (Utils/Sequence.pm
    :577-579). ``None`` for a coordinate before the stop codon's last base (``*-5``), whose
    exon part the pattern cannot read."""
    first, last = notation.group("start"), notation.group("end") or notation.group("start")
    if any(c.startswith("*-") for c in (first, last)):
        return None
    return coordinate_intron_offset(first), coordinate_intron_offset(last)


def shift_bound(inserted: int, towards_5prime: bool = False) -> int:
    """The steps ``perform_shift`` walks an insertion of ``inserted`` bases at most: its loop
    runs over the 1,000 bases of flank less the allele's length, inclusive, and over the whole
    flank for a longer allele (TranscriptVariationAllele.pm:303-306). Going 5', on a
    reverse-strand transcript, the loop counts from 1 to the flank less the allele's length
    plus 1, which is 0 for an allele of 1,001 bases and is left at that (:304 resets a
    negative limit alone): that allele is not walked at all."""
    if inserted <= SHIFT_FLANK_BASES:
        return SHIFT_FLANK_BASES + 1 - inserted
    return 0 if towards_5prime and inserted == SHIFT_FLANK_BASES + 1 else SHIFT_FLANK_BASES


def insertion_left_where_the_shift_stops(e: Entry, actual: Entry | None, name: str) -> bool:
    """``hgvs_insertion_left_where_the_shift_stops`` on one differing field of an insertion's
    row whose reference HGVS_OFFSET is, in size, exactly the bound of ``perform_shift`` for
    the allele's length (``shift_bound``) beside an ``ins`` or a ``dup`` of the row's
    transcript: the walk ended on its bound, not on a base that differs. Where that bound is
    0, an allele of 1,001 bases on a row whose STRAND is -1, the reference walked nothing and
    its row carries no HGVS_OFFSET; that row is read only where its Location names the two
    bases around an insertion's point (a VCF or JSON entry's the Location its record gives
    its allele, ``Entry.allele_location``), the record's alleles, where the entry carries
    them, hold no reference base, and the reference's notation is the ``ins`` of the allele's
    own bases or a ``dup`` of as many. The reference maps each end of the copy to a coordinate
    alone (TranscriptVariationAllele.pm:1524-1525), so two exonic coordinates count the copied
    bases less those of every intron between them: no more than the allele holds. A ``dup``
    that counts two fewer, the bases an allele that keeps the two bases of its span adds, and
    a ``dup`` of one coordinate stay differences. Two shapes of candidate row.

    Further 3': a ``dup`` or an ``ins`` in the same numbering beside an HGVS_OFFSET of the
    same sign, negative where the reference carries none, and a larger size. A ``dup`` whose coordinates are exonic names as many bases
    as the allele holds; an ``ins`` writes the allele's bases as the transcript's strand
    reads them, turned by the steps its HGVS_OFFSET states. The HGVSc, the HGVS_OFFSET and
    the plugin's two fields are the class.

    Past the last base: no HGVSc and neither plugin field, beside no HGVS_OFFSET or the
    reference's own, which a row with an HGVSp keeps. The HGVSc and the two fields are the
    class, and the HGVS_OFFSET where the candidate writes none.

    Where the rows carry the plugin's fields, each row's two are the intron parts of its own
    notation's coordinates in their printed order (``printed_intron_offsets``). A deletion's row, a reference
    HGVS_OFFSET of any other size, a candidate HGVS_OFFSET no larger than the bound and
    every other pair of values stay differences."""
    if actual is None or re.fullmatch(r"[ACGT]+", e.allele) is None:
        return False
    stopped = field_of(e, HGVS_OFFSET_FIELDS)
    reference = clipped_notation(field_of(e, HGVSC_FIELDS), e.feature)
    if reference is None or not (reference.group("kind") == "ins" or reference.group("plain") == "dup"):
        return False
    if blank(stopped):
        # The one walk that never starts: its row has no HGVS_OFFSET to read a direction from.
        # An allele that keeps bases of its span carries none either and is never walked
        # (:123), so the row must be an insertion's own: a Location of the two bases around
        # its point and, where the record spells its alleles, no reference base among them.
        span_start, _, span_end = (at_its_location(e).location or "").split(":", 1)[-1].partition("-")
        if not (span_start.isdigit() and span_end.isdigit() and int(span_end) == int(span_start) + 1) or e.insertion is False:
            return False
        towards_5prime = field_of(e, STRAND_FIELDS) == "-1"
        bound = shift_bound(len(e.allele), towards_5prime)
        if bound != 0:
            return False
        # Described where the record writes it, the insertion is its own bases as the
        # transcript's strand reads them, or a copy of as many. Each end of the copy is mapped
        # to a coordinate alone (:1524-1525), so two exonic coordinates count the copied bases
        # less those of every intron between them. An allele that keeps the two bases of its
        # span adds two fewer than it holds, and a ``dup`` that counts those stays a difference.
        if reference.group("kind") == "ins":
            if reference.group("bases") != reverse_complement(e.allele):
                return False
        else:
            copied = exonic_run(reference.group("start"), reference.group("end") or reference.group("start"))
            if copied is not None and (not 2 <= copied <= len(e.allele) or copied == len(e.allele) - 2):
                return False
    else:
        if NON_ZERO_INTEGER_RE.fullmatch(stopped) is None:
            return False
        towards_5prime = stopped.startswith("-")
        bound = shift_bound(len(e.allele), towards_5prime)
        if integer_of(stopped.lstrip("-"), "HGVS_OFFSET") != bound:
            return False
    reference_offsets = tuple(field_of(e, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS)
    candidate_offsets = tuple(field_of(actual, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS)
    # A run without the plugin writes the two fields on neither row.
    without_plugin = all(blank(v) for v in (*reference_offsets, *candidate_offsets))
    if not without_plugin and reference_offsets != printed_intron_offsets(reference):
        return False
    candidate_hgvsc, moved = field_of(actual, HGVSC_FIELDS), field_of(actual, HGVS_OFFSET_FIELDS)
    if blank(candidate_hgvsc):
        if any(not blank(v) for v in candidate_offsets) or not (blank(moved) or moved == stopped):
            return False
        if field_named(name, HGVS_OFFSET_FIELDS):
            return blank(moved)
        return field_named(name, HGVSC_FIELDS) or field_named(name, HGVS_INTRON_OFFSET_FIELDS)
    candidate = clipped_notation(candidate_hgvsc, e.feature)
    if candidate is None or candidate.group("numbering") != reference.group("numbering"):
        return False
    if NON_ZERO_INTEGER_RE.fullmatch(moved) is None or moved.startswith("-") != towards_5prime:
        return False
    further = integer_of(moved.lstrip("-"), "HGVS_OFFSET") - bound
    if further <= 0 or not (without_plugin or candidate_offsets == printed_intron_offsets(candidate)):
        return False
    if candidate.group("plain") == "dup":
        duplicated = exonic_run(candidate.group("start"), candidate.group("end") or candidate.group("start"))
        if duplicated is not None and duplicated != len(e.allele):
            return False
    elif candidate.group("kind") == "ins":
        placed = reverse_complement(e.allele) if moved.startswith("-") else e.allele
        if candidate.group("bases") != turned_allele(placed, bound + further, True):
            return False
    else:
        return False
    return any(field_named(name, fields) for fields in (HGVSC_FIELDS, HGVS_OFFSET_FIELDS, HGVS_INTRON_OFFSET_FIELDS))


def notation_intron_offsets(notation: re.Match) -> tuple[str, str] | None:
    """The two values the HGVSIntronOffset plugin writes beside a notation (a
    ``CLIPPED_NOTATION_RE`` match), the start's and the end's: the intron parts of the two
    coordinates as ``hgvs_transcript`` reads them (TranscriptVariationAllele.pm:1530-1533, an
    absent part read as 0; stored at :1543 and :1545) before it orders the pair for printing
    (:1535-1538). The routine holds an ``ins`` from its 3' flank to its 5' flank, the clip
    leaving the start one above the end (:2245-2248), and turns the two unless the end holds
    a star, so the start is the coordinate printed second, or the one printed first where the
    5' flank lies past the stop codon; it holds a ``dup`` and every other kind in its printed
    order, one coordinate standing for both. ``None`` for a pair the routine would have printed the other way
    round, an ``ins`` of one coordinate, and a coordinate before the stop codon's last base
    (``*-5``), whose exon part the pattern cannot read."""
    first, last = notation.group("start"), notation.group("end")
    insertion = notation.group("kind") == "ins"
    if (insertion and last is None) or any(c.startswith("*-") for c in (first, last or first)):
        return None
    first_parts, last_parts = (
        tuple(int(part or 0) for part in COORDINATE_PARTS_AS_READ_RE.search(c).groups()) for c in (first, last or first)
    )
    if insertion and last_parts > first_parts and "*" not in first:
        return str(last_parts[1]), str(first_parts[1])
    if first_parts > last_parts and "*" not in (last or first):
        return None
    return str(first_parts[1]), str(last_parts[1])


def insertion_before_the_first_base_against_duplication(
    e: Entry, actual: Entry, reference: re.Match | None, candidate: re.Match, pairs: list[tuple[str, ...]]
) -> bool:
    """Whether the two rows are the reference's ``ins`` from one base before a forward-strand
    transcript's first base against the candidate's duplication of the bases the allele adds,
    the fourth shape of ``clipped_notation_intron_offsets``, read from the two rows alone in
    four parts, all of which hold.

    The row: its allele holds more bases than its Location spans, ``A`` bases over ``S``
    positions; its STRAND is ``1``; its cDNA position is open at its start and closed at its
    end, ``?-N``, with ``N`` above 0 and below ``S``: the span's last base is the
    transcript's base ``N`` and ``S`` - ``N`` bases of the span lie before the transcript.

    The reference: an ``ins`` of the row's transcript whose two coordinates are one
    coordinate ``X`` carrying the intron part ``-1`` and the same ``X`` carrying none; it
    writes the allele's bases less the last ``N``; it carries no HGVS_OFFSET; the plugin's
    pair is 0 and -1. An ``ins`` the clip leaves has consumed every reference base
    ``hgvs_transcript`` handed it (TranscriptVariationAllele.pm:2231-2242, :2271), so it
    was handed ``N`` bases, fewer than the span holds, which the clamp alone does
    (:2732-2733). A clamp at the last base leaves the position open at its end, so the
    clamp is at the first base: the bases handed are the transcript's first ``N``, the
    insertion ends on the base before them, and ``X`` is the first base's coordinate.

    The candidate: a ``dup`` of the same transcript and numbering whose coordinates carry
    no intron part and cover as many bases as the allele adds, ``A`` - ``S``
    (``exonic_run``); its HGVS_OFFSET a whole number ``k``, an empty one read as 0; the
    plugin's pair 0 and 0.

    The tie: from ``X`` to the duplication's last coordinate are ``N`` + ``k`` exonic bases
    of one region (``exonic_run``). The duplication ends ``k`` bases 3' of the span's last
    base, the transcript's base ``N``, which lies that far from ``X`` only where ``X`` is
    base 1.

    The ``ins`` the clip leaves at a true intron's last base fails the count: nothing
    clamped, it writes ``A`` - ``S`` bases, and a clamped span that reaches across a whole
    intron hands the clip the intron's bases too, so it writes fewer than ``A`` - ``N``. A
    plain insertion at an exon's first base carries the pair -1 and 0 beside a closed
    position (Utils/Sequence.pm:592-593). A duplication across the start codon from ``X``
    stays a difference, ``exonic_run`` reading one region. A row that carries BAM_EDIT is
    refused: on an edited model ``hgvs_transcript`` respells the notation's alleles after
    the clip (:1510-1514), so the bases written are not the clip's. A row without a
    Location or without STRAND and a reverse-strand transcript's row, on which the
    reference writes no notation at such a base (:2783, :2834), are refused."""
    if reference is None or e.location is None or re.fullmatch(r"[ACGT]+", e.allele) is None:
        return False
    if field_of(e, STRAND_FIELDS) != "1" or any(not blank(field_of(entry, BAM_EDIT_FIELDS)) for entry in (e, actual)):
        return False
    start, _, end = e.location.split(":", 1)[-1].partition("-")
    end = end or start
    if not (start.isdigit() and end.isdigit()):
        return False
    span_bases = integer_of(end, "Location position") - integer_of(start, "Location position") + 1
    added = len(e.allele) - span_bases
    text_position_field, _json_start_field = CDNA_START_FIELDS
    closed = e.fields.get(text_position_field, "").partition("/")[0].partition("-")[2] or e.fields.get(CDNA_END_FIELD, "")
    if added < 1 or cdna_open_side(e) != "start" or not closed.isdigit():
        return False
    inside = integer_of(closed, "cDNA position")
    first, last = reference.group("start"), reference.group("end")
    if not 0 < inside < span_bases or reference.group("kind") != "ins" or last is None:
        return False
    first_base = next((m.group(1) for c, o in ((first, last), (last, first))
                       if (m := ONE_BEFORE_COORDINATE_RE.fullmatch(c)) and m.group(1) == o), None)
    if first_base is None or reference.group("bases") != e.allele[: len(e.allele) - inside]:
        return False
    if not blank(field_of(e, HGVS_OFFSET_FIELDS)) or pairs[0] != ("0", "-1"):
        return False
    if candidate.group("plain") != "dup" or candidate.group("numbering") != reference.group("numbering"):
        return False
    duplicated_last = candidate.group("end") or candidate.group("start")
    if exonic_run(candidate.group("start"), duplicated_last) != added or pairs[1] != ("0", "0"):
        return False
    carried = field_of(actual, HGVS_OFFSET_FIELDS)
    if not blank(carried) and re.fullmatch(r"[0-9]+", carried) is None:
        return False
    return exonic_run(first_base, duplicated_last) == inside + (0 if blank(carried) else integer_of(carried, "HGVS_OFFSET"))


def clipped_notation_intron_offsets(
    e: Entry, actual: Entry | None, name: str, ev: str, av: str, edges: TranscriptEdges | None = None
) -> bool:
    """``hgvs_clipped_reference_against_whole_alternate`` on one differing HGVSIntronOffset
    field, in four shapes. ``HGVSIntronOffset.pm:62-72`` write the two fields beside whatever
    notation ``hgvs_transcript`` returns, each the intron part of a coordinate as the plugin
    reads it (``notation_intron_offsets``), and a row without a notation gets neither. The
    reference alone fills both fields, on a row whose HGVSc the class takes
    (``clipped_reference_against_whole_alternate``) against a candidate that prints none, the
    reference's an ``ins`` or a ``dup``, or the ``del``,
    ``delins`` or substitution of an allele that deletes bases over the edge. The candidate alone fills both, on a row whose
    reference prints no HGVSc where the candidate prints an ``ins`` inside the span at a
    reverse-strand transcript's first base (``insertion_without_reference_notation``); the
    HGVSc of such a row is not the class's and stays a difference. Or both fill both, on a
    row whose two notations the class takes as a pair, each an ``ins``, and the field named
    differs. In each of these three a side's two values are those intron parts of its own
    notation. Or both fill both and the end's field alone is asked, on a forward-strand
    transcript's row whose reference prints an ``ins`` from one base before the transcript's
    first base, beside the pair 0 and -1, and whose candidate prints the duplication of the
    bases the allele adds, beside the pair 0 and 0 and the bases the HGVS 3' rule carried
    it in HGVS_OFFSET (``insertion_before_the_first_base_against_duplication``): the
    reference's -1 is the intron part of a coordinate written for a base that lies in no
    intron. That shape takes neither the start's field nor the row's HGVS_OFFSET, and the
    row's HGVSc is the class's only by its own test.
    Any other pair, a candidate notation of another kind beside an empty reference or a
    filled one, a reference ``dup`` beside the candidate's fields, and, outside the fourth
    shape, a row either side's HGVS_OFFSET is written on beside the candidate's fields and a
    row of two notations whose HGVSc pair the class does not take stay differences."""
    if actual is None:
        return False
    pairs = [tuple(field_of(entry, (f,)) for f in HGVS_INTRON_OFFSET_FIELDS) for entry in (e, actual)]
    filled = [all(not blank(value) for value in pair) for pair in pairs]
    empty = [all(blank(value) for value in pair) for pair in pairs]
    hgvsc, candidate_hgvsc = field_of(e, HGVSC_FIELDS), field_of(actual, HGVSC_FIELDS)
    reference, candidate = clipped_notation(hgvsc, e.feature), clipped_notation(candidate_hgvsc, e.feature)
    if filled[0] and empty[1]:
        if blank(ev) or not blank(av) or field_of(e, (name,)) != ev or blank(hgvsc) or not blank(candidate_hgvsc):
            return False
        if not clipped_reference_against_whole_alternate(e, hgvsc, candidate_hgvsc, edges):
            return False
        offsets = notation_intron_offsets(reference)
        return offsets is not None and offsets == pairs[0]
    if not filled[1] or blank(av) or field_of(actual, (name,)) != av:
        return False
    if candidate is not None and candidate.group("plain") == "dup":
        if not field_named(name, HGVS_INTRON_END_OFFSET_FIELDS) or field_of(e, (name,)) != ev:
            return False
        return insertion_before_the_first_base_against_duplication(e, actual, reference, candidate, pairs)
    if candidate is None or candidate.group("kind") != "ins":
        return False
    if any(not blank(field_of(entry, HGVS_OFFSET_FIELDS)) for entry in (e, actual)):
        return False
    if notation_intron_offsets(candidate) != pairs[1]:
        return False
    if empty[0]:
        return blank(ev) and blank(hgvsc) and insertion_without_reference_notation(e, candidate)
    if not clipped_reference_against_whole_alternate(e, hgvsc, candidate_hgvsc, edges):
        return False
    if not filled[0] or blank(ev) or field_of(e, (name,)) != ev or ev == av or reference is None:
        return False
    if reference.group("kind") != "ins":
        return False
    return notation_intron_offsets(reference) == pairs[0]


def whole_alternate_over_fewer_bases(
    e: Entry, reference: re.Match | None, candidate: re.Match | None, edges: TranscriptEdges | None
) -> bool:
    """Whether the two notations of a substitution's row show the clamp: the reference's
    ``delins`` holds the bases of the row's allele that lie outside the transcript and that the
    clamp cut from the candidate's. ``edges`` places the allele on its transcript from the exon
    spans of the candidate cache (``TranscriptEdges``): the bases outside at each end, the
    transcript's strand and the exon or intron each base lies in; without a cache, or where it
    cannot place the row, nothing is shown and the row stays a difference. Each engine clips
    its alleles at the ends alone (``_clip_alleles``, TranscriptVariationAllele.pm:2204-2242),
    so each ``delins`` writes one run of the allele read on the transcript's strand. The
    candidate's run starts as many bases into the allele as the clamp cut before the
    transcript's first base, then the unchanged bases its clip strips; the reference's, the
    whole alternate read against the bases inside (:1474-1483), starts at the bases its own
    clip strips. The offset of the candidate's run past the reference's, plus the signed
    distance from the candidate's start coordinate to the reference's, is therefore the count
    of bases cut before the first base, and the same holds past the last base at the other
    end. The row is the class when the two counts are the cache's, one of them positive, the
    candidate's run lies among the allele's bases inside the transcript, and each coordinate
    of either notation is written as the exon spans place its base (``coordinates_as_placed``):
    an exon base without an intron offset, an intron base with its own, and the exon bases
    the coordinates are numbered from as far apart along the exons as their numbers say. The
    candidate's two coordinates name the first and last base of its run; the reference's name
    the base its run starts at once the bases cut before the first base are passed, and the
    base its run ends at less the bases cut past the last. At an end the clamp did
    not cut, both engines clip the same pairs of bases, so the two notations name the same
    coordinate there. Coordinates are compared start with start and end with end on
    ``cdna_coordinate_index``'s scale, so an intron or the stop codon inside a notation's span
    does not enter a count; where the stop codon or the midpoint of an intron lies between the
    two start coordinates, or between the two end coordinates, the scale misstates their
    distance, the count is not the cache's and the row stays a difference. A coordinate
    written across the midpoint from the base its run places it at can leave the count the
    cache's on that scale, and ``coordinates_as_placed`` refuses it. So can a coordinate
    written across the stop codon from the other notation's: one not written with ``*``
    whose place on that scale, or its exon base's, is a billion or more lies where a ``*``
    coordinate does, and ``coordinates_as_placed`` refuses it by that place. Refused: an
    insertion's or deletion's row (``indel_row``) and an
    allele that is not a sequence in capitals; a cDNA position that is closed, its two ends
    mapping to a cDNA base, or not written (``cdna_position_empty``); a row whose STRAND is
    not the cache's; a row of an edited model (``BAM_EDIT`` ``OK``: its alleles are respelled
    whole after the clip, :1510-1514); a notation of another kind on either side; and a
    reference without an exon position (every coordinate of an intronic substitution names an
    intron base, ``intron_coordinate``). A candidate that clips a changed base at a clamped
    end and moves its coordinate with it reads as one that clipped an unchanged base, and
    where a run lies more than once in the allele any pair of offsets that fits is taken:
    telling these apart takes the reference sequence, which this rule does not read; the
    golden harness pins the candidate's bases on the corpus rows. A ``c.`` numbering shifted
    whole against the transcript's own is not told from it either, the cache's coding start
    not being read."""
    if edges is None or (cdna_open_side(e) is None and not cdna_position_empty(e)):
        return False
    if indel_row(e) or any(e.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS):
        return False
    if reference is None or candidate is None or reference.group("kind") != "delins" or candidate.group("kind") != "delins":
        return False
    if all(intron_coordinate(c) for c in (reference.group("start"), reference.group("end")) if c):
        return False
    allele = e.allele
    if not re.fullmatch(r"[ACGTN]+", allele):
        return False
    placed = edges.place(e)
    if placed is None:
        return False
    before_first, past_last = placed.cuts()
    strand = field_of(e, STRAND_FIELDS)
    if not (before_first or past_last) or (strand in ("1", "-1") and (strand == "-1") != placed.reverse):
        return False
    read = reverse_complement(allele) if placed.reverse else allele
    whole, cut = reference.group("bases"), candidate.group("bases")
    reference_span, candidate_span = notation_coordinate_indexes(reference), notation_coordinate_indexes(candidate)
    reference_end, candidate_end = (m.group("end") or m.group("start") for m in (reference, candidate))
    if (not before_first and reference.group("start") != candidate.group("start")) or (not past_last and reference_end != candidate_end):
        return False
    # The two counts fix how far the candidate's run starts past the reference's and how much
    # shorter it is, so one pass over the reference's offsets finds every pair that fits.
    run_shift = before_first - (reference_span[0] - candidate_span[0])
    if len(whole) - len(cut) != run_shift + past_last - (candidate_span[-1] - reference_span[-1]):
        return False
    cut_offsets = set(run_offsets(read, cut))
    for whole_at in run_offsets(read, whole):
        cut_at = whole_at + run_shift
        if cut_at not in cut_offsets or cut_at < before_first or cut_at + len(cut) > len(read) - past_last:
            continue
        written = (
            (candidate.group("start"), cut_at),
            (candidate_end, cut_at + len(cut) - 1),
            (reference.group("start"), before_first + whole_at),
            (reference_end, whole_at + len(whole) - past_last - 1),
        )
        if len(whole) > before_first + past_last and coordinates_as_placed(placed, written, reference.group("numbering")):
            return True
    return False


@dataclass(frozen=True)
class PlacedAllele:
    """An entry's allele beside its transcript: the allele's first and last genomic position,
    the transcript's exon spans in genomic order and whether it lies on the reverse strand."""

    start: int
    end: int
    exons: tuple[tuple[int, int], ...]
    reverse: bool

    def cuts(self) -> tuple[int, int]:
        """The bases of the allele before the transcript's first base and past its last base,
        on the transcript's strand."""
        first, last = min(start for start, _ in self.exons), max(end for _, end in self.exons)
        low, high = max(first - self.start, 0), max(self.end - last, 0)
        return (high, low) if self.reverse else (low, high)

    def exon_base(self, offset: int) -> tuple[int, int]:
        """The exon base that the allele's base ``offset`` bases along the transcript's strand
        is numbered from, as its place along the exons counted from the transcript's first
        base, and the intron offset from it: 0 in an exon; in an intron the distance past the
        last base of the exon before it or, where the exon after it is nearer, that distance
        below zero from that exon's first base (TranscriptVariationAllele.pm:2798-2829). A
        base outside the exons' span has place 0."""
        position = self.end - offset if self.reverse else self.start + offset
        place, before = 0, None
        for start, end in reversed(self.exons) if self.reverse else self.exons:
            if start <= position <= end:
                return place + (end - position if self.reverse else position - start) + 1, 0
            if before is not None and (end < position < before[0] if self.reverse else before[1] < position < start):
                past = before[0] - position if self.reverse else position - before[1]
                ahead = position - end if self.reverse else start - position
                return (place, past) if past <= ahead else (place + 1, -ahead)
            place, before = place + end - start + 1, (start, end)
        return 0, 0

    def intron_offset(self, offset: int) -> int:
        """The intron offset of the allele's base ``offset`` bases along the transcript's
        strand (``exon_base``)."""
        return self.exon_base(offset)[1]


class TranscriptEdges:
    """Where an entry's allele lies on its transcript, from the exon spans of the candidate
    JSON cache (``--json-cache``): spans and lengths alone, no sequence. A text entry's span is
    its Location; a VCF or JSON entry's is its input record's POS and REF length (``inputs``),
    less the shared first base where the parser drops it. A shard of the cache is read when an
    entry first needs it."""

    def __init__(self, cache: Path, inputs: list[str] | None = None) -> None:
        if not is_directory(cache / "transcripts"):
            raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
        self.cache = cache
        self.inputs = inputs
        self._chromosomes: dict[str, list[tuple[int, Path]]] = {}
        self._shards: dict[Path, dict[str, list[tuple[tuple[tuple[int, int], ...], bool]]]] = {}

    def allele_span(self, e: Entry) -> tuple[str, int, int] | None:
        """The chromosome and the first and last genomic position of the entry's allele. A
        record of several alternates, one of another length than REF, whose alleles all start
        with one base is read without that base, one position on (Parser/VCF.pm:301-321); an
        allele the record does not write so, as one --minimal trims, has no span here."""
        if e.location is not None:
            chrom, _, span = e.location.rpartition(":")
            start, _, end = span.partition("-")
            if not start.isdigit() or (end and not end.isdigit()):
                return None
            return chrom, integer_of(start, "Location position"), integer_of(end or start, "Location position")
        if self.inputs is None or e.record is None or e.record >= len(self.inputs):
            return None
        cells = self.inputs[e.record].split("\t")
        if len(cells) < 5 or not cells[1].isdigit():
            return None
        reference, alternates = cells[3].upper(), cells[4].upper().split(",")
        chopped = (
            len(alternates) > 1
            and any(len(a) != len(reference) for a in alternates)
            and len({a[:1] for a in (reference, *alternates) if "*" not in a}) == 1
        )
        written = [a[1:] if chopped and "*" not in a else a for a in alternates]
        if e.allele not in written or len(e.allele) != len(reference) - chopped:
            return None
        start = integer_of(cells[1], "input record POS") + chopped
        return cells[0], start, start + len(e.allele) - 1

    def chromosome_shards(self, chrom: str) -> list[tuple[int, Path]]:
        """The shards of a chromosome, each with the first position of the range its name
        gives, the highest first. The directory is the chromosome's name as written, without
        or with ``chr``, or the mitochondrion's other name."""
        if chrom not in self._chromosomes:
            bare = chrom[3:] if chrom[:3].lower() == "chr" else chrom
            other = {"M": "MT", "MT": "M"}.get(bare, bare)
            names = [n for n in dict.fromkeys((chrom, bare, f"chr{bare}", other, f"chr{other}")) if n and "/" not in n and n != ".."]
            directory = next((d for d in (self.cache / "transcripts" / n for n in names) if is_directory(d)), None)
            shards = []
            for shard in () if directory is None else [*directory.glob("*.json"), *directory.glob("*.json.gz")]:
                low, _, high = shard.name.partition(".")[0].partition("-")
                if low.isdigit() and high.isdigit():
                    shards.append((integer_of(low, f"{shard}: the range start in the shard's name"), shard))
            self._chromosomes[chrom] = sorted(shards, reverse=True)
        return self._chromosomes[chrom]

    def transcript(self, chrom: str, start: int, end: int, feature: str) -> tuple[tuple[tuple[int, int], ...], bool] | None:
        """The exon spans and strand of the transcript ``feature`` where its exons meet
        ``start`` to ``end``. A converted cache files a transcript in every shard its gene
        overlaps and a built one in the shard of its first base alone, so the shards are read
        from the allele's own toward the chromosome's start until one files the transcript. The
        id looked up is the Feature as written, then without its version."""
        for low, shard in self.chromosome_shards(chrom):
            if low > end:
                continue
            if shard not in self._shards:
                self._shards[shard] = self.shard_exons(shard)
            for stable_id in dict.fromkeys((feature, feature.partition(".")[0])):
                for exons, reverse in self._shards[shard].get(stable_id, ()):
                    if min(first for first, _ in exons) <= end and start <= max(last for _, last in exons):
                        return exons, reverse
        return None

    @staticmethod
    def shard_exons(shard: Path) -> dict[str, list[tuple[tuple[tuple[int, int], ...], bool]]]:
        """Each transcript of a shard by stable id, once per locus the shard files the id at:
        its exon spans in genomic order and whether its strand is -1. A transcript without an
        id or without exons is left out. A shard that is not JSON or is nested deeper than the
        interpreter reads, and a transcript whose id is no string, whose exons are no array of
        objects or whose exon bound is no integer, end the run with the comparator's error."""
        try:
            parsed = json.loads(read_text(shard))
        except RecursionError as exc:
            raise ComparatorError(f"{shard} is nested too deep to read") from exc
        except ValueError as exc:
            raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
        transcripts: dict[str, list[tuple[tuple[tuple[int, int], ...], bool]]] = defaultdict(list)
        for t in NearestSites.shard_transcripts(shard, parsed):
            stable_id, exons = t.get("stable_id"), t.get("exons")
            if not stable_id or not exons:
                continue
            if not isinstance(stable_id, str) or not isinstance(exons, list) or not all(x is None or isinstance(x, dict) for x in exons):
                raise ComparatorError(f"{shard}: a transcript has an id that is no string or exons that are no array of objects")
            bounds = [x.get(k) for x in exons if x is not None for k in ("start", "end")]
            if any(isinstance(b, bool) or not isinstance(b, (int, str)) or not str(b).isdecimal() for b in bounds):
                raise ComparatorError(f"{shard}: transcript {stable_id[:40]!r} has an exon without an integer start and end")
            bound_of = f"{shard}: an exon bound of transcript {stable_id[:40]!r},"
            spans = sorted(tuple(sorted(integer_of(str(b), bound_of) for b in pair)) for pair in zip(bounds[::2], bounds[1::2]))
            if spans:
                transcripts[stable_id].append((tuple(spans), t.get("strand") in (-1, "-1")))
        return dict(transcripts)

    def place(self, e: Entry) -> PlacedAllele | None:
        """The entry's allele beside its transcript; ``None`` where the allele's span or the
        transcript's exons cannot be read, the allele is not as long as its span, or the two do
        not meet."""
        span = self.allele_span(e)
        if span is None:
            return None
        chrom, start, end = span
        if end - start + 1 != len(e.allele):
            return None
        transcript = self.transcript(chrom, start, end, e.feature)
        if transcript is None:
            return None
        exons, reverse = transcript
        return PlacedAllele(start, end, exons, reverse)


def intron_coordinate(coordinate: str) -> bool:
    """Whether an HGVS cDNA coordinate names an intron base: it ends with an intron offset, or
    is ``*-M``, the distance before the stop codon's last base that Ensembl VEP prints without
    that base's coordinate (TranscriptVariationAllele.pm:2854-2859)."""
    return INTRON_OFFSET_RE.search(coordinate) is not None or coordinate.startswith("*-")


def written_as_placed(coordinate: str, intron_offset: int) -> bool:
    """Whether an HGVS cDNA coordinate is written as the exon spans place its base
    (``PlacedAllele.intron_offset``): without an intron offset in an exon, with that offset in
    an intron. Beside the stop codon's last base Ensembl VEP writes an intron base without that
    base's coordinate, ``*N`` after it and ``*-N`` before it
    (TranscriptVariationAllele.pm:2854-2859)."""
    if not intron_offset:
        return not intron_coordinate(coordinate)
    written = INTRON_OFFSET_RE.search(coordinate)
    if written is not None and integer_of(written.group(), "HGVS coordinate offset") == intron_offset:
        return True
    return coordinate == f"*{intron_offset}"


def coordinates_as_placed(placed: PlacedAllele, written: Collection[tuple[str, int]], numbering: str) -> bool:
    """Whether each coordinate of two notations is written as the exon spans place its base;
    ``written`` pairs a coordinate with the offset of its base along the allele on the
    transcript's strand. Each coordinate has the intron offset of its base
    (``written_as_placed``), and the exon bases the coordinates are numbered from
    (``PlacedAllele.exon_base``) lie as far apart along the exons as the numbers say: on
    ``cdna_coordinate_index``'s scale every coordinate before the stop codon gives one and the
    same first place of the numbering, and every ``*`` coordinate one and the same last base of
    the coding sequence. Under ``n.`` a number is the base's own place. A coordinate not
    written with ``*`` whose place on that scale, or its exon base's, is a billion or more
    lies where a ``*`` coordinate does, across the stop codon from it, and is refused. Where
    the coding sequence starts is not read from the cache, so a ``c.`` numbering shifted
    whole against the transcript's own is not told from it."""
    origins: dict[bool, set[int]] = {False: set(), True: set()}
    for coordinate, offset in written:
        place, intron_offset = placed.exon_base(offset)
        if not written_as_placed(coordinate, intron_offset):
            return False
        index = cdna_coordinate_index(coordinate) - intron_offset
        after_stop = coordinate.startswith("*")
        if not after_stop and max(index, index + intron_offset) >= 10**9:
            return False
        origins[after_stop].add(place - (index - 10**9 if after_stop else index))
    if numbering == "n":
        return origins == {False: {1}, True: set()}
    return all(len(found) <= 1 for found in origins.values())


def cdna_position_empty(e: Entry) -> bool:
    """Whether the row writes a cDNA position and it is empty: the text formats' column present
    and blank before the transcript's length that --total_length appends (``-/1500``), or a
    JSON object, which leaves an empty position out (JSON.pm:327), with neither ``cdna_start``
    nor ``cdna_end``. A text row whose columns leave the position out says nothing of it."""
    text_position_field, json_start_field = CDNA_START_FIELDS
    if text_position_field in e.fields:
        return blank(e.fields[text_position_field].partition("/")[0])
    _, json_hgvsc_field = HGVSC_FIELDS
    return json_hgvsc_field in e.fields and all(blank(e.fields.get(f, ABSENT)) for f in (json_start_field, CDNA_END_FIELD))


def run_offsets(allele: str, bases: str) -> list[int]:
    """Every offset at which ``bases`` lies in ``allele``, overlapping ones included."""
    return [m.start() for m in re.finditer(f"(?={re.escape(bases)})", allele)]


def stop_codon_protein_comparison(consequence: str, flags: str, ev: str, av: str) -> bool:
    """``hgvsp_stop_codon_protein_comparison`` on a differing HGVSp: both values name one
    protein, the row's FLAGS (``flags``) name no ``cds_end_NF``, and the two descriptions, an
    escaped ``=`` and a pair of prediction parentheses set aside, are one of the class's nine
    shapes; ``consequence`` is the reference row's terms, which the last five shapes read:
    ``stop_lost`` in the fifth to the eighth, ``stop_retained_variant`` in the ninth."""
    if blank(ev) or blank(av) or "cds_end_NF" in flags:
        return False
    reference_name, separator, reference = ev.partition(":p.")
    candidate_name, _, candidate = av.partition(":p.")
    if not separator or reference_name != candidate_name:
        return False
    reference, candidate = reference.replace("%3D", "="), candidate.replace("%3D", "=")
    if reference.startswith("(") and reference.endswith(")") and candidate.startswith("(") and candidate.endswith(")"):
        reference, candidate = reference[1:-1], candidate[1:-1]
    return (
        stop_replaced_by_a_stop(reference, candidate)
        or stop_replaced_by_first_residues(reference, candidate)
        or stop_replaced_through_an_inserted_stop(reference, candidate)
        or residues_inserted_before_the_stop(reference, candidate)
        or lost_stop_written_as_unchanged(consequence, reference, candidate)
        or range_to_the_stop_replaced(consequence, reference, candidate)
        or unknown_count_where_a_stop_is_met(consequence, reference, candidate)
        or stop_deleted_with_an_extension_count(consequence, reference, candidate)
        or residue_and_stop_inserted_before_the_last_residue(consequence, reference, candidate)
    )


def residue_names(text: str) -> list[str]:
    """The three-letter residue names a run of them holds, in order."""
    return [text[i:i + 3] for i in range(0, len(text), 3)]


def stop_extension(description: str) -> tuple[str, str, str] | None:
    """The residue number, first new residue and count (``?`` for an unknown one) of an
    extension named at the stop of a protein of at least one residue; ``None`` for any other
    description."""
    m = STOP_EXTENSION_RE.fullmatch(description)
    if m is None or m.group(1) == "1" or m.group(2) == "Ter":
        return None
    return m.group(1), m.group(2), m.group(3)


def stop_replaced_by_a_stop(reference: str, candidate: str) -> bool:
    """The stop replaced by a stop against the stop written as unchanged at the same residue."""
    replaced, unchanged = STOP_REPLACED_RE.fullmatch(reference), STOP_UNCHANGED_RE.fullmatch(candidate)
    return (replaced is not None and unchanged is not None and replaced.group(2) == "Ter"
            and replaced.group(1) == unchanged.group(1) != "1")


def stop_replaced_by_first_residues(reference: str, candidate: str) -> bool:
    """The stop alone replaced by the first residues that follow the protein, closed by ``Ter``
    or not, against an extension at the same residue with the first of them: its count unknown,
    above the residues the reference names, or their number where the reference closes with no
    ``Ter``."""
    replaced, extension = STOP_REPLACED_RE.fullmatch(reference), stop_extension(candidate)
    if replaced is None or extension is None:
        return False
    number, first, count = extension
    named = residue_names(replaced.group(2))
    closed = named[-1] == "Ter"
    if closed:
        named.pop()
    if not named or "Ter" in named or replaced.group(1) != number or named[0] != first:
        return False
    return count == "?" or int(count) > len(named) or (int(count) == len(named) and not closed)


def stop_replaced_through_an_inserted_stop(reference: str, candidate: str) -> bool:
    """The stop replaced by every added residue and then a stop against an extension at the
    same residue with the first of them, its count the number of residues the reference names."""
    replaced, extension = STOP_REPLACED_RE.fullmatch(reference), stop_extension(candidate)
    if replaced is None or extension is None:
        return False
    number, first, count = extension
    named = residue_names(replaced.group(2))
    if named.pop() != "Ter" or not named or "Ter" in named:
        return False
    return replaced.group(1) == number and named[0] == first and count == str(len(named))


def residues_inserted_before_the_stop(reference: str, candidate: str) -> bool:
    """Residues inserted between the last residue and the stop against an extension at the stop
    with the first of them, its count their number."""
    inserted, extension = STOP_INSERTION_RE.fullmatch(reference), stop_extension(candidate)
    if inserted is None or extension is None:
        return False
    number, first, count = extension
    named = residue_names(inserted.group(3))
    return (int(inserted.group(1)) + 1 == int(inserted.group(2)) and inserted.group(2) == number
            and "Ter" not in named and named[0] == first and count == str(len(named)))


def lost_stop_written_as_unchanged(consequence: str, reference: str, candidate: str) -> bool:
    """The stop written as unchanged against an extension at the same residue, on a row whose
    terms hold ``stop_lost``."""
    unchanged, extension = STOP_UNCHANGED_RE.fullmatch(reference), stop_extension(candidate)
    return (unchanged is not None and extension is not None and unchanged.group(1) == extension[0]
            and "stop_lost" in consequence_set(consequence).split(","))


def range_to_the_stop_replaced(consequence: str, reference: str, candidate: str) -> bool:
    """A range from a residue to the stop replaced by residues, none of them a stop, against a
    frameshift named at that residue with the first of them, on a row whose terms hold
    ``stop_lost``. The frameshift's count is unknown, or above their number and not the one
    that leaves the new stop at the stop's own residue: a protein that keeps its length is no
    frameshift of the candidate's."""
    replaced, frameshift = STOP_RANGE_REPLACED_RE.fullmatch(reference), STOP_FRAMESHIFT_RE.fullmatch(candidate)
    if replaced is None or frameshift is None or "stop_lost" not in consequence_set(consequence).split(","):
        return False
    residue, number, stop, named = replaced.group(1), replaced.group(2), replaced.group(3), residue_names(replaced.group(4))
    count = frameshift.group(4)
    return (residue != "Ter" and "Ter" not in named and int(number) < int(stop)
            and frameshift.group(1) == residue and frameshift.group(2) == number and frameshift.group(3) == named[0]
            and (count == "?" or (int(count) > len(named) and int(count) != int(stop) - int(number) + 1)))


def stop_named_on_incomplete_coding_sequence(consequence: str, flags: str, ev: str, av: str) -> bool:
    """``hgvsp_stop_named_on_incomplete_coding_sequence`` on a differing HGVSp: the candidate
    writes none; the reference's description, an escaped ``=`` and a pair of prediction
    parentheses set aside, names the stop alone as unchanged, deleted, replaced by one residue
    or extended; and the row shows the incomplete coding sequence: its FLAGS (``flags``) name
    ``cds_end_NF``, its terms (``consequence``) hold ``incomplete_terminal_codon_variant``, or
    the description is ``Ter1=`` and its FLAGS name ``cds_start_NF``."""
    if blank(ev) or not blank(av):
        return False
    _, separator, reference = ev.partition(":p.")
    reference = reference.replace("%3D", "=")
    if reference.startswith("(") and reference.endswith(")"):
        reference = reference[1:-1]
    named = STOP_NAMED_ALONE_RE.fullmatch(reference)
    if not separator or named is None:
        return False
    return ("cds_end_NF" in flags or "incomplete_terminal_codon_variant" in consequence_set(consequence).split(",")
            or (named.group(1) == "1" and named.group(2) == "=" and "cds_start_NF" in flags))


def unknown_count_where_a_stop_is_met(consequence: str, reference: str, candidate: str) -> bool:
    """An extension with an unknown count against the same extension, at the same residue and
    with the same first residue, with a count, on a row whose terms hold ``stop_lost``."""
    unknown, counted = stop_extension(reference), stop_extension(candidate)
    return (unknown is not None and counted is not None and unknown[2] == "?" and counted[2] != "?"
            and unknown[:2] == counted[:2] and "stop_lost" in consequence_set(consequence).split(","))


def stop_deleted_with_an_extension_count(consequence: str, reference: str, candidate: str) -> bool:
    """The stop written as deleted with an extension count and no residue in the stop's place,
    against an extension at the same residue with the same count, on a row whose terms hold
    ``stop_lost``."""
    deleted, extension = STOP_DELETED_EXTENSION_RE.fullmatch(reference), stop_extension(candidate)
    return (deleted is not None and extension is not None and deleted.group(1) == extension[0]
            and deleted.group(2) == extension[2] and "stop_lost" in consequence_set(consequence).split(","))


def residue_and_stop_inserted_before_the_last_residue(consequence: str, reference: str, candidate: str) -> bool:
    """One residue and a stop inserted between two neighbouring residues, the residue the right
    flank's, against the stop written as unchanged at the residue after the right flank, on a
    row whose terms hold ``stop_retained_variant``."""
    inserted, unchanged = RESIDUE_AND_STOP_INSERTED_RE.fullmatch(reference), STOP_UNCHANGED_RE.fullmatch(candidate)
    if inserted is None or unchanged is None or STOP_RETAINED_TERM not in consequence_set(consequence).split(","):
        return False
    left, right, residue = inserted.group(1, 3, 5)
    left_number, right_number, stop = (integer_of(n, "HGVSp residue number") for n in (*inserted.group(2, 4), unchanged.group(1)))
    return "Ter" not in (left, right) and residue == right and left_number + 1 == right_number == stop - 1


def edge_row_reach(e: Entry, side: str) -> int | None:
    """How far apart, in bases, the two HGVSc notations of an allele running over a transcript's
    edge can lie: the bases of the allele's span, which bound the bases inside the transcript,
    plus one base of adjacency on either side. The span is the row's Location; without one (the
    VCF and JSON formats) and with the start open, the closed cDNA end counts the bases inside;
    ``None`` with the end open and no Location."""
    if e.location is not None:
        start, _, end = e.location.split(":", 1)[-1].partition("-")
        if start.isdigit():
            span_bases = integer_of(end, "Location position") - integer_of(start, "Location position") + 1 if end.isdigit() else 1
            return span_bases + 1
    if side == "start":
        text_position_field, _ = CDNA_START_FIELDS
        cdna_end = e.fields.get(text_position_field, "").partition("/")[0].partition("-")[2]
        cdna_end = cdna_end or e.fields.get(CDNA_END_FIELD, "")
        if cdna_end.isdigit():
            return integer_of(cdna_end, "cDNA position") + 1
    return None


def notation_coordinate_indexes(notation: re.Match) -> list[int]:
    """The ``cdna_coordinate_index`` of a ``CLIPPED_NOTATION_RE`` match's one or two
    coordinates."""
    return [cdna_coordinate_index(c) for c in (notation.group("start"), notation.group("end")) if c]


def cdna_coordinate_index(coordinate: str) -> int:
    """An HGVS cDNA coordinate as a position on one integer scale along the transcript: the 5'
    UTR ``-N`` at ``-N``, the coding sequence ``N`` at ``N - 1`` (``c.-1`` and ``c.1`` are
    neighbours), the 3' UTR ``*N`` a billion further, since its distance from the coding
    sequence is the coding length the notation does not carry, each with its intron offset added
    (``-49-1`` is one base before ``-49``, ``83+1`` one after ``83``)."""
    region, position, offset = cdna_coordinate_order(coordinate)
    if region < 0:
        return position + offset
    if region == 0:
        return position - 1 + offset
    return 10**9 + position + offset


def clipped_notation(hgvsc: str, feature: str) -> re.Match | None:
    """The ``CLIPPED_NOTATION_RE`` match of an HGVSc that names the row's transcript (its
    stable id, with or without a version), or ``None`` for an empty value, another
    transcript's notation or a notation of another shape."""
    if blank(hgvsc):
        return None
    name, _, body = hgvsc.partition(":")
    if name.split(".")[0] != feature.split(".")[0]:
        return None
    return CLIPPED_NOTATION_RE.fullmatch(body)


def copy_of_the_intron_last_base(body: str) -> bool:
    """The copy form of ``hgvsp_intron_end_insertion``: an HGVSc body duplicating the one
    position before an exon's first coding base ``N``, ``c.N-1dup``."""
    return INTRON_END_COPY_RE.fullmatch(body) is not None


def duplication_of_the_intron_last_bases(body: str) -> bool:
    """The duplication form of ``hgvsp_intron_end_insertion``: an HGVSc body duplicating two or
    more intron positions that end on the one before an exon's first coding base ``N``, its
    first position counted back from ``N`` (``c.N-K_N-1dup``, ``K`` at least 2) or, in the
    intron's first half, on from the coding base before ``N`` (``c.M+J_N-1dup``, ``M`` one
    below ``N``)."""
    m = INTRON_END_DUPLICATION_RE.fullmatch(body)
    if m is None:
        return False
    first, sign, distance, last = m.groups()
    if sign == "-":
        return first == last and int(distance) >= 2
    return int(first) == int(last) - 1


def insertion_that_cannot_end_the_intron(body: str) -> bool:
    """The insertion form of ``hgvsp_intron_end_insertion``: an HGVSc body inserting, between
    the position before an exon's first coding base ``N`` and that base (``c.N-1_Nins``), one
    base, or bases that do not end in one of ``INTRON_END_PAIRS``."""
    m = INTRON_END_INSERTION_RE.fullmatch(body)
    if m is None or m.group(1) != m.group(2):
        return False
    bases = m.group(3)
    return len(bases) == 1 or bases[-2:] not in INTRON_END_PAIRS


def insertion_that_ends_in_an_intron_end_pair(body: str) -> bool:
    """The pair form of ``hgvsp_intron_end_insertion``: an HGVSc body inserting, between the
    position before an exon's first coding base ``N`` and that base (``c.N-1_Nins``), two or
    more bases that end in one of ``INTRON_END_PAIRS``."""
    m = INTRON_END_INSERTION_RE.fullmatch(body)
    if m is None or m.group(1) != m.group(2):
        return False
    bases = m.group(3)
    return len(bases) >= 2 and bases[-2:] in INTRON_END_PAIRS


def intron_end_insertion_form(body: str) -> str | None:
    """The form of ``hgvsp_intron_end_insertion`` an HGVSc body has, each form one helper and
    one block here, or ``None`` for any other notation."""
    if copy_of_the_intron_last_base(body):
        return "copy"

    if duplication_of_the_intron_last_bases(body):
        return "duplication"

    if insertion_that_cannot_end_the_intron(body):
        return "insertion"

    if insertion_that_ends_in_an_intron_end_pair(body):
        return "pair"

    return None


def hgvsp_intron_end_insertion(e: Entry, a: Entry | None, ev: str, av: str, intron_ends: IntronEnds | None) -> bool:
    """``hgvsp_intron_end_insertion`` on a differing HGVSp: the reference's value is a protein
    notation and the candidate's is empty; the row's HGVSc is one value on both entries, a
    notation of the row's transcript (its stable id, with or without a version) of one of the
    four forms (``intron_end_insertion_form``); neither entry carries an HGVS_OFFSET; and the
    record is written in a gap where an intron of the row's transcript ends before a coding exon
    base (``IntronEnds.holds``). A record that holds the insertion with reference bases it
    shares is placed by the bases an ``ins`` notation inserts, and a text row that gives such a
    record's span with a minimised allele by those bases or by the count of positions a ``dup``
    notation duplicates (``IntronEnds.record_gaps``). A VCF or JSON entry is placed from its
    input line by the same bases or count, however the entry spells its allele
    (``notation_insertion_gaps``). The gaps come from the candidate's JSON cache: without
    ``intron_ends`` nothing is accepted."""
    if intron_ends is None or a is None or not blank(av) or ":p." not in ev:
        return False
    hgvsc = field_of(e, HGVSC_FIELDS)
    if blank(hgvsc) or field_of(a, HGVSC_FIELDS) != hgvsc:
        return False
    if not blank(field_of(e, HGVS_OFFSET_FIELDS)) or not blank(field_of(a, HGVS_OFFSET_FIELDS)):
        return False
    name, _, body = hgvsc.partition(":")
    if name.split(".")[0] != e.feature.split(".")[0] or intron_end_insertion_form(body) is None:
        return False
    inserted = notation_inserted_bases(INTRON_END_INSERTION_RE, body)
    return intron_ends.holds(e, inserted, notation_duplicated_positions(body))


def insertion_after_the_exon_last_base(body: str) -> bool:
    """The form of ``hgvsp_intron_start_insertion``: an HGVSc body inserting bases between an
    exon's last coding base ``N`` and the first position of the intron after it,
    ``c.N_N+1ins``."""
    m = INTRON_START_INSERTION_RE.fullmatch(body)
    return m is not None and m.group(1) == m.group(2)


def hgvsp_intron_start_insertion(e: Entry, a: Entry | None, ev: str, av: str, intron_ends: IntronEnds | None) -> bool:
    """``hgvsp_intron_start_insertion`` on a differing HGVSp: the reference's value is a protein
    notation and the candidate's is empty; the row's HGVSc is one value on both entries, a
    notation of the row's transcript (its stable id, with or without a version) of the class's
    form (``insertion_after_the_exon_last_base``); the two entries carry no HGVS_OFFSET, or one
    non-zero integer, the bases the 3' rule moved the record; and the gap the record is written
    in, moved by that offset as signed, is one where an exon of the row's transcript ends inside
    the coding region before an intron (``IntronEnds.holds_start``). A text row that carries the
    offset is the class as well where the gap its Location names is such a gap unmoved: a run
    under ``--shift_genomic 1`` prints the Location the 3' rule carried the record to. A VCF or
    JSON entry is placed from its input line, which is the record as written, by the bases the
    notation inserts, however the entry spells its allele (``notation_insertion_gaps``), and is
    held to the moved gap alone. A text row whose record holds the insertion with reference
    bases it shares carries no offset and is placed by the bases the notation inserts, and a
    text row that gives such a record's span with a minimised allele by those bases, turned by
    the offset where the 3' rule moved the record (``IntronEnds.record_gaps``). A run under
    ``--shift_3prime 1`` prints a moved record's allele as the rule moved it, so on a row that
    carries the offset a minimised allele is the notation's bases as they stand or turned back
    by the offset, its span's gaps read as named or moved, and the bases an ALT of an entry's
    input line inserts are the notation's as they stand or turned back by the offset. The gaps
    come from the candidate's JSON cache: without ``intron_ends`` nothing is accepted."""
    if intron_ends is None or a is None or not blank(av) or ":p." not in ev:
        return False
    hgvsc = field_of(e, HGVSC_FIELDS)
    if blank(hgvsc) or field_of(a, HGVSC_FIELDS) != hgvsc:
        return False
    offset, candidate_offset = field_of(e, HGVS_OFFSET_FIELDS), field_of(a, HGVS_OFFSET_FIELDS)
    moved = 0
    if not (blank(offset) and blank(candidate_offset)):
        if candidate_offset != offset or NON_ZERO_INTEGER_RE.fullmatch(offset) is None:
            return False
        moved = int(offset)
    name, _, body = hgvsc.partition(":")
    if name.split(".")[0] != e.feature.split(".")[0] or not insertion_after_the_exon_last_base(body):
        return False
    return intron_ends.holds_start(e, moved, notation_inserted_bases(INTRON_START_INSERTION_RE, body))


def notation_inserted_bases(pattern: re.Pattern, body: str) -> str:
    """The bases an HGVSc body inserts where it is an ``ins`` notation of ``pattern``
    (``INTRON_END_INSERTION_RE``, ``INTRON_START_INSERTION_RE``); empty for any other body."""
    m = pattern.fullmatch(body)
    return m.group(3) if m is not None else ""


def forward_readings(inserted: str, strand: str) -> tuple[str, ...]:
    """The bases a row's HGVSc inserts, on the transcript's strand, read on the forward strand
    by the row's ``strand`` (``1``, ``-1``; any other value admits either reading)."""
    if strand == "1":
        return (inserted,)
    if strand == "-1":
        return (reverse_complement(inserted),)
    return (inserted, reverse_complement(inserted))


def held_insertion_gaps(first: int, length: int, allele: str, inserted: str, strand: str) -> list[int]:
    """The gaps in which an insertion can lie that a text row's record holds with the reference
    bases it shares, each the lower of the two positions either side of the gap. The record's
    reference spans ``length`` positions from ``first`` and ``allele`` is its alternate allele
    on the forward strand, those bases with the inserted ones; ``inserted`` are the bases the
    row's HGVSc inserts, read on the forward strand by ``strand`` (``forward_readings``). The
    allele holds as many bases more than the span as the notation inserts. For each count of
    span bases clipped from the allele's front, the rest clipped from its back, the bases
    between must be the inserted ones, and the gap's lower position is ``first`` plus that
    count less one. A text row states no reference bases, so every clipping that leaves the
    inserted bases names a gap."""
    if not inserted or length < 1 or len(allele) - length != len(inserted):
        return []
    readings = forward_readings(inserted, strand)
    gaps = []
    for front in range(length + 1):
        back = len(allele) - (length - front)
        if allele[front:back] in readings:
            gaps.append(first + front - 1)
    return gaps


def clipped_insertions(first: int, reference: str, alternate: str) -> list[tuple[int, str]]:
    """Every pure insertion an alternate allele clips to against its record's reference allele,
    which starts at ``first``: the lower of the two positions either side of the gap, and the
    bases inserted there. For each count of reference bases clipped from the alternate's front,
    the rest clipped from its back, the clipped bases must be the reference's own and the bases
    between are the inserted ones. An empty reference allele, an insertion written alone, names
    the gap before ``first``."""
    extra = len(alternate) - len(reference)
    if extra < 1:
        return []
    insertions = []
    for front in range(len(reference) + 1):
        back = front + extra
        if alternate[:front] == reference[:front] and alternate[back:] == reference[front:]:
            insertions.append((first + front - 1, alternate[front:back]))
    return insertions


def notation_insertion_gaps(
    first: int, reference: str, alternates: list[str], inserted: str, strand: str, duplicated: int = 0, moved: int = 0
) -> list[int]:
    """The gaps an input line's alleles give the insertion a row's HGVSc describes, each the
    lower of the two positions either side of the gap. ``reference`` is the record's reference
    allele from ``first`` and ``alternates`` its alternate alleles, on the forward strand; an
    alternate that clips against the reference to a pure insertion names a gap with its
    inserted bases (``clipped_insertions``), and the notation claims the gap whose bases are
    its own: the bases an ``ins`` notation inserts (``inserted``), read on the forward strand
    by ``strand`` (``forward_readings``), as they stand or, on a record the 3' rule moved
    ``moved`` positions, turned back by the offset (``turned_back``); or, beside a ``dup``
    notation that states its count, as many bases as the positions duplicated (``duplicated``)
    on an unmoved record. A ``dup`` notation that states no count names neither, and claims
    only an insertion written alone, whose bases follow a reference allele of one base or none.
    How an entry spells its own allele plays no part."""
    wanted = {bases for reading in forward_readings(inserted, strand) for bases in (reading, turned_back(reading, moved))}
    gaps = []
    for alternate in alternates:
        for gap, bases in clipped_insertions(first, reference, alternate):
            if inserted:
                claimed = bases in wanted
            elif duplicated > 0:
                claimed = moved == 0 and len(bases) == duplicated
            else:
                claimed = len(reference) < 2 and gap == first + len(reference) - 1
            if claimed:
                gaps.append(gap)
    return gaps


def notation_duplicated_positions(body: str) -> int:
    """The count of positions an HGVSc body duplicates where the notation states it: one for a
    copy of the intron's last base (``c.N-1dup``), ``K`` for a duplication counted back from
    the exon's first coding base (``c.N-K_N-1dup``). 0 for any other body, a duplication that
    starts in the intron's first half among them, whose count depends on the intron's length."""
    if INTRON_END_COPY_RE.fullmatch(body) is not None:
        return 1
    m = INTRON_END_DUPLICATION_RE.fullmatch(body)
    if m is not None and m.group(2) == "-" and m.group(1) == m.group(4):
        return int(m.group(3))
    return 0


def lone_insertion_gaps(first: int, length: int, allele: str, inserted: str, duplicated: int = 0, moved: int = 0) -> list[int]:
    """The gaps a text row's Location names where its Allele is the insertion alone, each the
    lower of the two positions either side of the gap: every gap inside and beside the span of
    ``length`` positions from ``first``, from ``first`` less one to the span's last position. A
    record that carries another allele beside the insertion keeps its reference span as the
    row's Location, and a minimised allele spells the inserted bases alone. ``allele`` is the
    insertion the row's HGVSc describes: the bases an ``ins`` notation inserts (``inserted``),
    on either strand, as they stand or turned back by the bases the 3' rule moved the record
    (``moved``), since a run prints a moved record's allele as the rule moved it or as the
    record writes it; or, for a ``dup`` notation, which names no bases, as many bases as the
    positions it duplicates (``duplicated``) on an unmoved record. None for any other allele."""
    if length < 1 or any(ord(base) not in COMPLEMENT for base in allele):
        return []
    if inserted:
        turns = dict.fromkeys((0, abs(moved) % len(inserted)))
        readings = (allele, reverse_complement(allele))
        if all(bases[turn:] + bases[:turn] != inserted for bases in readings for turn in turns):
            return []
    elif duplicated < 1 or moved != 0 or len(allele) != duplicated:
        return []
    return list(range(first - 1, first + length))


def turned_back(allele: str, moved: int) -> str:
    """An insertion's bases on the forward strand as the record writes them, from ``allele``,
    the bases as the 3' rule moved them ``moved`` positions along the chromosome: a run that
    prints the moved allele turns the inserted bases as it carries them, their first base to the
    end for each position moved up the chromosome. The allele itself for an unmoved record."""
    back = len(allele) - moved % len(allele) if allele else 0
    return allele[back:] + allele[:back]


def cdna_genomic_position(exons: list[tuple[int, int]], forward: bool, cdna: int) -> int | None:
    """The genomic position of the cDNA position ``cdna``, counted along ``exons`` (in genomic
    order) in the transcript's direction; ``None`` for a position outside them."""
    remaining = cdna
    for start, end in exons if forward else reversed(exons):
        length = end - start + 1
        if 1 <= remaining <= length:
            return start + remaining - 1 if forward else end - remaining + 1
        remaining -= length
    return None


def coding_transcript_exons(transcript: dict) -> tuple[bool, list[tuple[int, int]], int, int] | None:
    """Whether a cache transcript lies on the forward strand, its exons in genomic order, and
    the lower and the higher genomic bound of its coding region, stop codon included: the
    genomic positions of its ``cdna_coding_start`` and ``cdna_coding_end``
    (``cdna_genomic_position``). ``None`` for a transcript without a coding region and for one
    with a coordinate that is no integer."""
    try:
        forward = int(transcript.get("strand", 1)) >= 0
        exons = sorted((int(exon["start"]), int(exon["end"])) for exon in transcript.get("exons") or [])
        coding = [int(transcript[key]) for key in ("cdna_coding_start", "cdna_coding_end")]
    except (KeyError, TypeError, ValueError):
        return None
    bounds = [cdna_genomic_position(exons, forward, cdna) for cdna in coding]
    if None in bounds:
        return None
    return forward, exons, min(bounds), max(bounds)


def transcript_intron_end_gaps(transcript: dict) -> set[int]:
    """The gaps of one cache transcript in which an intron of more than ``SHORT_INTRON_BASES``
    ends before an exon base inside the coding region, each the lower of the two genomic
    positions either side of the gap: on the forward strand the intron's last base, the exon's
    first base following it; on the reverse strand the exon's first base in the transcript's
    direction, its highest genomic base, the intron's last base following it. Empty where
    ``coding_transcript_exons`` reads no coding region."""
    coding = coding_transcript_exons(transcript)
    if coding is None:
        return set()
    forward, exons, low, high = coding
    gaps: set[int] = set()
    for index in range(1, len(exons)):
        lower_end, upper_start = exons[index - 1][1], exons[index][0]
        if upper_start - lower_end - 1 <= SHORT_INTRON_BASES:
            continue
        exon_first = upper_start if forward else lower_end
        if low <= exon_first <= high:
            gaps.add(upper_start - 1 if forward else lower_end)
    return gaps


def transcript_intron_start_gaps(transcript: dict) -> set[int]:
    """The gaps of one cache transcript in which an exon ends inside the coding region before an
    intron of any length, each the lower of the two genomic positions either side of the gap: on
    the forward strand the exon's last base, the intron's first base following it; on the
    reverse strand the intron's first base in the transcript's direction, its highest genomic
    base, the exon's last base in that direction following it. Empty where
    ``coding_transcript_exons`` reads no coding region."""
    coding = coding_transcript_exons(transcript)
    if coding is None:
        return set()
    forward, exons, low, high = coding
    gaps: set[int] = set()
    for index in range(1, len(exons)):
        lower_end, upper_start = exons[index - 1][1], exons[index][0]
        exon_last = lower_end if forward else upper_start
        if upper_start - lower_end > 1 and low <= exon_last <= high:
            gaps.add(lower_end if forward else upper_start - 1)
    return gaps


class SpliceWindows:
    """The introns of a candidate JSON cache's transcripts, by transcript stable id, each
    transcript's in its own order with its strand: ``holds`` reads from them whether an
    entry's record lies in the two canonical windows ``maxentscan_overlapping_windows`` names.
    A chromosome's shards are read on the first lookup of an entry whose record names its
    directory, as written or without a leading ``chr``; a record on a chromosome the cache
    holds no directory of lies in no window. A stable id a chromosome's shards give two
    different intron lists is left out.
    ``input_lines`` are the input's data lines, VCF or Ensembl format, which place a VCF or
    JSON entry."""

    def __init__(self, cache: Path, input_lines: list[str] | None = None) -> None:
        self.root = cache / "transcripts"
        if not self.root.is_dir():
            raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
        self.input_lines = input_lines or []
        self._chromosomes = {directory.name for directory in self.root.iterdir() if directory.is_dir()}
        self._introns: dict[str, dict[str, tuple[int, tuple[tuple[int, int], ...]]]] = {}

    def _read(self, chromosome: str) -> dict[str, tuple[int, tuple[tuple[int, int], ...]]]:
        introns: dict[str, tuple[int, tuple[tuple[int, int], ...]]] = {}
        conflicting: set[str] = set()
        directory = self.root / chromosome
        for shard in sorted(directory.glob("*.json")) + sorted(directory.glob("*.json.gz")):
            try:
                parsed = json.loads(read_text(shard))
            except ValueError as exc:
                raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
            for t in NearestSites.shard_transcripts(shard, parsed):
                stable_id, exons, strand = t.get("stable_id"), t.get("exons"), cache_integer(t.get("strand"))
                if not isinstance(stable_id, str) or not isinstance(exons, list) or strand not in (1, -1):
                    continue
                if not all(isinstance(exon, dict) for exon in exons):
                    continue
                bounds = [(cache_integer(exon.get("start")), cache_integer(exon.get("end"))) for exon in exons]
                if any(None in pair for pair in bounds):
                    continue
                spans = sorted(bounds)
                gaps = [(spans[i][1] + 1, spans[i + 1][0] - 1) for i in range(len(spans) - 1)]
                read = (strand, tuple(gaps if strand == 1 else reversed(gaps)))
                if introns.setdefault(stable_id, read) != read:
                    conflicting.add(stable_id)
        for stable_id in conflicting:
            del introns[stable_id]
        return introns

    def place(self, e: Entry) -> tuple[str, int] | None:
        """The chromosome and position of an entry's record where the record is one base: a
        text entry's Location that names one position, and the input line of a VCF or JSON
        entry whose reference allele is one base. ``None`` for any other record."""
        if e.location is not None:
            chromosome, _, span = e.location.rpartition(":")
            return (chromosome, int(span)) if chromosome and span.isdigit() else None
        ordinal = e.record if e.record is not None else e.line
        if ordinal is None or ordinal >= len(self.input_lines):
            return None
        record = input_line_record(self.input_lines[ordinal])
        if record is None or len(record[2]) != 1:
            return None
        return record[0], record[1]

    def holds(self, e: Entry) -> bool:
        """Whether the entry's record is one base that lies in the acceptor window of one
        intron of the entry's transcript and in the donor window of the next intron in the
        transcript's order, the windows as MaxEntScan.pm:282-291 places them: the donor's
        from three bases before an intron to its sixth base, the acceptor's from an intron's
        twentieth base from its end to three bases after it. A Feature carrying a version the
        cache's stable id lacks is looked up without it."""
        placed = self.place(e)
        if placed is None:
            return False
        chromosome, position = placed
        directory = next((name for name in (chromosome, chromosome.removeprefix("chr")) if name in self._chromosomes), None)
        if directory is None:
            return False
        if directory not in self._introns:
            self._introns[directory] = self._read(directory)
        introns = self._introns[directory]
        read = introns.get(e.feature, introns.get(re.sub(r"\.\d+$", "", e.feature)))
        if read is None:
            return False
        strand, gaps = read
        for index in range(len(gaps) - 1):
            (upstream_start, upstream_end), (downstream_start, downstream_end) = gaps[index], gaps[index + 1]
            if strand == 1:
                acceptor, donor = (upstream_end - 19, upstream_end + 3), (downstream_start - 3, downstream_start + 5)
            else:
                acceptor, donor = (upstream_start - 3, upstream_start + 19), (downstream_end - 5, downstream_end + 3)
            if acceptor[0] <= position <= acceptor[1] and donor[0] <= position <= donor[1]:
                return True
        return False


class IntronEnds:
    """The gaps of a candidate JSON cache at the two ends of an intron, by chromosome and
    transcript: those in which an intron ends before a coding exon base
    (``transcript_intron_end_gaps``), which ``hgvsp_intron_end_insertion`` reads, and those in
    which a coding exon ends before an intron (``transcript_intron_start_gaps``), which
    ``hgvsp_intron_start_insertion`` reads. A cache files its transcripts under
    ``transcripts/<chromosome>/<first>-<last>.json`` or ``.json.gz``; a shard is read when a
    row asks for a transcript no shard read so far holds, the shard whose window holds the
    row's position first and the chromosome's others nearest first, and no shard is read
    twice. ``input_lines`` are the input's data lines, VCF or Ensembl format, which place a VCF
    or JSON entry."""

    def __init__(self, cache: Path, input_lines: list[str] | None = None) -> None:
        self.transcripts = cache / "transcripts"
        if not self.transcripts.is_dir():
            raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
        self.input_lines = input_lines or []
        # Per chromosome the unread shards, each with its window; per chromosome and stable
        # id the gaps of every transcript read, those at an intron's end and those at its start.
        self._unread: dict[str, list[tuple[int, int, Path]]] = {}
        self._gaps: dict[tuple[str, str], tuple[set[int], set[int]]] = {}

    def holds(self, e: Entry, inserted: str = "", duplicated: int = 0) -> bool:
        """Whether one of the gaps the entry's insertion can be written in (``record_gaps``, by
        ``inserted``, the bases the row's HGVSc inserts, or ``duplicated``, the count of
        positions it duplicates) is a gap of the entry's transcript in which an intron ends
        before a coding exon base."""
        gaps = self.record_gaps(e, inserted, duplicated)
        return any(position in self.gaps(chromosome, e.feature, position) for chromosome, position in gaps)

    def holds_start(self, e: Entry, moved: int = 0, inserted: str = "") -> bool:
        """Whether one of the gaps the entry's insertion can be written in (``record_gaps``, by
        ``inserted``, the bases the row's HGVSc inserts), moved ``moved`` positions along the
        chromosome, is a gap of the entry's transcript in which a coding exon ends before an
        intron. A text entry's gap is read as named as well: a run that prints the moved
        location names the gap the insertion is left in."""
        gaps = self.record_gaps(e, inserted, moved=moved)
        moves = (moved,) if e.location is None else tuple(dict.fromkeys((moved, 0)))
        return any(position + move in self.start_gaps(chromosome, e.feature, position + move)
                   for chromosome, position in gaps for move in moves)

    def record_gaps(self, e: Entry, inserted: str = "", duplicated: int = 0, moved: int = 0) -> list[tuple[str, int]]:
        """The chromosome and lower position of each gap the entry's insertion can be written
        in. A VCF or JSON entry is placed from its input line by the insertion its HGVSc
        describes (``notation_insertion_gaps``), however the entry spells its allele: an ALT of
        the line, the ALT less a base every allele of the record shares, a minimised allele or
        the allele as the 3' rule moved it. A text entry is placed by its Location: the gap of
        an insertion written alone first (``written_gap``); a record that carries another
        allele beside the insertion keeps the reference base they share and gives its reference
        span, one position or more, where an unmoved allele that holds the span's bases with
        the inserted ones names the gaps its clippings leave (``held_insertion_gaps``), as
        written or reverse-complemented, since the row spells it on the record's strand, and an
        allele that is the insertion alone, a minimised allele, names every gap of the span
        (``lone_insertion_gaps``), a span of two adjacent positions among them, which is an
        insertion written alone or a reference allele of two bases. ``inserted`` are the bases
        the row's HGVSc inserts, ``duplicated`` the count of positions it duplicates and
        ``moved`` the positions the 3' rule moved the record, by which a moved allele is read;
        the gaps are the record's own."""
        strand = field_of(e, STRAND_FIELDS)
        if e.location is None:
            record = self.input_record(e)
            if record is None:
                return []
            chromosome, first, reference, alternates = record
            positions = notation_insertion_gaps(first, reference, alternates, inserted, strand, duplicated, moved)
            return [(chromosome, position) for position in dict.fromkeys(positions)]
        gap = self.written_gap(e)
        gaps = [gap] if gap is not None else []
        allele = e.allele.upper()
        chromosome, _, span = e.location.rpartition(":")
        start, _, end = span.partition("-")
        end = end or start
        if re.fullmatch(r"[A-Z]+", allele) is None or not (chromosome and start.isdigit() and end.isdigit()):
            return gaps
        first, length = int(start), int(end) - int(start) + 1
        positions = lone_insertion_gaps(first, length, allele, inserted, duplicated, moved)
        if moved == 0:
            for bases in (allele, reverse_complement(allele)):
                positions += held_insertion_gaps(first, length, bases, inserted, strand)
        return gaps + [(chromosome, position) for position in dict.fromkeys(positions)]

    def input_record(self, e: Entry) -> tuple[str, int, str, list[str]] | None:
        """The chromosome, first position, reference allele and alternate alleles of a VCF or
        JSON entry's input line (``input_line_record``): the line of the entry's record, or the
        Ensembl-format line a VCF entry's own line is written from (``Entry.line``). ``None``
        for an entry without one and for a line of neither format."""
        ordinal = e.record if e.record is not None else e.line
        if ordinal is None or ordinal >= len(self.input_lines):
            return None
        return input_line_record(self.input_lines[ordinal])

    def written_gap(self, e: Entry) -> tuple[str, int] | None:
        """The chromosome and the lower of the two positions either side of the gap a text
        entry's record is written in: the first position of its Location where that names two
        adjacent positions. ``None`` for any other Location and for a VCF or JSON entry, which
        is placed from its input line."""
        if e.location is None:
            return None
        chromosome, _, span = e.location.rpartition(":")
        start, _, end = span.partition("-")
        if not (chromosome and start.isdigit() and end.isdigit() and int(end) == int(start) + 1):
            return None
        return chromosome, int(start)

    def gaps(self, chromosome: str, stable_id: str, position: int) -> Collection[int]:
        """The gaps of the transcript ``stable_id`` on ``chromosome`` in which an intron ends
        before a coding exon base (``transcript_gaps``)."""
        return self.transcript_gaps(chromosome, stable_id, position)[0]

    def start_gaps(self, chromosome: str, stable_id: str, position: int) -> Collection[int]:
        """The gaps of the transcript ``stable_id`` on ``chromosome`` in which a coding exon ends
        before an intron (``transcript_gaps``)."""
        return self.transcript_gaps(chromosome, stable_id, position)[1]

    def transcript_gaps(self, chromosome: str, stable_id: str, position: int) -> tuple[Collection[int], Collection[int]]:
        """The gaps of the transcript ``stable_id`` on ``chromosome``, those at an intron's end
        and those at an intron's start, none for a transcript the cache does not hold there:
        the transcript the cache files under that id or, where it files none, under the id
        without a trailing version (a row under ``--transcript_version``). The chromosome's
        directory is the one of that name or, failing that, of the name without a leading
        ``chr``; its unread shards are read nearest to ``position`` first until one holds the
        transcript."""
        named = (self.transcripts / chromosome).is_dir()
        directory = chromosome if named else chromosome.removeprefix("chr")
        if directory not in self._unread:
            self._unread[directory] = self.shards(directory)
        unread = self._unread[directory]
        key = (directory, stable_id)
        while key not in self._gaps and unread:
            nearest = min(unread, key=lambda shard: max(shard[0] - position, position - shard[1], 0))
            unread.remove(nearest)
            self.read(directory, nearest[2])
        if key not in self._gaps:
            key = (directory, VERSION_SUFFIX_RE.sub("", stable_id))
        return self._gaps.get(key, ((), ()))

    def shards(self, directory: str) -> list[tuple[int, int, Path]]:
        """The shards of the chromosome directory ``directory`` with their windows; none for a
        chromosome the cache does not hold."""
        path = self.transcripts / directory
        if not directory or not path.is_dir():
            return []
        named = ((CACHE_SHARD_NAME_RE.fullmatch(shard.name), shard) for shard in sorted(path.iterdir()))
        return [(int(m.group(1)), int(m.group(2)), shard) for m, shard in named if m is not None]

    def read(self, chromosome: str, shard: Path) -> None:
        """Reads one shard: the gaps of each of its transcripts, kept under ``chromosome`` and
        the transcript's stable id."""
        try:
            parsed = json.loads(read_text(shard))
        except ValueError as exc:
            raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
        for transcript in NearestSites.shard_transcripts(shard, parsed):
            stable_id = str(transcript.get("stable_id") or "")
            if stable_id:
                ends, starts = self._gaps.setdefault((chromosome, stable_id), (set(), set()))
                ends.update(transcript_intron_end_gaps(transcript))
                starts.update(transcript_intron_start_gaps(transcript))


def span_matches_cdna_position(e: Entry, numbering: str, start: str, end: str) -> bool:
    """Whether a candidate notation's coordinates agree with the row's cDNA position under a
    run that prints ``SHIFT_LENGTH``, where that position is the shifted span's. In ``n.``
    numbering, which is cDNA numbering, a plain coordinate on both sides must be equal; in
    ``c.`` numbering two plain coordinates of one region (the 5' UTR, the coding sequence or
    the 3' UTR) span as many bases as the row's two plain sides. A side the row leaves open
    (``?``), a coordinate with an intron offset and a ``c.`` span across a region boundary say
    nothing, as does a run without ``SHIFT_LENGTH``, which positions the row at the written
    span. The text formats write both sides in one ``cDNA_position`` (``23-43``, ``?-31``,
    ``23`` for one base, ``23-43/1500`` under --total_length); the JSON format writes
    ``cdna_start`` and ``cdna_end`` apart and leaves an open side out (JSON.pm:340-341)."""
    if all(f not in e.fields for f in SHIFT_LENGTH_FIELDS):
        return True
    text_position_field, json_start_field = CDNA_START_FIELDS
    if text_position_field in e.fields:
        # Under --total_length the position carries the transcript's length after a slash.
        row_start, separator, row_end = e.fields[text_position_field].partition("/")[0].partition("-")
        if not separator:
            row_end = row_start
    else:
        row_start = e.fields.get(json_start_field, "")
        row_end = e.fields.get(CDNA_END_FIELD, "")
    if numbering == "n":
        for row_position, coordinate in ((row_start, start), (row_end, end)):
            if row_position.isdigit() and coordinate.isdigit() and row_position != coordinate:
                return False
        return True
    if not (row_start.isdigit() and row_end.isdigit()):
        return True
    start_region, start_position, start_offset = cdna_coordinate_order(start)
    end_region, end_position, end_offset = cdna_coordinate_order(end)
    if start_offset or end_offset or start_region != end_region:
        return True
    return end_position - start_position == int(row_end) - int(row_start)


def acceptor_site_offsets(body: str) -> tuple[str, ...] | None:
    """The intron offsets of the coordinates that place an insertion's or deletion's shifted
    span, where an HGVSc body places it at an acceptor site's two bases or in the exon,
    touching no intron base away from the splice site: every coordinate is exonic or carries
    the intron offset -1 or -2, a duplication's first coordinate excepted, since the
    duplicated bases run back into the intron while the insertion sits after the last. An
    empty tuple is a notation whose coordinates are all exonic; ``None`` a coordinate at
    offset -3 or beyond, or on the donor side, or any other notation."""
    m = INDEL_NOTATION_RE.fullmatch(body)
    if m is None:
        return None
    first, last, kind = m.groups()
    if last is None:
        coordinates = [first]
    elif kind == "dup":
        coordinates = [last]
    else:
        coordinates = [first, last]
    offsets = []
    for coordinate in coordinates:
        offset = INTRON_OFFSET_RE.search(coordinate)
        if offset is not None:
            if offset.group() not in ACCEPTOR_OFFSETS:
                return None
            offsets.append(offset.group())
    return tuple(offsets)


def cdna_step(position: int, step: int) -> int:
    """The cDNA coordinate one base (``step`` 1 or -1) from ``position``: the ``c.`` numbering
    has no 0, so -1 is followed by 1."""
    moved = position + step
    return moved + step if moved == 0 else moved


def clamped_insertion_position(body: str) -> tuple[str, int, str, int] | None:
    """The numbering, the position P, the inserted bases and the count of bases stripped from
    them of the notation Ensembl VEP writes for an insertion clamped to one base of the
    transcript slice (``hgvs_shifted_insertion_clamped_end``): ``P X>Y`` and ``P delins``
    carry the bases whole, ``P-1_P ins`` all but the last; ``None`` for any other notation."""
    m = CLAMPED_NOTATION_RE.fullmatch(body)
    if m is not None:
        return m.group(1), int(m.group(2)), m.group(3) or m.group(4), 0
    m = INSERTION_NOTATION_RE.fullmatch(body)
    if m is not None and cdna_step(int(m.group(3)), -1) == int(m.group(2)):
        return m.group(1), int(m.group(3)), m.group(4), 1
    return None


def shifted_insertion_end(body: str) -> tuple[str, int, str | None, int] | None:
    """The numbering, the position P after the notation of an insertion at its shifted span,
    its inserted bases and their count: a duplication ending at P-1, whose notation carries
    no bases (``None``) and whose span is their count (a span crossing from c.-1 to c.1 holds
    one base fewer than its bounds say, since the ``c.`` numbering has no 0), or an insertion
    between P-1 and P; ``None`` for any other notation."""
    m = SHIFTED_DUPLICATION_RE.fullmatch(body)
    if m is not None:
        end = int(m.group(3))
        start = int(m.group(2)) if m.group(2) else end
        duplicated_count = end - start + 1
        if start < 0 < end:
            duplicated_count -= 1
        return m.group(1), cdna_step(end, 1), None, duplicated_count
    m = INSERTION_NOTATION_RE.fullmatch(body)
    if m is not None and cdna_step(int(m.group(3)), -1) == int(m.group(2)):
        return m.group(1), int(m.group(3)), m.group(4), len(m.group(4))
    return None


def deleted_span_term_pair(e: Entry, a: Entry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key have the
    shape of ``refseq_edited_model_deleted_span_term``: a Transcript entry of a sequence allele
    (bases, or ``-`` for a deletion's), ``BAM_EDIT`` ``OK`` on both, the reference's terms
    ``intergenic_variant`` alone and the candidate's exactly one of
    ``DELETED_SPAN_REGION_TERMS``, no cDNA position on either and one IMPACT on both."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript":
        return False
    if not re.fullmatch(r"[ACGTN-]+", e.allele):
        return False
    if consequence_set(e.consequence) != INTERGENIC_TERM or consequence_set(a.consequence) not in DELETED_SPAN_REGION_TERMS:
        return False
    for entry in (e, a):
        if not any(entry.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS):
            return False
        if any(not blank(entry.fields.get(f, ABSENT)) for f in CDNA_POSITION_FIELDS):
            return False
    expected_impact = next((e.fields[f] for f in IMPACT_FIELDS if f in e.fields), ABSENT)
    actual_impact = next((a.fields[f] for f in IMPACT_FIELDS if f in a.fields), ABSENT)
    return expected_impact == actual_impact


def tail_term_pair(e: Entry, a: Entry) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key have the
    shape of ``refseq_edited_model_tail_term``: a Transcript entry of a sequence allele,
    ``BAM_EDIT`` ``OK`` on both, the reference's terms ``intergenic_variant`` alone and the
    candidate's ``3_prime_UTR_variant`` alone, a cDNA position on both and one IMPACT on both."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript":
        return False
    if not re.fullmatch(r"[ACGTN-]+", e.allele):
        return False
    if consequence_set(e.consequence) != INTERGENIC_TERM or consequence_set(a.consequence) != THREE_PRIME_UTR_TERM:
        return False
    for entry in (e, a):
        if not any(entry.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS):
            return False
        if all(blank(entry.fields.get(f, ABSENT)) for f in CDNA_POSITION_FIELDS):
            return False
    expected_impact = next((e.fields[f] for f in IMPACT_FIELDS if f in e.fields), ABSENT)
    actual_impact = next((a.fields[f] for f in IMPACT_FIELDS if f in a.fields), ABSENT)
    return expected_impact == actual_impact


def cds_bounds(e: Entry) -> tuple[int, int] | None:
    """The first and last coding base of the entry's span: the ``CDS_position`` of the text
    formats (one position or a range) or the JSON ``cds_start`` and ``cds_end``, each bound read
    up to the ``/total`` suffix ``--total_length`` appends; ``None`` where a bound is missing or
    is not a number (``?``)."""
    position = e.fields.get("CDS_position")
    if position is not None:
        start, _, end = position.partition("/")[0].partition("-")
        end = end or start
    else:
        start = e.fields.get("cds_start", ABSENT).partition("/")[0]
        end = e.fields.get("cds_end", ABSENT).partition("/")[0]
    if not start.isdigit() or not end.isdigit():
        return None
    return int(start), int(end)


def cds_span_length(e: Entry) -> int | None:
    """The number of coding bases the entry's span holds (``cds_bounds``), ``None`` where a
    bound is missing or is not a number."""
    bounds = cds_bounds(e)
    return None if bounds is None else bounds[1] - bounds[0] + 1


def deleted_span_codons_pair(e: Entry, a: Entry, expected_value: str, actual_value: str) -> bool:
    """Whether the Codons ``expected_value`` of the reference entry ``e`` and ``actual_value``
    of the candidate entry ``a`` have the shape of ``refseq_edited_model_deleted_span_codons``:
    a Transcript entry of a sequence allele carrying ``BAM_EDIT`` ``OK`` on both whose USED_REF
    has more bases than its CDS span, the two values' alternate codons equal, their reference
    codons of one length agreeing at every lower-case base and differing in an upper-case one,
    the reference's upper-case bases one run that ends its codons and spells the leading bases
    of USED_REF read on the transcript's strand, as the splice leaves them."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript":
        return False
    if not re.fullmatch(r"[ACGTN-]+", e.allele) or blank(expected_value) or blank(actual_value):
        return False
    for entry in (e, a):
        if not any(entry.fields.get(f) == "OK" for f in BAM_EDIT_FIELDS):
            return False
    used_ref = next((e.fields[f] for f in USED_REF_FIELDS if f in e.fields), ABSENT)
    span = cds_span_length(e)
    if blank(used_ref) or span is None or len(used_ref) <= span:
        return False
    expected_reference_codons, _, expected_alternate_codons = expected_value.partition("/")
    actual_reference_codons, _, actual_alternate_codons = actual_value.partition("/")
    if expected_alternate_codons != actual_alternate_codons:
        return False
    if (
        len(expected_reference_codons) != len(actual_reference_codons)
        or expected_reference_codons == actual_reference_codons
    ):
        return False
    if not all(
        x == y for x, y in zip(expected_reference_codons, actual_reference_codons) if x.islower() or y.islower()
    ):
        return False
    # The splice writes the allele from the variant's first base to the window's end, so the
    # reference's upper-case bases are one run that ends the window and spells the allele's
    # leading bases on the transcript's strand.
    upper_start = next((i for i, x in enumerate(expected_reference_codons) if x.isupper()), None)
    if upper_start is None:
        return False
    strand = next((e.fields[f] for f in STRAND_FIELDS if f in e.fields), ABSENT)
    allele = reverse_complement(used_ref) if strand == "-1" else used_ref
    # A lower-case base inside the run fails the prefix test, since the allele is upper-cased.
    return allele.upper().startswith(expected_reference_codons[upper_start:])


def unshifted_span_term_pair(e: Entry, a: Entry, reference_term: str, candidate_term: str, either_way: bool) -> bool:
    """Whether the reference entry ``e`` and the candidate entry ``a`` of one key have the
    shape of a class whose one test Ensembl VEP reads at the unshifted span
    (``mature_mirna_unshifted_span``, ``stop_codon_unshifted_span``): a Transcript entry with a
    positive ``SHIFT_LENGTH``, which only an insertion or deletion is given, whose two
    Consequence sets are equal but for ``reference_term`` on the reference's side against
    ``candidate_term`` on the candidate's, or, under ``either_way``, the two terms swapped. A
    set differing by any other term, or by these two beside another, stays a difference."""
    if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "transcript" or not positive_shift(e):
        return False
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    reference_only = expected_terms - actual_terms
    candidate_only = actual_terms - expected_terms
    if reference_only == {reference_term} and candidate_only == {candidate_term}:
        return True
    return either_way and reference_only == {candidate_term} and candidate_only == {reference_term}


def mature_mirna_unshifted_span_pair(e: Entry, a: Entry) -> bool:
    """``mature_mirna_unshifted_span``: ``mature_miRNA_variant`` against
    ``non_coding_transcript_exon_variant`` in either direction."""
    return unshifted_span_term_pair(e, a, MATURE_MIRNA_TERM, NON_CODING_EXON_TERM, either_way=True)


def stop_codon_unshifted_span_pair(e: Entry, a: Entry) -> bool:
    """``stop_codon_unshifted_span``: the reference's ``coding_sequence_variant`` against the
    candidate's ``stop_retained_variant``."""
    return unshifted_span_term_pair(e, a, CODING_SEQUENCE_TERM, STOP_RETAINED_TERM, either_way=False)


def stop_codon_test_coordinates(e: Entry) -> bool:
    """``compare_vep_outputs.stop_codon_test_coordinates``, which defines the test for both comparators."""
    return _stop_codon_test_coordinates(e)


def stop_codon_test_not_made(e: Entry) -> bool:
    """``compare_vep_outputs.stop_codon_test_not_made``, which defines the test for both comparators."""
    return _stop_codon_test_not_made(e)


def stop_retained_without_stop_codon_test_pair(e: Entry, a: Entry) -> bool:
    """``compare_vep_outputs.stop_retained_without_stop_codon_test_pair``, which defines the class for both comparators."""
    return _stop_retained_without_stop_codon_test_pair(e, a)


def stop_codon_overlap_genomic_span_pair(e: Entry, a: Entry) -> bool:
    """``compare_vep_outputs.stop_codon_overlap_genomic_span_pair``, which defines the class for both comparators."""
    return _stop_codon_overlap_genomic_span_pair(e, a)


def stop_codon_overlap_insertion_range_pair(e: Entry, a: Entry) -> bool:
    """``compare_vep_outputs.stop_codon_overlap_insertion_range_pair``, which defines the class for both comparators."""
    return _stop_codon_overlap_insertion_range_pair(e, a)


def stop_retained_three_base_stop_codon_test_pair(e: Entry, a: Entry) -> bool:
    """``compare_vep_outputs.stop_retained_three_base_stop_codon_test_pair``, which defines the class for both comparators."""
    return _stop_retained_three_base_stop_codon_test_pair(e, a)


def stop_codon_removed_without_stop_codon_test_pair(e: Entry, a: Entry) -> bool:
    """``compare_vep_outputs.stop_codon_removed_without_stop_codon_test_pair``, which defines the class for both comparators."""
    return _stop_codon_removed_without_stop_codon_test_pair(e, a)


def at_its_location(e: Entry) -> Entry:
    """A VCF or JSON entry with the Location its record gives its allele
    (``Entry.allele_location``) as its ``location``, the entry a predicate over a default
    row's Location reads; any other entry as it is."""
    if e.location is not None or e.allele_location is None:
        return e
    return replace(e, location=e.allele_location)


def splice_term_past_the_record_span_pair(e: Entry, a: Entry) -> bool:
    """``compare_vep_outputs.splice_term_past_the_record_span_pair``, which defines the class
    for both comparators, on the two entries at their Locations (``at_its_location``): a VCF or
    JSON entry whose record gives its allele no Location is refused."""
    return _splice_term_past_the_record_span_pair(at_its_location(e), at_its_location(a))


def read_entry_references(expected: list[Entry], actual: list[Entry], path: Path) -> None:
    """Gives the default and tab entries of a key whose two Consequence texts differ, where
    they carry no reference bases (``Entry.reference``), the bases their allele stands in
    for, read by position from the records of the input VCF ``path`` as the tuple comparator
    reads them (``compare_vep_outputs.read_record_references``). An entry no record states,
    or two records read differently, keeps ``None``, and an entry whose Location names no
    span of whole positions (``location_span``), which no record gives, asks for none."""
    actual_by_key = {a.key: a for a in actual}
    asking: list[Entry] = []
    for e in expected:
        a = actual_by_key.get(e.key)
        if a is not None and e.consequence != a.consequence:
            asking += (entry for entry in (e, a) if location_span(entry) is not None)
    try:
        read_record_references(path, asking)
    except (OSError, EOFError, zlib.error) as exc:
        raise ComparatorError(f"cannot read {path}: {getattr(exc, 'strerror', None) or exc}") from exc


# ensembl-variation Utils/Constants.pm %OVERLAP_CONSEQUENCES at release/116, SO term to impact.
TERM_IMPACT = {
    **dict.fromkeys((
        "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant",
        "stop_lost", "start_lost", "transcript_amplification", "feature_elongation", "feature_truncation",
    ), "HIGH"),
    **dict.fromkeys(("inframe_insertion", "inframe_deletion", "missense_variant", "protein_altering_variant", "TFBS_ablation"), "MODERATE"),
    **dict.fromkeys((
        "splice_donor_5th_base_variant", "splice_region_variant", "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
        "stop_retained_variant", "synonymous_variant",
    ), "LOW"),
    **dict.fromkeys((
        "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant",
        "coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_amplification",
        "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
        "regulatory_region_variant", "intergenic_variant", "sequence_variant",
    ), "MODIFIER"),
}
# The plugin fields written only on a row whose terms let the plugin run. Downstream.pm:97-99
# asks for ``frameshift_variant`` and no term holding ``splice``; NMD.pm:65 and :83 and LOFTEE's
# LoF.pm:166-167 for one of ``stop_gained``, ``frameshift_variant`` and the two splice site terms.
DOWNSTREAM_FIELDS = ("DownstreamProtein", "ProteinLengthChange")
TRUNCATING_TERM_FIELDS = ("NMD", "LoF", "LoF_info", "LoF_filter", "LoF_flags")
TRUNCATING_TERMS = frozenset({STOP_GAINED_TERM, FRAMESHIFT_TERM, "splice_donor_variant", SPLICE_ACCEPTOR_TERM})


def terms_impact(terms: Collection[str]) -> str | None:
    """The IMPACT a row of these terms carries, its most severe term's (OutputFactory.pm:1337
    and :1344 read the first consequence by rank), or ``None`` for a term without a rank."""
    if not terms or any(t not in CONSEQUENCE_RANK for t in terms):
        return None
    return TERM_IMPACT[min(terms, key=CONSEQUENCE_RANK.__getitem__)]


def terms_write_field(name: str, terms: Collection[str]) -> bool | None:
    """Whether a row of these terms permits the plugin that writes the field ``name`` to run, or
    ``None`` for a field no term decides."""
    if field_named(name, DOWNSTREAM_FIELDS):
        return FRAMESHIFT_TERM in terms and not any("splice" in t for t in terms)
    if field_named(name, TRUNCATING_TERM_FIELDS):
        return bool(TRUNCATING_TERMS & set(terms))
    return None


def hgvsp_frameshift_form(value: str) -> bool:
    """Whether an HGVSp is a frameshift's notation: ``fs``, or ``ext`` where the frameshift
    begins at the stop codon."""
    change = value.partition(":p.")[2]
    return "fs" in change or "ext" in change


# A row whose two Consequence sets differ by the terms of a stop codon class or of
# ``splice_term_past_the_record_span`` differs in the cells
# that are read from the terms as well, and one rule reads them all: a cell follows the terms
# where its two values are the ones the two sets give. IMPACT is the impact of each side's
# most severe term. A field of Downstream, NMD or LOFTEE is written on the side whose terms
# let the plugin run and left out on the side whose terms do not; the value itself is the
# plugin's. An HGVSp is no frameshift's notation on the side that holds
# ``stop_retained_variant``, whose predicate keeps ``frameshift`` from holding
# (Utils/VariationEffect.pm:1549, which the notation's type reads at
# TranscriptVariationAllele.pm:2053), and a frameshift's on the side without it, both naming
# one protein. A cell that differs in any other way stays a
# difference: an IMPACT that is not its side's terms', a plugin field on a side whose terms
# keep the plugin out or on both sides with two values, a notation of another form.
# Every field no term decides stays a difference.
def term_pair_following_cell(e: Entry, a: Entry, name: str, ev: str, av: str) -> bool:
    """Whether the differing values ``ev`` and ``av`` of the field ``name`` on the reference
    entry ``e`` and the candidate entry ``a`` of one key are the ones the two entries'
    Consequence sets give: the IMPACT of each set, a Downstream, NMD or LOFTEE field present
    exactly where the set permits the plugin to run, or an HGVSp that is a frameshift's notation
    exactly on the side without ``stop_retained_variant``."""
    expected_terms = set(consequence_set(e.consequence).split(","))
    actual_terms = set(consequence_set(a.consequence).split(","))
    if expected_terms == actual_terms or ev == av:
        return False
    if field_named(name, IMPACT_FIELDS):
        return ev == terms_impact(expected_terms) and av == terms_impact(actual_terms)
    expected_written, actual_written = terms_write_field(name, expected_terms), terms_write_field(name, actual_terms)
    if expected_written is not None:
        return expected_written != actual_written and blank(ev) != expected_written and blank(av) != actual_written
    if field_named(name, HGVSP_FIELDS):
        expected_kept, actual_kept = STOP_RETAINED_TERM in expected_terms, STOP_RETAINED_TERM in actual_terms
        if expected_kept == actual_kept or blank(ev) or blank(av) or ev.partition(":p.")[0] != av.partition(":p.")[0]:
            return False
        return hgvsp_frameshift_form(ev) != expected_kept and hgvsp_frameshift_form(av) != actual_kept
    return False


# The built-in row classes that accept a transcript row's Consequence, each with the predicate
# deciding whether a reference entry and the candidate's of one key have its shape; a row is
# the first class's whose predicate holds, and the JSON record's ``most_severe_consequence``
# is left uncompared on a record carrying one.
TERM_PAIR_CLASSES = {
    REFSEQ_EDITED_MODEL_DELETED_SPAN_TERM: deleted_span_term_pair,
    REFSEQ_EDITED_MODEL_TAIL_TERM: tail_term_pair,
    MATURE_MIRNA_UNSHIFTED_SPAN: mature_mirna_unshifted_span_pair,
    STOP_CODON_UNSHIFTED_SPAN: stop_codon_unshifted_span_pair,
    STOP_RETAINED_WITHOUT_STOP_CODON_TEST: stop_retained_without_stop_codon_test_pair,
    STOP_CODON_OVERLAP_GENOMIC_SPAN: stop_codon_overlap_genomic_span_pair,
    STOP_RETAINED_THREE_BASE_STOP_CODON_TEST: stop_retained_three_base_stop_codon_test_pair,
    STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST: stop_codon_removed_without_stop_codon_test_pair,
    STOP_CODON_OVERLAP_INSERTION_RANGE: stop_codon_overlap_insertion_range_pair,
    SPLICE_TERM_PAST_THE_RECORD_SPAN: splice_term_past_the_record_span_pair,
}
# The row classes of ``TERM_PAIR_CLASSES`` whose rows differ in the cells read from the two
# sets of terms as well: on a row such a class takes, ``term_pair_following_cell`` reads each
# differing cell of ``TERM_PAIR_FOLLOWING_FIELDS``, and a cell it accepts is counted under the
# row's class, the row once. Such a class is a field class of those fields too, so a caller
# that applies a class by its name (``field_class_accepts``) reads the same rule;
# ``compare_entries`` asks a row's own class before every other field class and asks it of
# no other row.
TERM_PAIR_FOLLOWING_CLASSES: tuple[str, ...] = (
    STOP_RETAINED_WITHOUT_STOP_CODON_TEST,
    STOP_CODON_OVERLAP_GENOMIC_SPAN,
    STOP_RETAINED_THREE_BASE_STOP_CODON_TEST,
    STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST,
    STOP_CODON_OVERLAP_INSERTION_RANGE,
    SPLICE_TERM_PAST_THE_RECORD_SPAN,
)
TERM_PAIR_FOLLOWING_FIELDS = (*IMPACT_FIELDS, *DOWNSTREAM_FIELDS, *TRUNCATING_TERM_FIELDS, *HGVSP_FIELDS)
FIELD_CLASSES.update(dict.fromkeys(TERM_PAIR_FOLLOWING_CLASSES, TERM_PAIR_FOLLOWING_FIELDS))


# ``plugin_value_follows_splice_term_past_the_record_span``: a field of Downstream, NMD or
# LOFTEE that one engine alone writes on a row ``splice_term_past_the_record_span_pair``
# takes. The three plugins run on a row's terms. LOFTEE sets its splice gate on
# ``splice_acceptor_variant`` or ``splice_donor_variant`` (LoF.pm:167) and returns without a
# confidence for a row that passes neither that gate, the one on ``stop_gained`` and
# ``frameshift_variant`` nor its own splice prediction (:222-232); NMD runs on one of the
# same four terms (NMD.pm:65 and :83). Where the donor or acceptor term the reference alone
# holds is its row's only such term, the fields of both stand on the reference's side alone.
# Downstream writes no value beside a term that holds ``splice`` (Downstream.pm:97-99):
# where the reference holds ``frameshift_variant`` beside the term and the candidate's row
# holds no splice term, its fields stand on the candidate's side alone. ``intron_variant`` is
# a term of none of those tests, so the intron term of the class moves no field. The shape:
# the row's two entries are a pair of that predicate, read with the record's alleles where
# an entry carries them (``Entry.reference``) and without them where none does, and the field
# is a cell ``term_pair_following_cell`` reads from the two term sets. A row the predicate
# refuses, a field on a side whose terms keep the plugin out, a field both sides carry and
# every field no term decides stay differences.
PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN = "plugin_value_follows_splice_term_past_the_record_span"
FIELD_CLASSES[PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN] = DOWNSTREAM_FIELDS + TRUNCATING_TERM_FIELDS
# The classes whose predicate reads the reference bases an entry's allele stands in for
# (``Entry.reference``): under one of them a default or tab comparison reads the records of
# ``--input`` for its rows under an identifier (``read_entry_references``).
RECORD_REFERENCE_CLASSES = (SPLICE_TERM_PAST_THE_RECORD_SPAN, PLUGIN_VALUE_FOLLOWS_SPLICE_TERM_PAST_THE_RECORD_SPAN)


# ``plugin_value_follows_stop_codon_term_pair``: a field of Downstream, NMD or LOFTEE that one
# engine alone writes on a row whose two Consequence sets differ as a stop codon class admits.
# The three plugins run on a row's terms (Downstream.pm:97-99, NMD.pm:65 and :83,
# LoF.pm:166-167), so where the terms differ the field stands on the side whose terms let the
# plugin run and is left out on the other. The shape: the row's two entries are a pair of one
# predicate of ``STOP_CODON_PLUGIN_VALUE_PAIRS``, and the field is a cell
# ``term_pair_following_cell`` reads from the two term sets. A row of any other pair, a field
# on a side whose terms keep the plugin out, a field both sides carry and every field no term
# decides stay differences.
PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR = "plugin_value_follows_stop_codon_term_pair"
STOP_CODON_PLUGIN_VALUE_PAIRS = (
    stop_retained_without_stop_codon_test_pair,
    stop_codon_overlap_insertion_range_pair,
    stop_retained_three_base_stop_codon_test_pair,
)
FIELD_CLASSES[PLUGIN_VALUE_FOLLOWS_STOP_CODON_TERM_PAIR] = DOWNSTREAM_FIELDS + TRUNCATING_TERM_FIELDS


def freqs_items(value: str) -> dict[str, str] | None:
    """The ``population:allele`` keys of a ``FREQS`` value with their frequencies, or ``None``
    when an item is not a ``population:allele:value`` triple; the items are joined with ``,``
    in the default and tab formats and ``&`` in the VCF (``OutputFactory.pm:1152-1159``)."""
    items: dict[str, str] = {}
    for item in re.split(r"[,&]", value):
        population, _, rest = item.partition(":")
        allele, sep, frequency = rest.partition(":")
        if not sep:
            return None
        items[f"{population}:{allele}"] = frequency
    return items


def freqs_differ_by_zero_read_as_na(expected_value: str, actual_value: str) -> bool:
    """Whether two ``FREQS`` values carry the same items and differ only where the reference
    reads ``NA`` and the candidate ``0`` (``frequency_filter_zero_value``)."""
    if blank(expected_value) or blank(actual_value):
        return False
    expected_items, actual_items = freqs_items(expected_value), freqs_items(actual_value)
    if expected_items is None or actual_items is None or expected_items.keys() != actual_items.keys():
        return False
    changed = [k for k in expected_items if expected_items[k] != actual_items[k]]
    return bool(changed) and all(expected_items[k] == "NA" and as_number(actual_items[k]) == 0 for k in changed)


def freqs_all_zero(e: Entry) -> bool:
    """Whether the entry's ``FREQS`` carries at least one item and every item's frequency is
    0: the rows of a record ``--freq_filter include`` keeps for a zero frequency."""
    value = next((e.fields[f] for f in FREQS_FIELDS if f in e.fields), ABSENT)
    if blank(value):
        return False
    items = freqs_items(value)
    return bool(items) and all(as_number(v) == 0 for v in items.values())


def registry_class(expected_consequence: str, actual_consequence: str, allele: str, classes: Collection[str]) -> str | None:
    """The registry class (``REGISTRY_CLASSES``, restricted to ``classes``) whose excluding
    rule accepts the pair of consequence sets, or ``None``: the SNP/indel comparator's
    ``EXCLUDING_RULES`` decide, so the pairs set aside here are the pairs its adjusted figure
    sets aside."""
    if not classes:
        return None
    script_dir = str(Path(__file__).resolve().parent)
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
    import compare_vep_outputs as cvo

    for rule in cvo.EXCLUDING_RULES:
        if rule.taxonomy_class in classes and rule.matches(expected_consequence, actual_consequence, allele):
            return rule.taxonomy_class
    return None


def structural_variant_comparator():
    """The structural-variant comparator (``scripts/validation/compare_sv_concordance.py``),
    imported from beside this script's directory: its predicates decide the classes of
    ``STRUCTURAL_VARIANT_CLASSES`` here, so the rows set aside are the rows its adjusted
    figures set aside."""
    validation_dir = str(Path(__file__).resolve().parent.parent / "validation")
    if validation_dir not in sys.path:
        sys.path.insert(0, validation_dir)
    import compare_sv_concordance as csc

    return csc


# The structural-variant comparator's tuple: Location, Allele, Feature, Feature_type and the
# Consequence set.
StructuralTuple = tuple[str, str, str, str, str]


class StructuralVariants:
    """What the built-in structural-variant classes read beside the entries: the input records
    (the comparator's ``VcfVariant``, whose ALT and INFO the skip routes read), the candidate
    cache's three views (the transcript ids per chromosome, the spans and the models), and each
    consequence entry's tuple as the structural-variant comparator keys a row. A text entry's
    Location is its own; a VCF or JSON entry's is the one vep-rs prints for its input record
    under the entry's allele (``vep_rs_record_rows``), so an entry of an output line no input
    record matches has no tuple and is outside every class. ``classes`` is the subset of
    ``STRUCTURAL_VARIANT_CLASSES`` requested; the skipped-record class needs the input VCF, the
    two breakend classes the candidate JSON cache."""

    def __init__(self, classes: Collection[str], input_path: Path | None, cache: Path | None) -> None:
        self.classes = {k for k in STRUCTURAL_VARIANT_CLASSES if k in classes}
        if REFERENCE_SKIPPED_RECORD in self.classes and input_path is None:
            raise ComparatorError(f"--masks {REFERENCE_SKIPPED_RECORD} needs --input, the input VCF")
        cache_reading = [k for k in CACHE_READING_CLASSES if k in self.classes]
        if cache_reading and cache is None:
            raise ComparatorError(f"--masks {','.join(cache_reading)} needs --json-cache, the candidate JSON cache")
        self.csc = structural_variant_comparator()
        self.input_variants = []
        if input_path is not None:
            try:
                self.input_variants = self.csc.parse_vcf(input_path)
            except OSError as exc:
                raise ComparatorError(f"cannot read {input_path}: {exc.strerror or exc}") from exc
        self.by_chr: dict[str, set[str]] = {}
        self.spans: dict = {}
        self.models: dict = {}
        if cache is not None:
            if not (cache / "transcripts").is_dir():
                raise ComparatorError(f"--json-cache: {cache} has no transcripts/ directory")
            self.by_chr, self.spans, self.models = self.csc.load_cache_authority(cache)
        self._record_locations: dict[int, tuple[dict[str, str], str | None]] = {}

    def location(self, e: Entry) -> str | None:
        """The entry's Location as the structural-variant comparator spells it: the text
        formats' own, normalised; a VCF or JSON entry's from its input record, the Location of
        the record's row under the entry's allele, else the record's one Location, else none."""
        if e.location is not None:
            return self.csc.normalize_location(e.location)
        if e.record is None or e.record >= len(self.input_variants):
            return None
        if e.record not in self._record_locations:
            rows = sorted(self.csc.vep_rs_record_rows(self.input_variants[e.record]))
            by_allele = {allele: location for location, allele in rows}
            locations = {location for location, _ in rows}
            self._record_locations[e.record] = (by_allele, locations.pop() if len(locations) == 1 else None)
        by_allele, only = self._record_locations[e.record]
        return by_allele.get(self.csc.normalize_allele(dash(e.allele)), only)

    def tuples(self, entries: list[Entry]) -> dict[str, StructuralTuple]:
        """Each consequence entry's tuple by its key; a record entry has none."""
        out: dict[str, StructuralTuple] = {}
        for e in entries:
            if e.key.endswith("|record"):
                continue
            location = self.location(e)
            if location is None:
                continue
            feature_type = e.key.rsplit("|", 1)[-1].partition("#")[0]
            allele = self.csc.normalize_allele(dash(e.allele))
            out[e.key] = (location, allele, dash(e.feature), feature_type, consequence_set(e.consequence))
        return out


def list_items(value: str, separator: str) -> list[str] | None:
    """The items of a list-valued field: the ``separator``-joined parts of a text format's
    value, or the elements of the JSON array ``json_scalar`` joined (``JSON_ARRAY``), each
    as its compact text; under a ``PAIR_GROUPS`` separator the parts are groups of ``:``-joined
    pairs, each group read with its pairs sorted; ``None`` when the value is not such a list."""
    if separator.startswith(PAIR_GROUPS):
        groups = value.split(separator.removeprefix(PAIR_GROUPS))
        return [":".join(sorted(group.split(":"))) for group in groups]
    if separator != JSON_ARRAY:
        return value.split(separator)
    try:
        items = json.loads(f"[{value}]")
    except ValueError:
        return None
    return [json_compact(item) for item in items]


def same_items(expected: str, actual: str, separator: str) -> bool:
    """Whether two list-valued field texts carry the same items the same number of times,
    in any order; an empty or absent value agrees with itself alone."""
    if expected in ("", "-", ABSENT) or actual in ("", "-", ABSENT):
        return expected == actual
    expected_items = list_items(expected, separator)
    actual_items = list_items(actual, separator)
    if expected_items is None or actual_items is None:
        return False
    return Counter(expected_items) == Counter(actual_items)


def list_parts(expected: str, actual: str) -> tuple[list[str], list[str]]:
    """The items of two list-valued field texts of one format, split at the format's list
    separator: ``&`` where the reference's text carries one, as a VCF value does (VCF.pm:391),
    ``,`` otherwise; no items for a text carrying both, which neither separator splits into
    the items its engine wrote."""
    if any("," in text and "&" in text for text in (expected, actual)):
        return [], []
    separator = "&" if "&" in expected else ","
    return expected.split(separator), actual.split(separator)


def reordered_items(expected: str, actual: str) -> bool:
    """Whether two list-valued field texts carry the same items the same number of times in
    different orders, so two or more of them (``nearest_exon_jb_item_order``)."""
    expected_items, actual_items = list_parts(expected, actual)
    return expected_items != actual_items and Counter(expected_items) == Counter(actual_items)


def reordered_pair_groups(expected: str, actual: str) -> bool:
    """Whether two ``5UTR_annotation`` texts carry the same groups of the same ``key=value``
    pairs with the groups, or the pairs inside one, in different orders, the candidate's pairs
    in byte order of their keys inside every group (``utrannotator_annotation_pair_order``)."""
    expected_groups, actual_groups = ([group.split(":") for group in groups] for groups in list_parts(expected, actual))
    if expected_groups == actual_groups:
        return False
    if not all("=" in pair for group in (*expected_groups, *actual_groups) for pair in group):
        return False
    if any(group != sorted(group, key=lambda pair: pair.partition("=")[0]) for group in actual_groups):
        return False
    return Counter(tuple(sorted(group)) for group in expected_groups) == Counter(tuple(sorted(group)) for group in actual_groups)


def order_only_fields(fmt: str) -> dict[str, str]:
    """The fields the ``order_only`` class compares as multisets in ``fmt``, each with its
    list separator: ``,`` in the default and tab formats (``convert_arrayref``,
    Utils.pm:185-198), ``&`` in VCF (VCF.pm:391), and in JSON the arrays of
    ``ORDER_ONLY_JSON_FIELDS``."""
    if fmt == "json":
        return dict(ORDER_ONLY_JSON_FIELDS)
    separator = "&" if fmt == "vcf" else ","
    fields = dict.fromkeys(ORDER_ONLY_COLUMNS, separator)
    fields[ORDER_ONLY_PAIR_GROUP_COLUMN] = PAIR_GROUPS + separator
    return fields


POSITION_RE = re.compile(r"(\w+?)_position$", re.IGNORECASE)
POSITION_PART_RE = re.compile(r"^\d+(/\d)?\d*$")
TOOL_SCORE_RE = re.compile(r"([a-z_]+)?\(?([\d.]+)?\)?", re.IGNORECASE)
DOMAIN_RE = re.compile(r"(\w+):(\w+)")
PERL_NUMBER_RE = re.compile(r"^\s*[+-]?(\d+\.?\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?)")


def json_floor_fields(column: str) -> tuple[str, ...] | None:
    """The JSON fields a floor on the default-format ``column`` is counted under
    (JSON.pm:311-385): the name lower-cased and renamed as ``%RENAME_KEYS`` renames it,
    ``Feature`` under the per-kind id, SIFT and PolyPhen under either of their two keys, a
    position under either of its start and end keys (JSON.pm:340-341 delete the start when
    it is ``?`` and keep the end, which the default format prints as ``?-N``); ``None`` for
    a column with no row count in the JSON format: a record-level column, ``Feature_type``,
    which chooses the consequence list an entry joins and is deleted from the entry
    (JSON.pm:348-350), and ``Extra``, the text formats' container of the keyed fields
    (VEP_output.pm:114-133)."""
    key = column.lower()
    if column in RECORD_LEVEL_COLUMNS or RECORD_LEVEL_COLUMN_RE.search(column) or key in ("feature_type", "extra"):
        return None
    if key == "feature":
        return tuple(f"{kind}_id" for _, kind in JSON_CONSEQUENCE_LISTS if kind != "intergenic")
    pm = POSITION_RE.match(column)
    if pm:
        coord = pm.group(1).lower()
        return (f"{coord}_start", f"{coord}_end")
    if key in ("sift", "polyphen"):
        return (f"{key}_prediction", f"{key}_score")
    return (JSON_RENAME.get(key, key),)


def field_objects(text: str) -> list[dict]:
    """The objects of a field whose text is one JSON object or an array of them as
    ``json_scalar`` renders it, the items' compact texts joined with ``,`` (a plugin writing one
    object per association, ``opentargets``, reads as ``{..},{..}``); empty for a scalar or any
    other text."""
    if not text.startswith("{") or not text.endswith("}"):
        return []
    try:
        items = json.loads(f"[{text}]")
    except ValueError:
        return []
    return [item for item in items if isinstance(item, dict)]


def json_object_member(e: Entry, parent: str, key: str) -> str | None:
    """The value of ``key`` inside the entry's ``parent`` field when that field is a JSON
    object, or an array of objects, carrying the key in any case (a plugin that nests its fields
    under its own name, ``alphamissense`` holding ``am_class``; the JSON factory lower-cases the
    object's key and leaves the members as the plugin spelt them): the first object's value,
    else ``None``."""
    wanted = key.lower()
    for members in field_objects(e.fields.get(parent, "")):
        for name, value in members.items():
            if name.lower() == wanted:
                return json_scalar(value)
    return None


def populated_under(entries: list[Entry], fields: Collection[str]) -> int:
    """The consequence entries (never a record's own entry, which repeats ``strand`` and
    the coordinates) carrying any of ``fields`` with a value other than ``-``: the JSON
    factory drops a key whose default-format value is ``-`` and keeps the allele keys
    (JSON.pm:327), so a present key is a populated column, an empty ``domains`` array
    included (JSON.pm:363-371 keep the key when no item re-splits). An item ``parent.key``
    counts an entry whose ``parent`` object carries ``key`` (``json_object_member``)."""

    def populated(e: Entry, f: str) -> bool:
        if f in e.fields:
            return e.fields[f] != "-"
        parent, dot, key = f.partition(".")
        if not dot:
            return False
        value = json_object_member(e, parent, key)
        return value is not None and value != "-"

    return sum(1 for e in entries if not e.key.endswith("|record") and any(populated(e, f) for f in fields))


def json_floor_location(column: str, expected: list[Entry], actual: list[Entry]) -> tuple[str, ...] | str | None:
    """Where a JSON floor on ``column`` is counted, read against the parsed entries of both
    sides: the fields of ``json_floor_fields`` when any entry carries one of them (``None``
    for a column the format has no row count for); else, for a column outside
    ``CORE_COLUMNS``, a plugin object's member, ``<object>.<key>``, when a consequence entry of
    either side carries an object (or an array of objects) whose member is the column's key;
    else the custom-source note when a record's ``custom_annotations`` of either side names a
    source the column spells (``<source>`` or ``<source>_<field>``); else the fields of
    ``json_floor_fields``, which no entry carries, so the floor reads as unmet on both sides."""
    fields = json_floor_fields(column)
    if fields is None:
        return None
    entries = [*expected, *actual]
    if any(f in e.fields for e in entries for f in fields):
        return fields
    key = column.lower()
    if key not in CORE_COLUMNS:
        # A field whose ``{``-opening text parsed as no object on one entry is not parsed again
        # on the entries after it; an object field is parsed on every entry, since a plugin may
        # omit the member on some.
        not_objects: set[str] = set()
        for e in entries:
            if e.key.endswith("|record"):
                continue
            for parent, value in e.fields.items():
                if parent in not_objects or not value.startswith("{"):
                    continue
                objects = field_objects(value)
                if not objects:
                    not_objects.add(parent)
                elif any(name.lower() == key for members in objects for name in members):
                    return (f"{parent}.{key}",)
    sources: set[str] = set()
    for e in entries:
        if e.key.endswith("|record"):
            for annotations_by_source in field_objects(e.fields.get(CUSTOM_ANNOTATIONS_FIELD, "")):
                sources.update(annotations_by_source.keys())
    if any(column == source or column.startswith(f"{source}_") for source in sources):
        return FLOOR_CUSTOM_SOURCE_NOT_APPLICABLE
    return fields


def declared_columns(fmt: str, path: Path) -> set[str]:
    """The columns a tab or VCF output's header declares, read from the file's leading ``#``
    lines alone: the tab ``#`` header line's columns, the VCF CSQ ``Format:`` list
    (``csq_format_fields``) of the last ``##INFO=<ID=CSQ`` line, the one ``parse_vcf`` reads
    the entries by (``--keep_csq`` writes the kept list first); empty for a file without one.
    Every row of these formats carries every declared column, so a column no header declares
    appears in no row."""
    opener = gzip.open if str(path).endswith(".gz") else open
    declared: set[str] = set()
    try:
        with opener(path, "rt", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.rstrip("\r\n")
                if not line.startswith("#"):
                    break
                if fmt == "vcf" and line.startswith("##INFO=<ID=CSQ"):
                    declared = set(csq_format_fields(line))
                elif fmt == "tab" and not line.startswith("##"):
                    return set(line[1:].split("\t"))
    except (OSError, EOFError, UnicodeDecodeError, zlib.error) as exc:
        raise ComparatorError(f"cannot read {path}: {getattr(exc, 'strerror', None) or exc}") from exc
    return declared


# ``colocated_list_cut_at_cache_cap``: the Ensembl cache dump ends a stored ``pubmed`` or
# ``var_synonyms`` value of a known variant at 1,024 bytes, inside an id as readily as after a
# separator, and Ensembl VEP prints the stored value as it stands (OutputFactory.pm:1090,
# :1126 and :1141, OutputFactory/JSON.pm:452-457 and :469 at release/116.2), the cut item as
# an id. The candidate reads a value of exactly that length ending in none of its separators
# (``,`` in ``pubmed``; ``--``, ``::`` and ``,`` in ``var_synonyms``) as the bytes before its
# last separator, a ``var_synonyms`` group whose only id is the cut one going with its source
# name and the ``--`` before that, and writes no field for a value left with no id. The shape:
# the reference's list is the candidate's plus one last item, or plus one last group of one
# id, and the reference's own text shows the stored value at exactly 1,024 bytes and ending in
# no separator. A shorter list, a list ending in a separator, a pair whose earlier items
# differ and a candidate that keeps the source name of such a group, bare in a row or with an
# empty array in JSON, stay differences.
#
# The default, tab and VCF formats join the lists of a row's known variants, ``PUBMED`` ids
# with ``,`` and ``VAR_SYNONYMS`` stored values with ``--`` (VCF.pm:401 writes each ``,`` as
# ``&``), so a stored value is a run of the row's value. The candidate's value is the
# reference's without the last item of one or more runs, each exactly 1,024 bytes from the
# value's start or the end of a joining separator to the value's end or the start of one; a
# run left with no id goes whole with one joining separator, and the field with it when the
# row joins nothing else. A run of that length spanning the values of two known variants
# reads as one stored value: the row's text does not tell the two apart.
#
# The JSON format holds one ``colocated_variants`` entry per known variant. A ``pubmed`` array
# holds the stored value's items, so its items joined with ``,`` are the stored value. A
# ``var_synonyms`` object is rebuilt as its groups, ``source::id,id``, joined with ``--``; the
# length does not depend on the groups' order, which the object loses, and the cut item's
# group is the one whose last id the candidate lacks, or the source the candidate lacks, a
# name with one id or with no ``::``. The object does not show which group the stored value
# ends in: a candidate lacking the last id of any one group, or any one group of one id, is
# read as the shape on the rebuilt length alone, where a row's text holds it to the last
# group of a run. Three shapes do not show the stored length there and stay
# differences: a second group listing no id (its stored text is the name with or without
# ``::``), a source name the stored value repeats (JSON.pm:457 keeps the ids of its last
# group alone) and an id whose JSON number is shorter than its stored text (``numberify``,
# Utils.pm:241, prints ``999001.0010`` as ``999001.001``).
#
# An accepted field is counted under the class once per entry: every row of the record in the
# default, tab and VCF formats, the record's object in JSON. The class is no member of
# ``FIELD_CLASSES`` or of any other class list of this module: ``compare_entries`` applies it
# when ``Documented.masks`` names it.
COLOCATED_LIST_CUT_AT_CACHE_CAP = "colocated_list_cut_at_cache_cap"
CACHE_LIST_CAP = 1024
# The separators of each list column's stored value, the first of them the one that joins the
# values of a row's known variants.
CUT_LIST_SEPARATORS = {"PUBMED": (b",",), "VAR_SYNONYMS": (b"--", b"::", b",")}
COLOCATED_VARIANTS_FIELD = "colocated_variants"


def cut_list_whole_items(stored: bytes, column: str) -> bytes | None:
    """The whole items of a stored ``column`` value the cache cut inside its last item: the
    bytes before its last separator, or before the ``--`` of the last group when that
    separator is the group's ``::`` (a ``var_synonyms`` source whose only id is cut goes with
    it), none of them when they hold no id (a ``pubmed`` value without a ``,``, a
    ``var_synonyms`` value without an id after a ``::``). ``None`` for a value not cut so:
    one of another length than ``CACHE_LIST_CAP`` bytes, or one ending in a separator."""
    separators = CUT_LIST_SEPARATORS[column]
    if len(stored) != CACHE_LIST_CAP or stored.endswith(separators):
        return None
    whole = stored[: max(0, *(stored.rfind(separator) for separator in separators))]
    if column == "VAR_SYNONYMS":
        if stored.startswith(b"::", len(whole)):
            whole = whole.rpartition(b"--")[0]
        listed = (group.partition(b"::")[2] for group in whole.split(b"--"))
        if not any(ids.strip(b",") for ids in listed):
            return b""
    return whole


def row_list_bytes(value: str, vcf: bool) -> bytes:
    """A default, tab or VCF list value as the bytes of the stored values it joins: none for
    an empty or absent value, ``;`` for the ``%3B`` the default and VCF formats write
    (VEP_output.pm:118, VCF.pm:402) and, in the VCF format, ``,`` for ``&`` (VCF.pm:401)."""
    text = "" if blank(value) else value.replace("%3B", ";")
    return (text.replace("&", ",") if vcf else text).encode()


def row_cut_lists_pair(
    reference: bytes,
    candidate: bytes,
    column: str,
    start: int = 0,
) -> bool:
    """Whether ``candidate`` is ``reference``, a row's ``column`` value (``row_list_bytes``),
    without the cut last item of one or more of the stored values it joins
    (``colocated_list_cut_at_cache_cap``). A run of ``CACHE_LIST_CAP`` bytes that begins at the
    value's start or after a joining separator, at ``start`` or later, and ends at the
    value's end or before a joining separator is read as ``cut_list_whole_items`` reads it,
    a run left with no id going with one joining separator; the bytes up to the run's last
    whole item are the candidate's, and the rest is the candidate's as it stands or with a
    later run read the same way."""
    join = CUT_LIST_SEPARATORS[column][0]
    size = len(reference)
    joined = (at for at in range(size) if reference.startswith(join, at))
    for begin in (0, *(at + len(join) for at in joined)):
        end = begin + CACHE_LIST_CAP
        if begin < start or not (end == size or reference.startswith(join, end)):
            continue
        whole = cut_list_whole_items(reference[begin:end], column)
        if whole is None:
            continue
        if whole:
            read, rest = reference[:begin] + whole, reference[end:]
        elif end < size:
            read, rest = reference[:begin], reference[end + len(join) :]
        else:
            read, rest = reference[: max(0, begin - len(join))], b""
        if not candidate.startswith(read):
            continue
        value = read + rest
        later = len(read) + (len(join) if whole else 0)
        if value == candidate or row_cut_lists_pair(value, candidate, column, later):
            return True
    return False


def json_cut_pubmed(ids: object, kept: object) -> bool:
    """Whether a JSON ``pubmed`` array ``ids`` of the reference and the candidate's ``kept``
    (``None`` for no member) have the shape of ``colocated_list_cut_at_cache_cap``: the candidate
    lists the reference's ids but the last, and has no member when none is left; the
    reference's ids joined with ``,``, the stored value, are cut inside the last
    (``cut_list_whole_items``)."""
    if not isinstance(ids, list) or not ids:
        return False
    left = ids[:-1] or None
    if (left is None) != (kept is None) or json_compact(left) != json_compact(kept):
        return False
    return cut_list_whole_items(json_scalar(ids).encode(), "PUBMED") is not None


def json_cut_synonyms(groups: object, kept: object) -> bool:
    """Whether a JSON ``var_synonyms`` object ``groups`` of the reference and the
    candidate's ``kept`` (``None`` for no member) have the shape of
    ``colocated_list_cut_at_cache_cap``. One source of the reference is the stored value's last
    group: the candidate lists its ids but the last, or lacks the source when the reference
    lists one id under it or none (a name with no ``::``), and has no member when no id is
    left under any source. Every other source lists an id, so the groups, ``source::id,id``,
    joined with ``--`` are the stored value, cut inside its last item
    (``cut_list_whole_items``). The object's keys are not in the stored order, so any source
    the candidate's object fits this way is taken as the last."""
    if not isinstance(groups, dict):
        return False
    if not all(isinstance(listed, list) for listed in groups.values()):
        return False
    for last, ids in groups.items():
        others = {source: listed for source, listed in groups.items() if source != last}
        if not all(others.values()):
            continue
        left = {**others, last: ids[:-1]} if len(ids) > 1 else others
        if not any(left.values()):
            left = None
        if (left is None) != (kept is None) or json_compact(left) != json_compact(kept):
            continue
        texts = [f"{name}::{json_scalar(listed)}" for name, listed in others.items()]
        texts.append(f"{last}::{json_scalar(ids)}" if ids else last)
        if cut_list_whole_items("--".join(texts).encode(), "VAR_SYNONYMS") is not None:
            return True
    return False


# The members of a JSON ``colocated_variants`` entry that hold a stored list, each with the
# predicate of its cut shape.
JSON_CUT_LISTS = {"pubmed": json_cut_pubmed, "var_synonyms": json_cut_synonyms}


def json_cut_lists_pair(expected: str, actual: str) -> bool:
    """Whether two JSON ``colocated_variants`` values have the shape of
    ``colocated_list_cut_at_cache_cap``: the same known variants in the same order, every member
    of each equal but the ``pubmed`` (``json_cut_pubmed``) or the ``var_synonyms``
    (``json_cut_synonyms``) of at least one."""
    expected_variants, actual_variants = field_objects(expected), field_objects(actual)
    if len(expected_variants) != len(actual_variants):
        return False
    cut = False
    for reference, candidate in zip(expected_variants, actual_variants):
        for key in reference.keys() | candidate.keys():
            ours, theirs = reference.get(key), candidate.get(key)
            in_both = key in reference and key in candidate
            if in_both and json_compact(ours) == json_compact(theirs):
                continue
            cut_shape = JSON_CUT_LISTS.get(key)
            if cut_shape is None or not cut_shape(ours, theirs):
                return False
            cut = True
    return cut


def cut_list_reading(e: Entry, name: str, ev: str, av: str) -> str | None:
    """The reference's differing field ``name`` of the entry ``e`` as the candidate reads the
    cut lists in it (``colocated_list_cut_at_cache_cap``), ``None`` for a field holding none:
    the candidate's value ``av`` for a ``PUBMED`` or ``VAR_SYNONYMS`` value
    (``row_cut_lists_pair``; an entry without a Location is a VCF format's) and for a JSON
    ``colocated_variants`` value (``json_cut_lists_pair``) of the class's shape; for the
    default format's ``Extra`` text, the reference's ``ev`` with each such pair carrying the
    candidate's value, or left out where the candidate writes none, the text's other pairs
    as they stand."""
    if name == COLOCATED_VARIANTS_FIELD:
        return av if json_cut_lists_pair(ev, av) else None
    vcf = e.location is None

    def cut_pair(column: str, reference: str, candidate: str) -> bool:
        if column not in CUT_LIST_SEPARATORS:
            return False
        values = row_list_bytes(reference, vcf), row_list_bytes(candidate, vcf)
        return row_cut_lists_pair(*values, column)

    if name != "Extra":
        return av if cut_pair(name, ev, av) else None
    written = extra_pairs(av)
    pairs: list[str] = []
    read = False
    for pair in ev.split(";"):
        key, _, value = pair.partition("=")
        if cut_pair(key, value, written.get(key, "-")):
            read = True
            if key not in written:
                continue
            pair = f"{key}={written[key]}"
        pairs.append(pair)
    return (";".join(pairs) or "-") if read else None


def cut_list_example(e: Entry, name: str, ev: str, av: str) -> str:
    """The ``mask_examples`` line of a field ``colocated_list_cut_at_cache_cap`` accepts: the two
    values' lengths, since a value at the cap runs to a thousand characters of ids."""
    written = 0 if blank(av) else len(av)
    return f"{e.key}: {name} expected {len(ev)} characters, got {written}"


def compare_entries(
    expected: list[Entry],
    actual: list[Entry],
    documented: Documented,
    floors: dict[str, int] | None = None,
    field_pairs: dict[str, set[tuple[str, str]]] | None = None,
    check_order: bool = True,
    examples: int = EXAMPLES_PER_GROUP,
    nearest_sites: NearestSites | None = None,
    order_only: dict[str, str] | None = None,
    floor_fields: dict[str, tuple[str, ...] | str | None] | None = None,
    structural_variants: StructuralVariants | None = None,
    gwas_pairs: GwasCatalogPairs | None = None,
    input_lines: list[str] | None = None,
    utrannotator_sweep: UtrannotatorSeedSweep | None = None,
    sequence_edges: SequenceEdges | None = None,
    codon_tables: OrganelleCodonTables | None = None,
    downstream_sequences: DownstreamSequences | None = None,
    transcript_edges: TranscriptEdges | None = None,
    intron_ends: IntronEnds | None = None,
    kept: KeptBack | None = None,
    revel_table: RevelTable | None = None,
    splice_windows: SpliceWindows | None = None,
    spliceai_table: SpliceAiTable | None = None,
) -> dict:
    """Compares expected against actual entries and returns the format's report: the
    one-sided key counts, every field with its differing count and examples grouped by the
    expected consequence set, the stale documented divergences, the order of the shared
    keys, the floors, and the ``verdict``. ``order_only`` names the fields compared as
    multisets with their list separators (``order_only_fields``), each agreement counted in
    ``documented.applied`` under the class; ``floor_fields`` maps a floor's column to the
    fields ``populated_under`` counts it by in a format whose keys differ from the column
    names, to ``None`` for a column the JSON format has no row count for, or to the note of
    a floor the format cannot count for another reason; a column absent from the map, and
    every column of a format without one, is counted by its own name over every entry.
    ``codon_tables`` is the candidate cache's codon tables, which
    ``organelle_peptide_paths_translation_table`` reads.
    ``downstream_sequences`` is the candidate cache's transcripts, which the classes of
    ``TRANSCRIPT_READING_CLASSES`` read.
    ``transcript_edges`` is the candidate cache's exon spans as
    ``hgvs_clipped_reference_against_whole_alternate`` counts them (``TranscriptEdges``).
    ``INDIVIDUAL_ALLELE_ORDER`` in
    ``documented.masks`` compares an individual's alleles as a set: the alternate rows of one
    copy may come in any order within a feature, its allele strings agree on the reference
    allele and the set of alternate alleles, and its ``ALLELE_NUM`` values differ by one
    permutation of the alternates' numbers applied to every row of the copy; each agreement
    is counted under the class, the reordered entries under ``counts``.
    ``sequence_edges`` is what ``HGVS_SHIFT_SEQUENCE_EDGE_WINDOW`` reads as a field class, the
    reference sequence and the candidate cache's transcripts; without it the class applies
    from a manifest alone.
    ``MOTIF_SCORE_REVERSE_STRAND_MNV`` in ``documented.masks`` accepts a differing
    MOTIF_SCORE_CHANGE on a reverse-strand motif entry whose allele is two or more bases,
    counted under the class. ``MERGED_SAME_STABLE_ID_COLLAPSE`` in ``documented.masks``
    accepts a Transcript entry the candidate alone writes, of a model the reference wrote no row
    for on a record it wrote a Transcript row for, its SOURCE RefSeq or absent, once per record
    and allele, counted under the class, and leaves such a record's ``most_severe_consequence``
    uncompared. Each class of ``TERM_PAIR_CLASSES`` in ``documented.masks`` accepts the
    Consequence of a pair its predicate admits, the row counted once and its other fields
    compared, and leaves such a record's ``most_severe_consequence`` uncompared. On a row of a
    class of ``TERM_PAIR_FOLLOWING_CLASSES`` the cells ``term_pair_following_cell`` accepts are
    counted under that class as well. The classes
    of ``structural_variants.classes`` (``STRUCTURAL_VARIANT_CLASSES``) are decided by the
    structural-variant comparator's predicates over the entries' tuples: a
    ``breakend_mate_local_read`` pair is skipped whole as a registry pair is, a
    ``reference_skipped_record`` row and a ``bnd_own_chromosome_rows_lost`` row are
    candidate-only rows counted under their class, a skipped record's record entry is counted
    as a manifest's skipped record's is, and the ``most_severe_consequence`` of a record
    carrying a pair or a lost own-end row is left uncompared. ``intron_ends`` is the candidate
    cache's intron ends, without which ``hgvsp_intron_end_insertion`` and
    ``hgvsp_intron_start_insertion`` accept nothing. ``splice_windows`` is the candidate
    cache's introns, without which ``maxentscan_overlapping_windows`` accepts nothing. A key
    the manifest states as equal on
    both engines in a run with a reference sequence (``Documented.equal_with_reference``) is
    left to no built-in row class, and a difference of its Consequence is accepted by no
    class, pair or rule. Where no entry of either side names a Location or an input record
    the lists number, the manifest items are taken out of the applied counts
    (``Documented.leave_unmet_items_uncounted``).
    ``kept`` is what a front end kept back (``hand_over``), the entries then being the
    remainder's: its rows are added to the entry counts, the compared pairs and the populated
    columns, and the place where the order of the shared keys differs counts the kept-back
    rows before it."""
    floors = floors or {}
    field_pairs = field_pairs or {}
    order_only = order_only or {}
    floor_fields = floor_fields or {}
    if spliceai_table is not None:
        spliceai_table.refusals.clear()
    actual_by_key = {e.key: e for e in actual}
    expected_by_key = {e.key: e for e in expected}
    if (expected or actual) and all(
        e.location is None and (e.record is None or documented.foreign_ordinals) for e in (*expected, *actual)
    ):
        documented.leave_unmet_items_uncounted()

    def manifest_reads(e: Entry) -> bool:
        """Whether the manifest speaks for the entry's key, so that no built-in row class
        reads the row: a documented key, or one stated as equal with a reference sequence."""
        return documented.divergent(e) is not None or documented.equal_with_reference(e)

    missing: list[Entry] = []
    unexpected: list[Entry] = []
    mismatches: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    mismatch_totals: Counter[str] = Counter()
    pairs_accepted: Counter[str] = Counter()
    nearest_accepted: Counter[str] = Counter()
    order_only_accepted: Counter[str] = Counter()
    order_only_examples: dict[str, list[str]] = defaultdict(list)
    allele_order = INDIVIDUAL_ALLELE_ORDER in documented.masks
    allele_order_accepted: Counter[str] = Counter()
    allele_order_examples: dict[str, list[str]] = defaultdict(list)
    reverse_motif = MOTIF_SCORE_REVERSE_STRAND_MNV in documented.masks
    reverse_motif_accepted: Counter[str] = Counter()
    reverse_motif_examples: dict[str, list[str]] = defaultdict(list)
    # The built-in classes of the corrected behaviours: the registry pairs are found before
    # the entries are compared, since a record entry's summary may precede its pair; a field
    # class is credited once per entry it accepts a field of, whatever the format splits the
    # row into, so its applied count is the same in every format (a class that reads a record's
    # alleles leaves a default or tab row under an identifier, ``RecordRows.one_record``).
    registry = {k for k in REGISTRY_CLASSES if k in documented.masks}
    field_classes = [k for k in FIELD_CLASSES if k in documented.masks]
    record_index = RecordRows(expected, input_lines) if any(k in RECORD_READING_CLASSES for k in field_classes) else None
    candidate_index = RecordRows(actual, input_lines) if any(k in SIBLING_ALLELE_CLASSES for k in field_classes) else None
    # A VCF or JSON consequence entry's span is on its record's entry (``indel_span``).
    record_entries = {
        e.record: e for e in expected if e.record is not None and e.key.partition("#")[0].endswith("|record")
    } if sequence_edges is not None else {}
    registry_pairs: dict[str, str] = {}
    registry_records: set[int] = set()
    if registry:
        for e in expected:
            a = actual_by_key.get(e.key)
            if a is None or e.key.endswith("|record") or manifest_reads(e):
                continue
            klass = registry_class(consequence_set(e.consequence), consequence_set(a.consequence), e.allele, registry)
            if klass is not None:
                registry_pairs[e.key] = klass
                if e.record is not None:
                    registry_records.add(e.record)
    class_examples: dict[str, list[str]] = defaultdict(list)
    class_entries: dict[str, set[str]] = defaultdict(set)
    class_columns: dict[str, Counter[str]] = defaultdict(Counter)
    # The built-in row class ``merged_same_stable_id_collapse``: a Transcript entry the candidate
    # alone writes, its SOURCE RefSeq or absent, on a record the reference wrote a Transcript
    # row for, of a model the reference wrote no row for on that record (a record it skipped
    # whole, or wrote an intergenic row for because every model it overlaps is a lost copy, is
    # another shape; so are another allele's row and a duplicate of a row of a model it kept,
    # since the collapse drops a model's rows and never adds one). Found before the entries
    # are compared: only the pass over the actual entries reaches a candidate-only row, and the
    # summary of a record carrying one is left uncompared where the record entries are compared,
    # before that pass.
    collapsed_rows: set[str] = set()
    collapsed_records: set[int] = set()
    if MERGED_SAME_STABLE_ID_COLLAPSE in documented.masks:

        def transcript_entry(e: Entry) -> bool:
            return feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) == "transcript"

        kept_models = {(record_label(e), e.feature) for e in expected if transcript_entry(e)}
        transcript_labels = {label for label, _ in kept_models}
        # A lost copy's row once per record and allele: a repeat of one is a duplicate, whatever
        # the format's key spells (the text key omits the Uploaded_variation the label carries).
        accepted_rows: set[tuple[str, str, str]] = set()
        for a in actual:
            if a.key in expected_by_key or documented.extra(a) or not transcript_entry(a):
                continue
            label = record_label(a)
            if label not in transcript_labels or (label, a.feature) in kept_models:
                continue
            source = next((a.fields[f] for f in SOURCE_FIELDS if f in a.fields), ABSENT)
            row = (label, a.allele, a.feature)
            if (blank(source) or source == REFSEQ_SOURCE) and row not in accepted_rows:
                accepted_rows.add(row)
                collapsed_rows.add(a.key)
                if a.record is not None:
                    collapsed_records.add(a.record)
    # The built-in row classes of ``TERM_PAIR_CLASSES``: the pairs are found before the entries
    # are compared, as the registry pairs are, since a JSON record's summary precedes its
    # consequence entries and its ``most_severe_consequence`` follows the term.
    term_pair_classes = {k: predicate for k, predicate in TERM_PAIR_CLASSES.items() if k in documented.masks}
    term_pair_rows: dict[str, str] = {}
    term_pair_records: set[int] = set()
    if term_pair_classes:
        for e in expected:
            a = actual_by_key.get(e.key)
            if a is None or e.key.endswith("|record") or manifest_reads(e):
                continue
            klass = next((k for k, predicate in term_pair_classes.items() if predicate(e, a)), None)
            if klass is not None:
                term_pair_rows[e.key] = klass
                if e.record is not None:
                    term_pair_records.add(e.record)
    # The built-in structural-variant classes: the structural-variant comparator's predicates
    # over the entries read as its tuples, found before the entries are compared as the registry
    # pairs and the collapsed rows are. A mate-read pair joins the registry pairs, being skipped
    # whole in the same way; a skipped record's rows and a lost own-end row are candidate-only
    # rows, set aside where the collapsed rows are. A key a manifest documents is the manifest's.
    skipped_rows: set[str] = set()
    skipped_records: set[int] = set()
    own_rows_lost: set[str] = set()
    own_rows_records: set[int] = set()
    if structural_variants is not None and structural_variants.classes:
        csc = structural_variants.csc
        expected_tuples = structural_variants.tuples(expected)
        actual_tuples = structural_variants.tuples(actual)
        perl_tuples, rust_tuples = set(expected_tuples.values()), set(actual_tuples.values())

        def candidate_only_rows(excluded: set[StructuralTuple], rows: set[str], records: set[int]) -> None:
            for a in actual:
                if a.key in expected_by_key or documented.extra(a) or actual_tuples.get(a.key) not in excluded:
                    continue
                rows.add(a.key)
                if a.record is not None:
                    records.add(a.record)

        if REFERENCE_SKIPPED_RECORD in structural_variants.classes:
            excluded, _records = csc.filter_reference_skipped_records(
                perl_tuples, rust_tuples, structural_variants.input_variants
            )
            candidate_only_rows(excluded, skipped_rows, skipped_records)
        if BREAKEND_MATE_LOCAL_READ in structural_variants.classes:
            excluded_rust, excluded_perl, _kinds = csc.filter_breakend_mate_local_read_pairs(
                perl_tuples, rust_tuples, structural_variants.models
            )
            for e in expected:
                a = actual_by_key.get(e.key)
                if a is None or e.key in registry_pairs or manifest_reads(e):
                    continue
                if expected_tuples.get(e.key) in excluded_perl and actual_tuples.get(a.key) in excluded_rust:
                    registry_pairs[e.key] = BREAKEND_MATE_LOCAL_READ
                    if e.record is not None:
                        registry_records.add(e.record)
        if BND_OWN_CHROMOSOME_ROWS_LOST in structural_variants.classes:
            excluded = csc.filter_bnd_own_rows_lost(
                perl_tuples, rust_tuples, structural_variants.by_chr, structural_variants.spans
            )
            candidate_only_rows(excluded, own_rows_lost, own_rows_records)
    stale: list[str] = []
    compared = skipped_documented = annotated_without_reference = registry_set_aside = 0

    def field_class_accepting(e: Entry, a: Entry, name: str, ev: str, av: str) -> set[str]:
        """The built-in field classes that accept the differing field, a field two classes
        accept being the first's in ``FIELD_CLASSES`` order: in the ``Extra`` text, one value
        carrying every pair, each differing pair must be accepted by one of them, agree
        under ``order_only`` or, a ``NEAREST`` pair, agree by its start sites
        (``nearest_sites``), and each class used is credited. The class an accepted
        Downstream value needs named beside its own (``downstream_class``) is recorded in
        ``documented.needed``."""
        if not field_classes:
            return set()
        record_rows = record_index.of(e) if record_index is not None else None
        candidate_rows = candidate_index.of(a) if candidate_index is not None else None
        needed: set[str] = set()
        # A row a class of ``TERM_PAIR_FOLLOWING_CLASSES`` took reads its cells under that
        # class before every other, and such a class reads the cells of no other row.
        row_class = term_pair_rows.get(e.key)
        following = row_class if row_class in TERM_PAIR_FOLLOWING_CLASSES else None

        def first_accepting(key: str, expected_value: str, actual_value: str) -> str | None:
            if following is not None and field_class_accepts(following, e, key, expected_value, actual_value, a):
                return following
            for klass in field_classes:
                if klass in TERM_PAIR_FOLLOWING_CLASSES:
                    continue
                if klass == HGVS_SHIFT_SEQUENCE_EDGE_WINDOW:
                    record = record_entries.get(e.record)
                    if sequence_edges is not None and sequence_edges.accepts(e, a, record, key):
                        return klass
                    continue
                if field_class_accepts(
                    klass, e, key, expected_value, actual_value, a,
                    record_rows=record_rows,
                    gwas_pairs=gwas_pairs,
                    utrannotator_sweep=utrannotator_sweep,
                    codon_tables=codon_tables,
                    downstream_sequences=downstream_sequences,
                    transcript_edges=transcript_edges,
                    candidate_rows=candidate_rows,
                    declared=field_classes,
                    intron_ends=intron_ends,
                    revel_table=revel_table,
                    spliceai_table=spliceai_table,
                    assembly=documented.assembly,
                    splice_windows=splice_windows,
                ):
                    if klass == PLUGIN_HGVS_SHIFTED_CDS_COORDINATES and not field_named(key, NMD_FIELDS):
                        decided = downstream_class(e, a, key, expected_value, actual_value, downstream_sequences, record_rows)
                        if decided is not None and decided[1] is not None:
                            needed.add(decided[1])
                    return klass
            return None

        if name == "Extra":
            used: set[str] = set()

            def accept_pair(key: str, expected_value: str, actual_value: str) -> bool:
                klass = first_accepting(key, expected_value, actual_value)
                if klass is not None:
                    used.add(klass)
                    return True
                if key == "NEAREST" and nearest_sites is not None and nearest_sites.agree(expected_value, actual_value):
                    return True
                return key in order_only and same_items(expected_value, actual_value, order_only[key])

            accepted = used if extra_class_accepted(ev, av, accept_pair) else set()
        else:
            klass = first_accepting(name, ev, av)
            accepted = set() if klass is None else {klass}
        if accepted:
            documented.needed |= needed
        return accepted

    def carries_individual(e: Entry) -> bool:
        """Whether the entry itself carries one of ``INDIVIDUAL_FIELDS`` with a value."""
        return any(dash(e.fields.get(f, "")) != "-" for f in INDIVIDUAL_FIELDS)

    individual_records = {
        e.key.split("|", 1)[0] for e in expected if e.key.endswith("|record") and carries_individual(e)
    }

    def is_individual(e: Entry) -> bool:
        """An entry of an individual's copy of a record: it carries IND itself (a default, tab
        or CSQ entry) or its record object carries ind (a JSON consequence entry)."""
        return carries_individual(e) or e.key.split("|", 1)[0] in individual_records

    def copy_of(e: Entry) -> tuple[str, ...]:
        """The copy an entry belongs to: its record label, name and individual."""
        individual = tuple(e.fields.get(f, "") for f in INDIVIDUAL_FIELDS)
        return (e.key.split("|", 1)[0], e.fields.get("Uploaded_variation", ""), *individual)

    def copy_feature(e: Entry) -> tuple[str, ...]:
        """The copy and feature an entry belongs to: ``copy_of`` with its feature and its kind."""
        return (*copy_of(e), e.feature, e.key.rsplit("|", 1)[-1].partition("#")[0])

    def allele_number(e: Entry) -> str | None:
        """The entry's allele number, under either field name, or ``None``."""
        return next((e.fields[name] for name in ALLELE_NUMBER_FIELDS if name in e.fields), None)

    def allele_numbers(entries: list[Entry]) -> dict[tuple[str, ...], dict[str, set[str]]]:
        """Per copy, the allele numbers each of its alleles carries across the copy's rows."""
        numbers: dict[tuple[str, ...], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
        for e in entries:
            if not is_individual(e):
                continue
            number = allele_number(e)
            if number is not None:
                numbers[copy_of(e)][e.allele].add(number)
        return numbers

    expected_numbers = allele_numbers(expected) if allele_order else {}
    actual_numbers = allele_numbers(actual) if allele_order else {}

    def renumbered_alternates(e: Entry, a: Entry, expected_value: str, actual_value: str) -> bool:
        """Whether the row's allele numbers, ``expected_value`` and ``actual_value``, differ by
        one permutation of the alternates' numbers applied to every row of the copy: each
        allele carries one number on each side, the row's own among them, the two sides number
        the same alleles with the same numbers, and the allele numbered 0, the reference, is
        the same allele on both (Parser/VCF.pm:639 keeps the reference first)."""
        expected_map = expected_numbers.get(copy_of(e), {})
        actual_map = actual_numbers.get(copy_of(a), {})
        if not expected_map or expected_map.keys() != actual_map.keys():
            return False
        if any(len(numbers) != 1 for numbers in (*expected_map.values(), *actual_map.values())):
            return False
        expected_number_of = {allele: next(iter(numbers)) for allele, numbers in expected_map.items()}
        actual_number_of = {allele: next(iter(numbers)) for allele, numbers in actual_map.items()}
        if expected_number_of.get(e.allele) != expected_value or actual_number_of.get(a.allele) != actual_value:
            return False
        if sorted(expected_number_of.values()) != sorted(actual_number_of.values()):
            return False
        expected_reference = [allele for allele, number in expected_number_of.items() if number == "0"]
        actual_reference = [allele for allele, number in actual_number_of.items() if number == "0"]
        return expected_reference == actual_reference

    def reference_row(e: Entry) -> bool:
        """Whether the entry is identifiably the copy's reference-allele row: numbered 0, or
        its allele is the reference allele it prints (``-`` for an insertion's); a row
        carrying neither field is not."""
        number = allele_number(e)
        if number is not None:
            return number == "0"
        reference = next((dash(e.fields[name]) for name in REFERENCE_ALLELE_FIELDS if name in e.fields), None)
        return reference is not None and reference == e.allele

    def same_alleles(expected_value: str, actual_value: str) -> bool:
        """The same reference allele followed by the same alternate alleles in any order."""
        if expected_value in ("", "-", ABSENT) or actual_value in ("", "-", ABSENT):
            return expected_value == actual_value
        ref, _, alts = expected_value.partition("/")
        actual_ref, _, actual_alts = actual_value.partition("/")
        return ref == actual_ref and same_items(alts, actual_alts, "/")

    def allele_order_accepts(e: Entry, a: Entry, name: str, ev: str, av: str) -> bool:
        """Whether a field of an individual's entry differs in the order of the copy's
        alleles alone; in the Extra text every differing pair must."""
        if not (allele_order and is_individual(e)):
            return False

        def reordered_only(key: str, expected_value: str, actual_value: str) -> bool:
            if key in ALLELE_STRING_FIELDS:
                return same_alleles(expected_value, actual_value)
            if key in ALLELE_NUMBER_FIELDS:
                return renumbered_alternates(e, a, expected_value, actual_value)
            return key in order_only and same_items(expected_value, actual_value, order_only[key])

        if name == "Extra":
            return extra_accepted(ev, av, reordered_only)
        return name in ALLELE_STRING_FIELDS + ALLELE_NUMBER_FIELDS and reordered_only(name, ev, av)

    def reverse_motif_score_accepts(e: Entry, name: str, ev: str, av: str) -> bool:
        """Whether the differing field is a MOTIF_SCORE_CHANGE the candidate writes on a
        reverse-strand MotifFeature entry whose allele is two or more bases
        (``MOTIF_SCORE_REVERSE_STRAND_MNV``); a score the candidate leaves out where the
        reference writes one is a difference, and in the Extra text that pair alone may
        differ."""
        # The Allele column alone does not tell an insertion's bases from a substitution's, so
        # an insertion's row reads as a multi-base allele here.
        if not (reverse_motif and re.fullmatch(r"[ACGT]{2,}", e.allele)):
            return False
        if feature_kind(e.key.rsplit("|", 1)[-1].partition("#")[0]) != "motif_feature":
            return False
        if not any(e.fields.get(f) == "-1" for f in STRAND_FIELDS):
            return False

        def well_formed(score: str) -> bool:
            return MOTIF_SCORE_TEXT.fullmatch(score) is not None and abs(float(score)) <= 1

        if name == "Extra":
            if extra_accepted(ev, av, lambda k, v, a: k == "MOTIF_SCORE_CHANGE" and well_formed(a)):
                return True
            # Ensembl VEP writes no score where the allele runs past the motif's end; the candidate writes one.
            actual_pairs = av.split(";")
            without_score = [kv for kv in actual_pairs if not kv.startswith("MOTIF_SCORE_CHANGE=")]
            scores = [kv.partition("=")[2] for kv in actual_pairs if kv.startswith("MOTIF_SCORE_CHANGE=")]
            expected_scored = "MOTIF_SCORE_CHANGE=" in ev
            return len(scores) == 1 and well_formed(scores[0]) and not expected_scored and ";".join(without_score) == ev
        return name in MOTIF_SCORE_FIELDS and well_formed(av)

    def count_difference(e: Entry, name: str, ev: str, av: str) -> None:
        """Counts the differing field ``name`` of the entry, an example kept under the
        reference's consequence set."""
        mismatch_totals[name] += 1
        group = mismatches[name][e.consequence]
        if len(group) < examples:
            group.append(f"{e.key}: expected {ev!r}, got {av!r}")

    def order_key(e: Entry) -> str:
        """The key the order comparison reads: an individual's consequence entry by its copy
        and feature (``copy_feature``) and whether it is the reference row, without its
        allele and duplicate counter, since the copy's alternate rows come in the order of
        its alleles while the reference row keeps its place."""
        if not allele_order or e.key.endswith("|record") or not is_individual(e):
            return e.key
        return "|".join([*copy_feature(e), "reference" if reference_row(e) else "alternate"])

    for e in expected:
        a = actual_by_key.get(e.key)
        sets = documented.divergent(e)
        if sets is not None:
            skipped_documented += 1
            if a is None:
                if sets:
                    stale.append(f"{e.key}: documented candidate sets {sets}, but the candidate emits no row")
            else:
                got = consequence_set(a.consequence)
                # A VCF line and a JSON object print no consequence set of their own, so a
                # record-level row the named run documents for the record (the one row of a
                # `--summary` run) has nothing to observe on its record entry.
                record_entry = documented.run is not None and e.key.endswith("|record") and not e.consequence and not got
                if not record_entry and not any(consequence_set(s) == got for s in sets):
                    stale.append(
                        f"{e.key}: documented candidate sets {sets}, the candidate now emits "
                        f"{a.consequence!r} (expected {e.consequence!r})"
                    )
            continue
        if a is None:
            missing.append(e)
            continue
        registry_klass = registry_pairs.get(e.key)
        if registry_klass is not None:
            registry_set_aside += 1
            documented.applied[registry_klass] += 1
            if len(class_examples[registry_klass]) < examples:
                class_examples[registry_klass].append(f"{e.key}: expected {e.consequence!r}, got {a.consequence!r}")
            continue
        compared += 1
        record_summary_differs = (
            e.key.endswith("|record")
            and e.record is not None
            and (
                documented.record_divergent(e.record)
                or e.record in registry_records
                or e.record in collapsed_records
                or e.record in term_pair_records
                or e.record in own_rows_records
            )
        )
        # The default format's Extra is the raw text of the fields compared one by one
        # below, so it differs wherever one of them is documented to.
        extra_differs = documented.has_divergent_fields(e)
        stated_equal = documented.equal_with_reference(e)
        for name in sorted(set(e.fields) | set(a.fields)):
            if record_summary_differs and name == "most_severe_consequence":
                continue
            if extra_differs and name == "Extra":
                continue
            ev = e.fields.get(name, ABSENT)
            av = a.fields.get(name, ABSENT)
            documented_value = documented.divergent_field(e, name)
            if documented_value is not None:
                # The JSON format writes a numeric plugin value as a number, so the
                # documented default-format text is also matched with the trailing zeros
                # of its fractional part dropped (`1311092.0000` is written `1311092`).
                printed, recorded = documented_form(av), documented_form(documented_value)
                shortened = recorded.rstrip("0").rstrip(".") if re.fullmatch(r"-?\d+\.\d+", recorded) else recorded
                if printed not in (recorded, shortened):
                    stale.append(
                        f"{e.key}: field {name} documented as {documented_value!r}, the candidate "
                        f"now prints {av!r} (expected {ev!r})"
                    )
                continue
            if ev == av:
                continue
            if stated_equal and field_named(name, CONSEQUENCE_FIELDS):
                count_difference(e, name, ev, av)
                continue
            # ``colocated_list_cut_at_cache_cap``: the reference's value is read as the candidate
            # reads its cut lists and compared on as that reading, so an ``Extra`` text with
            # such a pair beside a pair another class accepts meets both.
            if COLOCATED_LIST_CUT_AT_CACHE_CAP in documented.masks:
                reading = cut_list_reading(e, name, ev, av)
                if reading is not None:
                    klass = COLOCATED_LIST_CUT_AT_CACHE_CAP
                    class_columns[klass][name] += 1
                    class_entries[klass].add(e.key)
                    if len(class_examples[klass]) < examples:
                        class_examples[klass].append(cut_list_example(e, name, ev, av))
                    ev = reading
                    if ev == av:
                        continue
            if (ev, av) in field_pairs.get(name, ()) or (
                name == "Extra" and extra_pairs_accepted(ev, av, field_pairs)
            ):
                pairs_accepted[name] += 1
                continue
            if nearest_sites is not None and (
                (name in NEAREST_FIELDS and nearest_sites.agree(ev, av))
                or (name == "Extra" and extra_nearest_accepted(ev, av, nearest_sites))
            ):
                nearest_accepted[name] += 1
                continue
            if (name in order_only and same_items(ev, av, order_only[name])) or (
                name == "Extra" and order_only and extra_order_only_accepted(ev, av, order_only)
            ):
                order_only_accepted[name] += 1
                documented.applied[ORDER_ONLY] += 1
                if len(order_only_examples[name]) < examples:
                    order_only_examples[name].append(f"{e.key}: expected {ev!r}, got {av!r}")
                continue
            if allele_order_accepts(e, a, name, ev, av):
                allele_order_accepted[name] += 1
                documented.applied[INDIVIDUAL_ALLELE_ORDER] += 1
                if len(allele_order_examples[name]) < examples:
                    allele_order_examples[name].append(f"{e.key}: expected {ev!r}, got {av!r}")
                continue
            if reverse_motif_score_accepts(e, name, ev, av):
                reverse_motif_accepted[name] += 1
                documented.applied[MOTIF_SCORE_REVERSE_STRAND_MNV] += 1
                if len(reverse_motif_examples[name]) < examples:
                    reverse_motif_examples[name].append(f"{e.key}: expected {ev!r}, got {av!r}")
                continue
            if e.key in term_pair_rows and field_named(name, CONSEQUENCE_FIELDS):
                accepting = {term_pair_rows[e.key]}
            else:
                accepting = field_class_accepting(e, a, name, ev, av)
            if accepting:
                for klass in accepting:
                    class_columns[klass][name] += 1
                    class_entries[klass].add(e.key)
                    if len(class_examples[klass]) < examples:
                        class_examples[klass].append(f"{e.key}: {name} expected {ev!r}, got {av!r}")
                continue
            count_difference(e, name, ev, av)
    # The built-in row shape of ``frequency_filter_zero_value``: a record the candidate alone
    # writes, every consequence row of which carries a FREQS of zeros, is set aside with its
    # record entry, as a manifest's documented tuples of the class are.
    zero_frequency_records: set[str] = set()
    if FREQUENCY_FILTER_ZERO_VALUE in documented.masks:
        reference_labels = {record_label(e) for e in expected}
        candidate_only: dict[str, list[Entry]] = defaultdict(list)
        for a in actual:
            if a.key not in expected_by_key and not documented.extra(a) and not a.key.endswith("|record"):
                candidate_only[record_label(a)].append(a)
        # A record the reference wrote any row for is kept by both engines, so an extra
        # candidate row on it is a difference whatever its FREQS says.
        zero_frequency_records = {
            label for label, rows in candidate_only.items()
            if label not in reference_labels and rows and all(freqs_all_zero(a) for a in rows)
        }
    for a in actual:
        if a.key in expected_by_key or documented.extra(a):
            continue
        # Only the record-level line or object is without a counterpart; a consequence
        # entry on a skipped record is a candidate-only key, documented, built in or unexpected.
        if (
            a.key.endswith("|record")
            and a.record is not None
            and (documented.reference_skipped(a.record) or a.record in skipped_records)
        ):
            annotated_without_reference += 1
            continue
        if record_label(a) in zero_frequency_records:
            if a.key.endswith("|record"):
                annotated_without_reference += 1
            else:
                class_entries[FREQUENCY_FILTER_ZERO_VALUE].add(a.key)
                if len(class_examples[FREQUENCY_FILTER_ZERO_VALUE]) < examples:
                    freqs = next((a.fields[f] for f in FREQS_FIELDS if f in a.fields), None)
                    class_examples[FREQUENCY_FILTER_ZERO_VALUE].append(f"{a.key}: candidate-only row, FREQS {freqs!r}")
            continue
        if a.key in collapsed_rows:
            class_entries[MERGED_SAME_STABLE_ID_COLLAPSE].add(a.key)
            if len(class_examples[MERGED_SAME_STABLE_ID_COLLAPSE]) < examples:
                class_examples[MERGED_SAME_STABLE_ID_COLLAPSE].append(f"{a.key}: candidate-only row, {a.consequence!r}")
            continue
        if a.key in skipped_rows:
            class_entries[REFERENCE_SKIPPED_RECORD].add(a.key)
            if len(class_examples[REFERENCE_SKIPPED_RECORD]) < examples:
                class_examples[REFERENCE_SKIPPED_RECORD].append(
                    f"{a.key}: candidate-only row on a record the reference skipped, {a.consequence!r}"
                )
            continue
        if a.key in own_rows_lost:
            class_entries[BND_OWN_CHROMOSOME_ROWS_LOST].add(a.key)
            if len(class_examples[BND_OWN_CHROMOSOME_ROWS_LOST]) < examples:
                class_examples[BND_OWN_CHROMOSOME_ROWS_LOST].append(
                    f"{a.key}: candidate-only own-end row, {a.consequence!r}"
                )
            continue
        unexpected.append(a)

    shared_expected = [e for e in expected if e.key in actual_by_key and documented.divergent(e) is None]
    shared_actual = [a for a in actual if a.key in expected_by_key and documented.divergent(a) is None]
    order = None
    reordered = 0
    if check_order:
        expected_order = [order_key(e) for e in shared_expected]
        actual_order = [order_key(a) for a in shared_actual]
        if expected_order != actual_order:
            first = next((i for i, (x, y) in enumerate(zip(expected_order, actual_order)) if x != y), 0)
            # The front end keeps back no row between two rows of one handed-over group, on
            # either side, and keeps its rows in one order on both, so a shared key of the
            # remainder has the same kept-back rows before it on both sides: the remainder's
            # first difference is the outputs' first, moved by those rows.
            before = kept.before(shared_expected[first].key) if kept is not None and first < len(shared_expected) else 0
            order = {
                "first_divergence": first + before,
                "expected": shared_expected[first].key if first < len(shared_expected) else None,
                "actual": shared_actual[first].key if first < len(shared_actual) else None,
            }
        else:
            reordered = sum(1 for e, a in zip(shared_expected, shared_actual) if e.key != a.key)
            if reordered:
                documented.applied[INDIVIDUAL_ALLELE_ORDER] += reordered

    for klass, keys in class_entries.items():
        documented.applied[klass] += len(keys)

    expected_populated = populated_counts(expected)
    actual_populated = populated_counts(actual)
    kept_entries = 0
    if kept is not None:
        kept_entries = kept.entries
        expected_populated.update(kept.populated)
        actual_populated.update(kept.populated)
    columns: dict[str, dict] = {}
    floors_not_applicable: list[str] = []
    for name in sorted(set(expected_populated) | set(actual_populated) | set(mismatch_totals) | set(floors)):
        expected_count = expected_populated.get(name, 0)
        actual_count = actual_populated.get(name, 0)
        floor = floors.get(name)
        floor_met = None
        floor_detail: dict[str, str | list[str]] = {}
        if floor is not None:
            counted = floor_fields.get(name)
            if name not in floor_fields:
                floor_met = expected_count >= floor and actual_count >= floor
            elif counted is None or isinstance(counted, str):
                floors_not_applicable.append(name)
                floor_detail["floor_note"] = FLOOR_NOT_APPLICABLE if counted is None else counted
            else:
                expected_count = populated_under(expected, counted)
                actual_count = populated_under(actual, counted)
                floor_detail["floor_fields"] = list(counted)
                floor_met = expected_count >= floor and actual_count >= floor
        col = {
            "expected_populated": expected_count,
            "actual_populated": actual_count,
            "differing": mismatch_totals.get(name, 0),
            "pairs_accepted": pairs_accepted.get(name, 0),
            "nearest_accepted": nearest_accepted.get(name, 0),
            "order_only": order_only_accepted.get(name, 0),
            "individual_allele_order": allele_order_accepted.get(name, 0),
            "motif_score_reverse_strand_mnv": reverse_motif_accepted.get(name, 0),
            "floor": floor,
            "floor_met": floor_met,
            **floor_detail,
        }
        if name in mismatches:
            col["examples"] = dict(mismatches[name])
        if name in order_only_examples:
            col["order_only_examples"] = order_only_examples[name]
        if name in allele_order_examples:
            col["individual_allele_order_examples"] = allele_order_examples[name]
        if name in reverse_motif_examples:
            col["motif_score_reverse_strand_mnv_examples"] = reverse_motif_examples[name]
        masked = {klass: counts[name] for klass, counts in sorted(class_columns.items()) if counts.get(name)}
        if masked:
            col["masked"] = masked
        columns[name] = col

    reasons: list[str] = []
    if missing:
        reasons.append(f"missing in the candidate: {len(missing)} entries")
    if unexpected:
        reasons.append(f"unexpected in the candidate: {len(unexpected)} entries")
    for name in sorted(mismatch_totals):
        reasons.append(f"field {name}: {mismatch_totals[name]} differing entries")
    if stale:
        reasons.append(f"documented divergences not observed: {len(stale)}")
    if order is not None:
        reasons.append(f"order differs at shared position {order['first_divergence']}")
    for name in sorted(floors):
        col = columns[name]
        if col["floor_met"] is False:
            reasons.append(
                f"floor {name}={floors[name]}: expected populated {col['expected_populated']}, "
                f"actual populated {col['actual_populated']}"
            )

    def grouped(entries: list[Entry]) -> dict[str, dict]:
        by_cons: dict[str, list[Entry]] = defaultdict(list)
        for e in entries:
            by_cons[e.consequence].append(e)
        return {
            cons: {"count": len(items), "examples": [e.key for e in items[:examples]]}
            for cons, items in sorted(by_cons.items())
        }

    # The rows the class of a record of several entries was asked about under the score files.
    refusals = spliceai_table.refusals.counts() if spliceai_table is not None else None
    return {
        "verdict": "red" if reasons else "green",
        "reasons": reasons,
        "counts": {
            "expected_entries": len(expected) + kept_entries,
            "actual_entries": len(actual) + kept_entries,
            "compared": compared + kept_entries,
            "documented_divergences_skipped": skipped_documented,
            "record_entries_on_reference_skipped_records": annotated_without_reference,
            "individual_entries_reordered": reordered,
            "registry_pairs_set_aside": registry_set_aside,
        },
        "one_sided": {"expected_only": len(missing), "actual_only": len(unexpected)},
        "missing": grouped(missing),
        "unexpected": grouped(unexpected),
        "columns": columns,
        "mask_examples": dict(sorted(class_examples.items())),
        "floors_not_applicable": floors_not_applicable,
        "stale_divergences": stale[:20],
        "order": order,
    } | ({"spliceai_refusals": refusals} if refusals and refusals["asked"] else {})


def is_regulatory_row(e: Entry) -> bool:
    """A RegulatoryFeature or MotifFeature consequence entry in any format: the key's last
    segment is the Feature_type (default, tab, VCF) or the JSON list's kind."""
    last_segment = e.key.rsplit("|", 1)[-1].partition("#")[0]
    return feature_kind(last_segment) in REGULATORY_KINDS


def filter_consequence_rows(
    entries: list[Entry], keep: Callable[[Entry], bool], keep_bare_records: bool
) -> list[Entry]:
    """The entries whose consequence rows ``keep`` accepts, in their order, with the VCF/JSON
    record entries: every one when ``keep_bare_records``, else those of a record that kept a
    row. A record that lost a row has its ``most_severe_consequence`` recomputed from the
    rows it kept, as the candidate computes it from the rows it writes. Mutates the kept
    record entries; raises ``ComparatorError`` on a term ``CONSEQUENCE_RANK`` lacks."""
    terms_by_record: dict[str, set[str]] = defaultdict(set)
    lost_rows: set[str] = set()
    for e in entries:
        if e.key.endswith("|record"):
            continue
        record = e.key.split("|", 1)[0]
        if keep(e):
            terms_by_record[record].update(t for t in re.split(r"[&,]", e.consequence) if t)
        else:
            lost_rows.add(record)
    kept: list[Entry] = []
    for e in entries:
        if not e.key.endswith("|record"):
            if keep(e):
                kept.append(e)
            continue
        record = e.key.split("|", 1)[0]
        if not keep_bare_records and record not in terms_by_record:
            continue
        if record in lost_rows and "most_severe_consequence" in e.fields:
            terms = terms_by_record.get(record, set())
            unranked = sorted(terms.difference(CONSEQUENCE_RANK))
            if unranked:
                raise ComparatorError(f"{e.key}: no rank for the consequence terms {unranked}")
            e.fields = {**e.fields, "most_severe_consequence": most_severe(terms)}
        kept.append(e)
    return kept


def without_regulatory_rows(expected: list[Entry]) -> list[Entry]:
    """The reference entries without their regulatory rows (``regulatory_feature_rows``): a
    record the reference annotates against a regulatory feature and no transcript compares
    as its intergenic row, every record entry stays (``filter_consequence_rows``)."""
    return filter_consequence_rows(expected, lambda e: not is_regulatory_row(e), True)


def only_regulatory_rows(entries: list[Entry]) -> list[Entry]:
    """The regulatory rows of either side and the record entries of the records that have
    one (``--keep-only regulatory_feature_rows``; ``filter_consequence_rows``)."""
    return filter_consequence_rows(entries, is_regulatory_row, False)


# The hand-over of a default or tab comparison to the front end ``tools/vep-formats-remainder``.
# A record group is every row of one side that shares a Location. ``compare_entries`` pairs
# entries by a key that starts with the Location, groups a record's rows under a label that
# starts with it (``record_label``) and reaches a rule only through a field that differs, so
# a group whose rows are the same text in the same order on both sides adds its entries to
# the counts, its fields to the populated columns and its place to the order of the shared
# keys, and nothing else. Three readings take an equal group all the same, and each hands
# its groups over whole: a key the manifest documents is checked for staleness on equal
# values too, ``regulatory_feature_rows`` drops reference rows before the comparison, and
# ``conservation_short_name_spelling`` rewrites the reference text.
FRONT_END_VARIABLE = "VEP_FORMATS_REMAINDER"
FRONT_END_FORMATS = ("default", "tab")
# The classes compared on the whole outputs: the structural-variant classes read both sides'
# rows as one set of tuples, and ``individual_allele_order`` reads the order of the shared
# keys by copy and feature.
FRONT_END_WHOLE_CLASSES = (*STRUCTURAL_VARIANT_CLASSES, INDIVIDUAL_ALLELE_ORDER)
# The Feature_type values whose rows ``regulatory_feature_rows`` is asked about: a reference
# row of a value it keeps stays where it is, a row of any other value hands its group over.
FRONT_END_FEATURE_TYPES = ("Transcript", "-", "")


class FrontEnd:
    """The front end ``FRONT_END_VARIABLE`` names, and ``note``, what the comparison did
    with it."""

    def __init__(self, executable: str) -> None:
        self.executable = executable
        self.note = "the whole outputs were compared: only a default or tab comparison is handed over"

    @classmethod
    def named(cls, value: str | None) -> FrontEnd | None:
        """The front end ``value`` names, a path or a command on ``PATH``; ``None`` for an
        unset or empty value. A value naming no executable is the comparator's error."""
        if not value:
            return None
        executable = shutil.which(value)
        if executable is None:
            raise ComparatorError(f"{FRONT_END_VARIABLE}={value} names no executable")
        return cls(executable)


@dataclass
class KeptBack:
    """What the front end kept back of one comparison: ``expected`` and ``actual`` are the
    remainder files under ``directory``, ``entries`` the rows kept back on each side,
    ``populated`` their count per field (``populated_counts``), and ``runs`` the reference
    remainder's data rows in file order as (rows, kept-back rows before them). ``keys`` are
    the keys of those rows, set once the remainder is parsed."""

    directory: tempfile.TemporaryDirectory
    expected: Path
    actual: Path
    entries: int
    populated: Counter
    runs: list[tuple[int, int]]
    keys: list[str] = field(default_factory=list)

    def before(self, key: str) -> int:
        """The kept-back rows that precede the remainder's row of ``key`` in the reference."""
        row = self.keys.index(key)
        for rows, before in self.runs:
            if row < rows:
                return before
            row -= rows
        return 0


def literal_head(pattern: re.Pattern) -> str | None:
    """The text every match of ``pattern`` starts with, read off the pattern's own text: its
    leading letters, digits and underscores, less the last one when a quantifier follows it.
    ``None`` for a pattern with an alternation or a flag, or without such a head."""
    source = pattern.pattern
    if "|" in source or pattern.flags != re.compile("").flags:
        return None
    head = re.match(r"[A-Za-z0-9_]*", source).group()
    if source[len(head) : len(head) + 1] in ("*", "+", "?", "{"):
        head = head[:-1]
    return head or None


def hand_over(
    front_end: FrontEnd,
    fmt: str,
    expected_path: Path,
    actual_path: Path,
    documented: Documented,
    unordered: bool,
    keep_only: str | None,
) -> KeptBack | None:
    """Runs the front end over one comparison and returns what it kept back, or ``None`` for
    a comparison made on the whole outputs; ``front_end.note`` says which. The front end is
    given the Locations of every key the manifest documents; under ``regulatory_feature_rows``
    the Feature_type values of the rows that class keeps; under
    ``conservation_short_name_spelling`` the text every misspelt key starts with
    (``literal_head``; a pattern that shows none is compared on the whole outputs); and
    whether the order of the shared keys is read. The front end hands every row over where
    that text stands in the reference's column line or in a key column, which the rewrite
    would rename or move, and where a Location holds a space, the separator of
    ``record_label``."""
    whole: str | None = None
    if fmt not in FRONT_END_FORMATS:
        whole = f"the {fmt} format"
    elif keep_only is not None:
        whole = "--keep-only"
    else:
        whole = next((k for k in FRONT_END_WHOLE_CLASSES if k in documented.masks), None)
    misspelt = None
    if whole is None and CONSERVATION_SHORT_NAME_SPELLING in documented.masks:
        misspelt = literal_head(MISSPELT_CONSERVATION_KEY)
        if misspelt is None:
            whole = CONSERVATION_SHORT_NAME_SPELLING
    if whole is not None:
        front_end.note = f"the whole outputs were compared: {whole} is outside the front end"
        return None
    try:
        directory = tempfile.TemporaryDirectory(prefix=f"{COMPONENT}.")
    except OSError as exc:
        raise ComparatorError(f"cannot make the front end's directory: {exc.strerror or exc}") from exc
    work = Path(directory.name)

    def listed(name: str, values: Collection[str]) -> str:
        """The file of ``values``, one a line; a value no line of an output can hold (a tab,
        a line end, a character UTF-8 does not write) is left out."""
        lines = []
        for value in sorted(values):
            try:
                value.encode("utf-8")
            except UnicodeEncodeError:
                continue
            if not re.search(r"[\t\n\r]", value):
                lines.append(value + "\n")
        try:
            (work / name).write_text("".join(lines), encoding="utf-8", newline="")
        except OSError as exc:
            raise ComparatorError(f"cannot write {work / name}: {exc.strerror or exc}") from exc
        return str(work / name)

    tables = (documented.by_location, documented.extra_by_location, documented.fields_by_location)
    command = [
        front_end.executable,
        "--expected", str(expected_path),
        "--actual", str(actual_path),
        "--out-dir", str(work),
        "--force-locations", listed("locations.txt", {key[0] for table in tables for key in table}),
    ]
    if misspelt is not None:
        command += ["--force-text", listed("texts.txt", [misspelt])]
    if REGULATORY_FEATURE_ROWS in documented.masks:
        kept_types = [t for t in FRONT_END_FEATURE_TYPES if feature_kind(t) not in REGULATORY_KINDS]
        command += ["--keep-values", "Feature_type", listed("feature_types.txt", kept_types)]
    if unordered:
        command.append("--unordered")
    try:
        done = subprocess.run(command, capture_output=True, text=True, errors="replace", check=False)
    except OSError as exc:
        raise ComparatorError(f"cannot run {front_end.executable}: {exc.strerror or exc}") from exc
    if done.returncode != 0:
        lines = done.stderr.strip().splitlines()
        raise ComparatorError(f"{front_end.executable} exited {done.returncode}: {lines[-1] if lines else 'no message'}")
    try:
        counters = json.loads((work / "counters.json").read_text(encoding="utf-8"))
        kept = KeptBack(
            directory,
            work / "expected.txt",
            work / "actual.txt",
            counters["kept_entries"],
            Counter(counters["kept_populated"]),
            [(rows, before) for rows, before in counters["kept_before"]],
        )
        rows, groups, handed = counters["expected"]["lines"], counters["groups"], counters["groups_handed_over"]
        reason = counters["whole"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise ComparatorError(f"{front_end.executable} left no counters this comparator reads: {exc}") from exc
    if reason is not None:
        front_end.note = f"every row was handed over by {front_end.executable}: {reason}"
    else:
        front_end.note = (
            f"{kept.entries} of {rows} reference rows were kept back by {front_end.executable}, "
            f"{handed} of {groups} record groups handed over"
        )
    return kept


def compare_files(
    fmt: str,
    expected_path: Path,
    actual_path: Path,
    inputs: list[str] | None,
    documented: Documented,
    masks: set[str],
    floors: dict[str, int],
    field_pairs: dict[str, set[tuple[str, str]]],
    examples: int,
    nearest_sites: NearestSites | None = None,
    unordered: bool = False,
    keep_only: str | None = None,
    structural_variants: StructuralVariants | None = None,
    gwas_pairs: GwasCatalogPairs | None = None,
    utrannotator_sweep: UtrannotatorSeedSweep | None = None,
    sequence_edges: SequenceEdges | None = None,
    codon_tables: OrganelleCodonTables | None = None,
    downstream_sequences: DownstreamSequences | None = None,
    transcript_edges: TranscriptEdges | None = None,
    intron_ends: IntronEnds | None = None,
    front_end: FrontEnd | None = None,
    input_path: Path | None = None,
    revel_table: RevelTable | None = None,
    splice_windows: SpliceWindows | None = None,
    spliceai_table: SpliceAiTable | None = None,
) -> dict:
    """The comparison-mode report for one format. The classes applied are ``documented.masks``
    alone (``masks`` is the same set the caller built ``documented`` from); the
    structural-variant classes among them read ``structural_variants``, which the recoder
    format, a comparison of no annotation rows, leaves out, and
    ``hgvsp_intron_end_insertion`` and ``hgvsp_intron_start_insertion`` read ``intron_ends``,
    and ``maxentscan_overlapping_windows`` reads ``splice_windows``.
    ``unordered``
    compares the shared entries without their order (a run whose reference rows Ensembl VEP
    emits in hash order); a Parquet actual is compared without order in any case.
    ``hgvs_shift_sequence_edge_window`` reads ``sequence_edges``, which the recoder format
    leaves out as well.
    ``organelle_peptide_paths_translation_table`` reads ``codon_tables``, which the recoder
    format leaves out as well.
    The classes of ``TRANSCRIPT_READING_CLASSES`` read ``downstream_sequences``, which the
    recoder format leaves out as well.
    ``keep_only`` (``regulatory_feature_rows``, the one class) keeps the regulatory rows of
    both sides and the record entries of the records that have one (``only_regulatory_rows``),
    and is refused beside the ``regulatory_feature_rows`` mask, which drops them. A floor
    is read where the format carries its column: under ``conservation_short_name_spelling``
    a floor spelt with the misspelt key is read under the renamed one; in the JSON format
    each floor is placed by ``json_floor_location``; in the tab and VCF formats a floor on a
    column neither side's header declares (``declared_columns``) is not applicable
    (``FLOOR_UNDECLARED_NOT_APPLICABLE``), since every row of those formats carries every
    column its header names.
    ``front_end`` hands a default or tab comparison over (``hand_over``): the entries are then
    parsed from the remainder and ``compare_entries`` adds what was kept back.
    ``input_path`` names the input VCF under a class of ``RECORD_REFERENCE_CLASSES``: the
    default and tab entries the class can read take their record's reference bases from it
    (``read_entry_references``)."""
    if keep_only is not None and keep_only != REGULATORY_FEATURE_ROWS:
        raise ComparatorError(f"unknown --keep-only class: {keep_only}")
    if keep_only == REGULATORY_FEATURE_ROWS and REGULATORY_FEATURE_ROWS in documented.masks:
        raise ComparatorError("--keep-only regulatory_feature_rows conflicts with --masks regulatory_feature_rows")
    kept_back = None
    if front_end is not None:
        kept_back = hand_over(front_end, fmt, expected_path, actual_path, documented, unordered, keep_only)
    if fmt == "recoder":
        if keep_only is not None or REGULATORY_FEATURE_ROWS in documented.masks:
            raise ComparatorError("the recoder format has no regulatory rows to keep or drop")
        expected = parse_recoder(read_text(expected_path), documented)
        actual = parse_recoder(read_text(actual_path), documented)
    else:
        expected_format = "tab-canonical" if fmt == "parquet" else fmt
        expected = parse_format(
            expected_format,
            expected_path if kept_back is None else kept_back.expected,
            inputs,
            lambda text: reference_text_masked(text, documented),
        )
        actual = parse_format(fmt, actual_path if kept_back is None else kept_back.actual, inputs)
    if kept_back is not None:
        kept_back.keys = [e.key for e in expected]
        kept_back.directory.cleanup()
    if REGULATORY_FEATURE_ROWS in documented.masks:
        kept = without_regulatory_rows(expected)
        # A class that dropped nothing is unused, as a manifest class matching no item is.
        if len(kept) < len(expected):
            documented.applied[REGULATORY_FEATURE_ROWS] += len(expected) - len(kept)
        expected = kept
    if keep_only == REGULATORY_FEATURE_ROWS:
        expected = only_regulatory_rows(expected)
        actual = only_regulatory_rows(actual)
    if input_path is not None:
        read_entry_references(expected, actual, input_path)

    def renamed_key(name: str) -> str:
        """The column name as the parsed reference spells it: under
        ``conservation_short_name_spelling`` the misspelt conservation key renamed."""
        if CONSERVATION_SHORT_NAME_SPELLING in documented.masks:
            return MISSPELT_CONSERVATION_KEY.sub(CONSERVATION_KEY_SPELLING, name)
        return name

    floors = {renamed_key(name): minimum for name, minimum in floors.items()}
    floor_fields: dict[str, tuple[str, ...] | str | None] | None = None
    if fmt == "json":
        floor_fields = {name: json_floor_location(name, expected, actual) for name in floors}
    elif fmt in ("tab", "vcf"):
        declared = {renamed_key(name) for name in declared_columns(fmt, expected_path) | declared_columns(fmt, actual_path)}
        floor_fields = {name: FLOOR_UNDECLARED_NOT_APPLICABLE for name in floors if name not in declared}
    report = compare_entries(
        expected,
        actual,
        documented,
        floors=floors,
        field_pairs=field_pairs,
        check_order=fmt != "parquet" and not unordered,
        examples=examples,
        nearest_sites=nearest_sites,
        order_only=order_only_fields(fmt) if ORDER_ONLY in documented.masks else None,
        floor_fields=floor_fields,
        structural_variants=structural_variants if fmt != "recoder" else None,
        gwas_pairs=gwas_pairs,
        input_lines=inputs,
        utrannotator_sweep=utrannotator_sweep,
        sequence_edges=sequence_edges if fmt != "recoder" else None,
        codon_tables=codon_tables if fmt != "recoder" else None,
        downstream_sequences=downstream_sequences if fmt != "recoder" else None,
        transcript_edges=transcript_edges,
        intron_ends=intron_ends,
        kept=kept_back,
        revel_table=revel_table if fmt != "recoder" else None,
        splice_windows=splice_windows if fmt != "recoder" else None,
        spliceai_table=spliceai_table if fmt != "recoder" else None,
    )
    report.update(
        {
            "mode": "compare",
            "format": fmt,
            "expected": str(expected_path),
            "actual": str(actual_path),
            "masks": {
                "requested": sorted(documented.masks),
                "applied": dict(sorted(documented.applied.items())),
                "unused": sorted(documented.masks - set(documented.applied) - documented.needed),
                "available": dict(sorted(documented.available.items())),
            },
            "field_pairs": {name: sorted(pairs) for name, pairs in sorted(field_pairs.items())},
            "nearest_sites": None if nearest_sites is None else len(nearest_sites.sites),
            "json_cache_transcripts": None if structural_variants is None else len(structural_variants.spans),
            "unordered": unordered,
            "keep_only": keep_only,
            "manifest_run": documented.run,
            "reference_sequence": documented.reference_sequence,
        }
    )
    return report


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def perl_numify(s: str) -> str:
    """``0 + $s``: the leading number of a string, 0 when it has none."""
    m = PERL_NUMBER_RE.match(s)
    if not m:
        return "0"
    text = m.group(1)
    try:
        value = float(text)
    except ValueError:
        return "0"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    if value == int(value) and "e" not in text.lower():
        return str(int(value))
    return _json_float(value)


def perl_numeric_or_text(word: str) -> int | str:
    """A ``\\w+`` capture as ``numberify`` leaves it: ASCII digits alone become an integer (Perl's
    ``+ 0`` on a string ``looks_like_number`` accepts); anything else stays text."""
    return int(word) if word.isascii() and word.isdigit() else word


def format_coords(start, end) -> str:
    """``Utils.pm:141-166``: ``start``, ``start-end``, or the pair low to high."""
    if start is not None and end is not None:
        try:
            s, e = float(start), float(end)
        except (TypeError, ValueError):
            return f"{start}-{end}"
        if s > e:
            return f"{end}-{start}"
        if s == e:
            return f"{start}"
        return f"{start}-{end}"
    if start is not None:
        return f"{start}-?"
    if end is not None:
        return f"?-{end}"
    return "-"


def uploaded_variation(obj: dict) -> str:
    """``Uploaded_variation`` for a record of VCF input as ``OutputFactory.pm:883-930``
    derives it from the JSON object's ``id``, ``seq_region_name``, ``start``,
    ``allele_string`` and ``input``: the ID when it is one, else ``chr_pos_ref/alt``."""
    chrom = json_scalar(obj.get("seq_region_name", ""))
    start = obj.get("start")
    allele_string = json_scalar(obj.get("allele_string", ""))
    try:
        var_name = f"{chrom}_{int(float(start))}_{allele_string}"
    except (TypeError, ValueError):
        var_name = f"{chrom}_{start}_{allele_string}"
    vn = obj.get("id")
    vn = json_scalar(vn) if vn is not None else None
    if vn not in (None, "", ".") and vn != var_name:
        return vn
    line = obj.get("input")
    cols = line.split("\t") if isinstance(line, str) else []
    if len(cols) >= 5:
        chrom_in, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
        alt = alt.replace(",", "/")
        return "_".join([chrom_in, perl_numify(pos), f"{ref}/{alt}" if alt != "" else ref])
    if vn not in (None, ".") and vn != var_name:
        return vn
    return f"{chrom}_{start}_{allele_string}"


def record_location(obj: dict) -> str:
    """``Location`` as ``OutputFactory.pm:931`` prints it from the JSON object's
    ``seq_region_name``, ``start`` and ``end``."""
    return f"{json_scalar(obj.get('seq_region_name', ''))}:{format_coords(obj.get('start'), obj.get('end'))}"


def vcf_escape(s: str) -> str:
    """``VCF.pm:401-404`` on a value already joined for the default format."""
    s = s.replace(",", "&").replace(";", "%3B")
    s = re.sub(r"\s+", "_", s)
    return s.replace("|", "&")


def reverse_complement(allele: str) -> str:
    """The reverse complement of a sequence allele; any other allele unchanged."""
    if allele and re.fullmatch(r"[ACGTNacgtn-]+", allele):
        return allele.translate(COMPLEMENT)[::-1]
    return allele


def values_equal(a: str, b: str) -> bool:
    """Equality of two scalar texts, numerically when both are numbers."""
    if a == b:
        return True
    try:
        return float(a) == float(b)
    except ValueError:
        return False


def lists_equal(a: list[str], b: list[str]) -> bool:
    return len(a) == len(b) and all(values_equal(x, y) for x, y in zip(a, b))


def default_list(value: str) -> list[str] | None:
    """A default-format value as the list ``convert_arrayref`` joined, percent-encoding
    undone; ``None`` for ``-``."""
    if value in ("-", ABSENT):
        return None
    return decode_extra(value).split(",")


def json_list(value) -> list[str] | None:
    """A JSON value as a list of scalar texts, a string split on ``,`` as the default
    format's list would be; ``None`` when absent."""
    if value is None:
        return None
    if isinstance(value, list):
        return [json_scalar(x) for x in value]
    return json_scalar(value).split(",")


def most_severe(terms: set[str]) -> str:
    """``JSON.pm:387-388``: the lowest-ranked term, ``?`` when there is none."""
    ranked = sorted((t for t in terms if t != "?"), key=lambda t: (CONSEQUENCE_RANK.get(t, 10**6), t))
    return ranked[0] if ranked else "?"


class SelfParity:
    """One engine's outputs on one input, every consequence entry keyed by the record's
    ``Uploaded_variation`` and ``Location``, its allele, its feature and its feature kind,
    and compared column by column against the default format."""

    def __init__(self, default_text: str, tab_text: str, vcf_text: str, json_text: str,
                 parquet_text: str | None = None, inputs: list[str] | None = None,
                 examples: int = EXAMPLES_PER_GROUP) -> None:
        self.examples = examples
        objects = json_objects(json_text)
        if inputs is None:
            inputs = [obj["input"] for _, obj in objects if isinstance(obj.get("input"), str)]
        index = InputIndex(inputs, vcf_line_key)
        self.records: dict[str, dict] = {}
        self.objects_by_name: dict[tuple[str, str], list[dict]] = defaultdict(list)
        self.parts: dict[str, tuple[str, str, str]] = {}
        for n, obj in objects:
            _, label = json_record_label(index, n, obj)
            self.records[label] = obj
            self.objects_by_name[(uploaded_variation(obj), record_location(obj))].append(obj)
        self.default = self._rekey_rows(parse_default(default_text))
        self.tab = self._rekey_rows(parse_tab(tab_text))
        self.parquet = self._rekey_rows(parse_tab(parquet_text)) if parquet_text is not None else None
        vcf_entries = parse_vcf(vcf_text, inputs)
        self.vcf_records = [e for e in vcf_entries if e.key.endswith("|record")]
        self.vcf = self._rekey_by_record(vcf_entries)
        self.json = self._rekey_by_record(parse_json(json_text, inputs))
        self.tab_columns = header_columns(tab_text)
        self.parquet_columns = header_columns(parquet_text) if parquet_text is not None else []
        self.csq_columns = self._csq_columns(vcf_text)
        self.default_columns: set[str] = set()
        for e in self.default:
            self.default_columns.update(e.fields)

    @staticmethod
    def _csq_columns(text: str) -> list[str]:
        for line in text_lines(text):
            if line.startswith("##INFO=<ID=CSQ"):
                return csq_format_fields(line)
        return []

    def _rekeyed(self, keyed: list[tuple[Entry, str, str, str]]) -> list[Entry]:
        """Entries keyed by (name, Location, allele, feature, kind), duplicates numbered,
        each key's parts kept in ``self.parts`` (a name may itself carry ``|``)."""
        out = [
            Entry(
                key=f"{name}|{location}|{e.allele}|{dash(e.feature)}|{kind}",
                record=e.record,
                location=e.location,
                allele=e.allele,
                feature=e.feature,
                consequence=e.consequence,
                fields=e.fields,
            )
            for e, name, location, kind in keyed
        ]
        dedup_keys(out)
        for entry, (_, name, location, kind) in zip(out, keyed):
            self.parts[entry.key] = (name, location, kind)
        return out

    def _rekey_rows(self, rows: list[Entry]) -> list[Entry]:
        """Default-format rows keyed by their own ``Uploaded_variation`` and ``Location``."""
        return self._rekeyed(
            [
                (
                    row,
                    row.fields.get("Uploaded_variation", ""),
                    row.fields.get("Location", ""),
                    feature_kind(row.fields.get("Feature_type", "")),
                )
                for row in rows
            ]
        )

    def _rekey_by_record(self, entries: list[Entry]) -> list[Entry]:
        """VCF and JSON consequence entries keyed by the name and Location of the JSON
        object of their record; an entry whose record has no object keeps its label."""
        keyed = []
        for e in entries:
            if e.key.endswith("|record"):
                continue
            label, _, _, kind = e.key.split("|", 3)
            kind = feature_kind(kind.partition("#")[0])
            obj = self.records.get(label)
            if obj is None:
                keyed.append((e, label, "", kind))
            else:
                keyed.append((e, uploaded_variation(obj), record_location(obj), kind))
        return self._rekeyed(keyed)

    @staticmethod
    def _by_key(entries: list[Entry]) -> dict[str, Entry]:
        return {e.key: e for e in entries}

    def _pair(self, name: str, pivot: dict[str, Entry], other: dict[str, Entry], compare_fields) -> dict:
        pivot_only = sorted(set(pivot) - set(other))
        other_only = sorted(set(other) - set(pivot))
        shared = sorted(set(pivot) & set(other))
        differing: dict[str, list[str]] = defaultdict(list)
        totals: Counter[str] = Counter()
        compared_columns: set[str] = set()
        for key in shared:
            for column, detail in compare_fields(pivot[key], other[key], compared_columns):
                totals[column] += 1
                if len(differing[column]) < self.examples:
                    differing[column].append(f"{key}: {detail}")
        reasons = []
        if pivot_only:
            reasons.append(f"{name}: {len(pivot_only)} entries of the default format are absent")
        if other_only:
            reasons.append(f"{name}: {len(other_only)} entries have no default-format row")
        for column in sorted(totals):
            reasons.append(f"{name}: column {column} differs on {totals[column]} entries")
        return {
            "shared_entries": len(shared),
            "one_sided": {"default_only": len(pivot_only), f"{name}_only": len(other_only)},
            "default_only_examples": pivot_only[: self.examples],
            f"{name}_only_examples": other_only[: self.examples],
            "columns_compared": sorted(compared_columns),
            "differing": {c: {"count": totals[c], "examples": ex} for c, ex in sorted(differing.items())},
            "reasons": reasons,
        }

    def _compare_tab_like(self, columns: list[str]):
        declared = [c for c in columns if c != "Extra"]

        def compare(p: Entry, o: Entry, compared: set[str]):
            for column in declared:
                compared.add(column)
                raw = p.fields.get(column, ABSENT)
                pv = raw.replace("%3B", ";")
                ov = o.fields.get(column, ABSENT)
                if pv == ABSENT and ov == "-":
                    continue
                if pv != ov:
                    yield column, f"default {raw!r}, other {ov!r}"

        return compare

    def _compare_vcf(self):
        columns = [c for c in self.csq_columns if c not in VCF_IMPLIED_COLUMNS or c in self.default_columns]

        def compare(p: Entry, o: Entry, compared: set[str]):
            strand = self._object_for(p).get("strand")
            for column in columns:
                compared.add(column)
                pv = p.fields.get(column, ABSENT)
                ov = o.fields.get(column, "")
                if column == "Allele":
                    expected = pv if pv != ABSENT else "-"
                    if strand == -1:
                        expected = reverse_complement(expected)
                    expected = vcf_escape(expected.replace("%3B", ";"))
                elif pv in (ABSENT, "-"):
                    expected = ""
                else:
                    expected = vcf_escape(pv.replace("%3B", ";"))
                if expected != ov:
                    yield column, f"default {pv!r} (as CSQ {expected!r}), vcf {ov!r}"

        return compare

    def _compare_json(self):
        def scalar_differs(want: str | None, got: str | None) -> bool:
            if want is None or got is None:
                return (want is None) != (got is None)
            return not values_equal(want, got)

        def compare(p: Entry, o: Entry, compared: set[str]):
            kind = self.parts[p.key][2]
            remaining = set(o.fields)
            for column, raw in p.fields.items():
                if column in ("Extra", "Feature_type") or column in RECORD_LEVEL_COLUMNS:
                    continue
                if RECORD_LEVEL_COLUMN_RE.search(column):
                    continue
                compared.add(column)
                key = column.lower()
                pm = POSITION_RE.match(column)
                if pm:
                    coord = pm.group(1).lower()
                    start_key, end_key = f"{coord}_start", f"{coord}_end"
                    remaining.discard(start_key)
                    remaining.discard(end_key)
                    want_start = want_end = None
                    if raw not in ("-", ""):
                        parts = raw.split("-")
                        start_text = parts[0]
                        end_text = parts[1] if len(parts) > 1 and parts[1] != "" else None
                        want_start = start_text if POSITION_PART_RE.match(start_text) else None
                        if end_text is not None and (end_text.isdigit() or "?" in end_text):
                            want_end = end_text
                        else:
                            want_end = start_text
                        if end_text is not None and not POSITION_PART_RE.match(want_end):
                            want_end = None
                    for jkey, want in ((start_key, want_start), (end_key, want_end)):
                        got = o.fields.get(jkey)
                        if scalar_differs(want, got):
                            yield column, f"default {raw!r}, json {jkey}={got!r}"
                    continue
                if key in ("sift", "polyphen"):
                    pred_key, score_key = f"{key}_prediction", f"{key}_score"
                    remaining.discard(pred_key)
                    remaining.discard(score_key)
                    want_pred = want_score = None
                    if raw not in ("-", ""):
                        tm = TOOL_SCORE_RE.search(raw)
                        if tm:
                            want_pred, want_score = tm.group(1) or None, tm.group(2)
                    for jkey, want in ((pred_key, want_pred), (score_key, want_score)):
                        got = o.fields.get(jkey)
                        if scalar_differs(want, got):
                            yield column, f"default {raw!r}, json {jkey}={got!r}"
                    continue
                if key == "domains":
                    remaining.discard("domains")
                    want = None
                    if raw not in ("-", ""):
                        # OutputFactory/JSON.pm:363-371 re-splits each entry with ``(\w+):(\w+)`` and
                        # pushes ``{db => $1, name => $2}`` if ``$1 && $2``. A failed match leaves ``$1``
                        # and ``$2`` at the list's last successful one, so such an entry repeats the
                        # previous object and is dropped only before the first match. A capture of ``0``
                        # is false to Perl, so its entry and every unmatched one after it push nothing.
                        # ``numberify`` then runs over the record (JSON.pm:288), so a name of digits
                        # alone is a number.
                        want = []
                        db = name = None
                        perl_false = (None, "", "0")
                        for item in decode_extra(raw).split(","):
                            dm = DOMAIN_RE.search(item)
                            if dm:
                                db, name = dm.group(1), dm.group(2)
                            if db not in perl_false and name not in perl_false:
                                domain = {"db": perl_numeric_or_text(db), "name": perl_numeric_or_text(name)}
                                want.append(json_compact(domain))
                    # json_scalar joins the list's compact objects with ",", which each object
                    # also contains, so the comparison is on the joined text.
                    got = o.fields.get("domains")
                    if (want is None) != (got is None) or (want is not None and ",".join(want) != got):
                        yield column, f"default {raw!r}, json domains={got!r}"
                    continue
                if column in JSON_OBJECT_FIELDS:
                    remaining.discard(key)
                    want = None if raw in ("-", "") else numberified_json_text(raw)
                    got = o.fields.get(key)
                    if want != got:
                        yield column, f"default {raw!r}, json {key}={got!r}"
                    continue
                jkey = f"{kind}_id" if key == "feature" else JSON_RENAME.get(key, key)
                remaining.discard(jkey)
                want = default_list(raw)
                if want is not None:
                    want = ["1" if x == "YES" else x for x in want]
                if key in ("allele", "used_ref", "given_ref") and want is None:
                    want = ["-"]
                got = o.fields.get(jkey)
                got_list = got.split(",") if got is not None else None
                if (want is None) != (got_list is None) or (want is not None and not lists_equal(want, got_list)):
                    yield column, f"default {raw!r}, json {jkey}={got!r}"
            for jkey in sorted(remaining):
                compared.add(jkey)
                yield jkey, f"json carries {jkey}={o.fields[jkey]!r} and the default format no column for it"

        return compare

    def _object_for(self, e: Entry) -> dict:
        """The JSON object whose name and Location an entry's key carries, else empty."""
        name, location, _ = self.parts[e.key]
        return self.objects_by_name.get((name, location), [{}])[0]

    def _record_checks(self) -> dict:
        """Record-level values the JSON format carries once: ``most_severe_consequence``
        against the rank order of the object's own terms (``JSON.pm:387-388``),
        ``colocated_variants`` ids against the rows' ``Existing_variation``, the ``zyg``
        value against the rows' ``ZYG`` (the list of ``--individual_zyg`` as a multiset, since
        its order carries no meaning; the one zygosity of ``--individual`` as itself), and the
        VCF ``ID`` column against the record's name."""
        problems: list[str] = []
        existing_by_name: dict[tuple[str, str], set[str]] = defaultdict(set)
        zyg_by_name: dict[tuple[str, str], set[str]] = defaultdict(set)
        for row in self.default:
            name, location, _ = self.parts[row.key]
            existing_by_name[(name, location)].add(row.fields.get("Existing_variation", "-"))
            zyg = row.fields.get("ZYG", "-")
            if zyg != "-":
                zyg_by_name[(name, location)].add(",".join(sorted(zyg.split(","))))
        for (name, location), objects in sorted(self.objects_by_name.items()):
            want_ids: set[str] = set()
            for value in existing_by_name.get((name, location), set()):
                want_ids.add(",".join(sorted(default_list(value) or [])))
            got_ids: set[str] = set()
            got_zyg: set[str] = set()
            for obj in objects:
                terms = set()
                for list_name, _ in JSON_CONSEQUENCE_LISTS:
                    for item in obj.get(list_name) or []:
                        if isinstance(item, dict):
                            terms.update(json_list(item.get("consequence_terms")) or [])
                # An object with no consequence has no row to hold its zyg to: the
                # reference-homozygous copy of --individual (OutputFactory.pm:511 takes the
                # alternate alleles alone without --process_ref_homs) and a record carried
                # with no allele to annotate.
                if "zyg" in obj and terms:
                    got_zyg.add(",".join(sorted(json_scalar(obj["zyg"]).split(","))))
                want = most_severe(terms)
                got = obj.get("most_severe_consequence")
                if got is None and not terms and "id" not in obj:
                    continue
                if got != want:
                    problems.append(
                        f"{name}|{location}: most_severe_consequence {got!r}, the object's terms give {want!r}"
                    )
                colocated = obj.get("colocated_variants")
                ids: list[str] = []
                if isinstance(colocated, list):
                    ids = [json_scalar(c.get("id", "")) for c in colocated if isinstance(c, dict)]
                got_ids.add(",".join(sorted(ids)))
            if (name, location) in existing_by_name and want_ids != got_ids:
                problems.append(
                    f"{name}|{location}: Existing_variation {sorted(want_ids)} "
                    f"against colocated_variants ids {sorted(got_ids)}"
                )
            want_zyg = zyg_by_name.get((name, location), set())
            if (name, location) in existing_by_name and want_zyg != got_zyg:
                problems.append(f"{name}|{location}: ZYG {sorted(want_zyg)} against the objects' zyg {sorted(got_zyg)}")
        for e in self.vcf_records:
            label = e.key.split("|", 1)[0]
            obj = self.records.get(label)
            vcf_id = e.fields.get("ID", ".")
            if obj is None or vcf_id == "." or "id" not in obj:
                continue
            name = uploaded_variation(obj)
            if name != vcf_id:
                problems.append(f"{label}: VCF ID {vcf_id!r} against the name {name!r}")
        return {"problems": problems}

    def run(self) -> dict:
        pivot = self._by_key(self.default)
        pairs = {
            "tab": self._pair("tab", pivot, self._by_key(self.tab), self._compare_tab_like(self.tab_columns)),
            "vcf": self._pair("vcf", pivot, self._by_key(self.vcf), self._compare_vcf()),
            "json": self._pair("json", pivot, self._by_key(self.json), self._compare_json()),
        }
        if self.parquet is not None:
            pairs["parquet"] = self._pair(
                "parquet", pivot, self._by_key(self.parquet), self._compare_tab_like(self.parquet_columns)
            )
        records = self._record_checks()
        reasons = [r for pair in pairs.values() for r in pair["reasons"]] + records["problems"]
        keyed_or_record_checked = {"Uploaded_variation", "Location", "Existing_variation", "ZYG"}
        not_compared = {
            "default_record_level": sorted(
                c for c in self.default_columns
                if (c in RECORD_LEVEL_COLUMNS and c not in keyed_or_record_checked)
                or RECORD_LEVEL_COLUMN_RE.search(c)
            ),
            "vcf_implied": sorted(
                c for c in self.csq_columns if c in VCF_IMPLIED_COLUMNS and c not in self.default_columns
            ),
            "tab_undeclared": sorted(c for c in self.default_columns if c != "Extra" and c not in self.tab_columns),
            "vcf_undeclared": sorted(
                c for c in self.default_columns
                if c != "Extra" and c not in self.csq_columns and c not in RECORD_LEVEL_COLUMNS
            ),
        }
        return {
            "mode": "self-parity",
            "verdict": "red" if reasons else "green",
            "reasons": reasons,
            "records": len(self.records),
            "default_entries": len(pivot),
            "pairs": pairs,
            "record_checks": records,
            "not_compared": not_compared,
        }


WARNING_LINE_RE = re.compile(r"^WARNING:\s*(?:line (\d+) skipped|variant skipped)\s*(?:\((.*?)\): )?(.*)$")
WARNING_OTHER_RE = re.compile(r"^WARNING:\s*(.*?)(?:\s+(?:on )?line (\d+))?$")
SKIPPED_LINE_RE = re.compile(r"^Line\s+(\S+)\s*\t(.*)\t(.*)$")
PATH_RE = re.compile(r"(?<!\S)(?:~|\.{1,2})?/[^\s]+")


def mask_description(text: str) -> str:
    """A warning description with every path token replaced by ``<path>``."""
    return PATH_RE.sub("<path>", text).strip()


def parse_warning_lines(text: str) -> set[tuple[int | None, str]]:
    """The (input line number, description) set of a ``_warnings.txt`` or a
    ``--skipped_variants_file`` output; ``BaseVEP.pm:802-863`` and
    ``BaseRunner.pm:364-380`` are the two line shapes."""
    out: set[tuple[int | None, str]] = set()
    for line in text_lines(text):
        if not line.strip() or line.startswith("[VEP skipped the following variants from"):
            continue
        m = SKIPPED_LINE_RE.match(line)
        if m:
            number = int(m.group(1)) if m.group(1).isdigit() else None
            out.add((number, mask_description(m.group(3))))
            continue
        m = WARNING_LINE_RE.match(line)
        if m and (m.group(1) is not None or line.startswith("WARNING: variant skipped")):
            number = int(m.group(1)) if m.group(1) is not None else None
            out.add((number, mask_description(m.group(3))))
            continue
        m = WARNING_OTHER_RE.match(line)
        if m:
            number = int(m.group(2)) if m.group(2) is not None else None
            out.add((number, mask_description(m.group(1))))
            continue
        out.add((None, mask_description(line)))
    return out


# The Perl run-time warnings a plugin's undefined allele string raises before the plugin fails
# (Utils/Sequence.pm:1193 and :1252 at release/116), written to the warnings file as they are.
UNDEFINED_ALLELE_WARNING = re.compile(r"Use of uninitialized value \$(?:orig_a_alt|alt) in .*Sequence\.pm line \d+")
EXCEPTION_RULE = re.compile(r"^-{10,}(?: EXCEPTION -{10,})?$")


def without_plugin_failure(text: str, failure: str, undefined_allele_lines: bool = False) -> tuple[str, int]:
    """The reference's warnings text without the lines reporting one plugin's failure: each
    ``WARNING: Plugin '<name>' went wrong`` line (``failure`` is its text) with the exception
    block Ensembl VEP writes after it, the dashed ``EXCEPTION`` header through the closing
    dashed line (``BaseVEP.pm`` writes the exception's own text; a block the next ``WARNING``
    line interrupts ends there), and with ``undefined_allele_lines`` the Perl warnings about
    the undefined allele a CADD failure raises, wherever they stand. Returns the text and the
    number of lines dropped, blocks aside."""
    kept: list[str] = []
    dropped = 0
    lines = text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if failure in line:
            dropped += 1
            # An exception block follows a failure whose text the plugin threw: its dashed
            # header, then its lines up to and including the closing dashed rule.
            if i < len(lines) and "EXCEPTION" in lines[i] and EXCEPTION_RULE.match(lines[i].strip()):
                i += 1
                while i < len(lines) and not lines[i].startswith("WARNING"):
                    closing = EXCEPTION_RULE.match(lines[i].strip()) and "EXCEPTION" not in lines[i]
                    i += 1
                    if closing:
                        break
            continue
        if undefined_allele_lines and UNDEFINED_ALLELE_WARNING.search(line):
            dropped += 1
            continue
        kept.append(line)
    return "\n".join(kept), dropped


def without_forked_child_prefix(text: str) -> tuple[str, int]:
    """The reference's warnings text with the ``<pid> : `` a forked parent wrote before a
    child's first line removed from every line carrying one, and the number of lines changed.
    A child's own warning then reads as the serial run writes it; a child's interpreter
    warning, bare after the prefix, stays a line of its own."""
    rewritten: list[str] = []
    changed = 0
    for line in text.split("\n"):
        stripped = FORKED_CHILD_RELAY_RE.sub("", line)
        if stripped != line:
            changed += 1
        rewritten.append(stripped)
    return "\n".join(rewritten), changed


def compare_warning_files(expected_path: Path, actual_path: Path, masks: Collection[str] = ()) -> dict:
    """The two files as one-sided entry sets. Under ``regulatory_feature_rows`` the candidate's
    notice that it annotates no regulatory feature is required and left out of the sets; under
    ``fork_warning_pid_prefix`` the ``<pid> : `` before a relayed line of the reference is
    removed; under a class of ``PLUGIN_WARNING_CLASSES`` the reference's lines naming that
    plugin's failure are left out; each line changed or left out is counted under its class in
    ``masks_applied``."""
    expected_text = read_text(expected_path)
    applied: dict[str, int] = {}
    if FORK_WARNING_PID_PREFIX in masks:
        expected_text, relayed = without_forked_child_prefix(expected_text)
        if relayed:
            applied[FORK_WARNING_PID_PREFIX] = relayed
    for klass, text in PLUGIN_WARNING_CLASSES.items():
        if klass in masks:
            expected_text, dropped = without_plugin_failure(
                expected_text, text, undefined_allele_lines=klass == CADD_STRUCTURAL_ALLELE_MATCH
            )
            if dropped:
                applied[klass] = dropped
    expected = parse_warning_lines(expected_text)
    actual = parse_warning_lines(read_text(actual_path))
    notice_missing = False
    if REGULATORY_FEATURE_ROWS in masks:
        notice = {
            entry for entry in actual
            if REGULATORY_NOTICE in entry[1] or MISSING_REGULATORY_NOTICE in entry[1]
        }
        notice_missing = not notice
        actual -= notice

    def rows(items):
        numbered_first = sorted(items, key=lambda x: (x[0] is None, x[0] or 0, x[1]))
        return [{"line": n, "description": d} for n, d in numbered_first]

    return {
        "expected": str(expected_path),
        "actual": str(actual_path),
        "expected_entries": len(expected),
        "actual_entries": len(actual),
        "shared": len(expected & actual),
        "expected_only": rows(expected - actual),
        "actual_only": rows(actual - expected),
        "regulatory_notice_missing": notice_missing,
        "masks_applied": applied,
    }


def compare_warnings(
    expected: Path, actual: Path, expected_skipped: Path | None, actual_skipped: Path | None, masks: Collection[str] = ()
) -> dict:
    """The warnings-mode report: the ``_warnings.txt`` pair and, when given, the
    skipped-variants pair, each as one-sided sets; ``masks`` applies to the warnings pair."""
    sections = {"warnings": compare_warning_files(expected, actual, masks)}
    if expected_skipped is not None or actual_skipped is not None:
        if expected_skipped is None or actual_skipped is None:
            raise ComparatorError("--expected-skipped and --actual-skipped go together")
        sections["skipped"] = compare_warning_files(expected_skipped, actual_skipped)
    reasons = []
    for name, section in sections.items():
        if section["expected_only"]:
            reasons.append(f"{name}: {len(section['expected_only'])} entries only in the expected file")
        if section["actual_only"]:
            reasons.append(f"{name}: {len(section['actual_only'])} entries only in the actual file")
        if section["regulatory_notice_missing"]:
            reasons.append(f"{name}: the regulatory notice is missing from the actual file")
    return {
        "mode": "warnings",
        "format": "warnings",
        "verdict": "red" if reasons else "green",
        "reasons": reasons,
        **sections,
    }


def parse_floors(spec: str | None) -> dict[str, int]:
    """``COL=MIN,...`` into a column-to-minimum map."""
    floors: dict[str, int] = {}
    if not spec:
        return floors
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        name, sep, value = item.partition("=")
        if not sep or not name or not value.isdigit():
            raise ComparatorError(f"--floors expects COL=MIN with a non-negative integer, got {item!r}")
        floors[name] = int(value)
    return floors


def parse_field_pairs(specs: list[str] | None) -> dict[str, set[tuple[str, str]]]:
    """``COL=EXPECTED=>ACTUAL`` items, several per argument separated by ``;``."""
    pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for spec in specs or []:
        for item in spec.split(";"):
            if not item.strip():
                continue
            name, sep, rest = item.partition("=")
            expected, arrow, actual = rest.partition("=>")
            if not sep or not name or not arrow:
                raise ComparatorError(f"--field-pairs expects COL=EXPECTED=>ACTUAL, got {item!r}")
            pairs[name].add((expected, actual))
    return dict(pairs)


def parse_masks(spec: str | None) -> set[str]:
    return {m.strip() for m in (spec or "").split(",") if m.strip()}


def load_manifest(path: Path | None) -> dict | None:
    if path is None:
        return None
    try:
        manifest = json.loads(read_text(path))
    except ValueError as exc:
        raise ComparatorError(f"{path} is not JSON: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ComparatorError(f"{path} is not a manifest object")
    return manifest


def write_report(report: dict, path: Path | None) -> None:
    if path is None:
        return
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    except OSError as exc:
        raise ComparatorError(f"cannot write {path}: {exc.strerror or exc}") from exc


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="compare_vep_formats.py",
        description="Per-column comparison of VEP output files; see the module docstring for the rules.",
        epilog=f"environment: {FRONT_END_VARIABLE}, when it names an executable (a path, or a command on PATH), "
        "hands a default or tab comparison over to that build of tools/vep-formats-remainder: the front end keeps "
        "back the record groups both outputs write alike and the rest is compared here, the report and the "
        "standard output unchanged; one line on the standard error says what was handed over. Unset or empty, "
        "the whole outputs are compared. A value naming no executable, and a front end that fails, exit 2.",
    )
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--self-parity", action="store_true", help="one engine's formats against each other")
    mode.add_argument("--warnings", action="store_true", help="two engines' warnings files against each other")
    ap.add_argument(
        "--format",
        choices=FORMATS,
        help="the format of --expected and --actual; recoder compares two variant_recoder outputs",
    )
    ap.add_argument("--expected", type=Path, help="the Ensembl VEP output, or its warnings file")
    ap.add_argument("--actual", type=Path, help="the candidate output, or its warnings file")
    ap.add_argument(
        "--input", type=Path,
        help="the input VCF (needed for --format vcf and json; beside --manifest its name selects the "
        "flag set documented over that input in every format, the text formats reading nothing else of it "
        "but under reference_skipped_record, the classes that read a record's alleles and "
        "splice_term_past_the_record_span with its plugin value class, which read from it the reference "
        "bases of a row under an identifier)",
    )
    ap.add_argument("--manifest", type=Path, help="the corpus manifest.json documenting divergences")
    ap.add_argument(
        "--run",
        help=f"the run of --manifest the pair belongs to: {OWN_RUN}, the manifest's own run, or the name of one of "
        "its flag sets; the lists of that run alone are read, and an --input of another name than the run's is "
        "refused; without it every run over the input is read as one list",
    )
    ap.add_argument(
        "--reference-sequence",
        choices=(REFERENCE_PRESENT, REFERENCE_ABSENT),
        help=f"whether the candidate's run had a reference sequence: under {REFERENCE_PRESENT} a divergences item "
        f"that carries {WITH_REFERENCE_SETS} documents those sets, and none where they are Ensembl VEP's own; "
        f"under {REFERENCE_ABSENT}, and without the option, an item documents its vep_rs_consequence_sets",
    )
    ap.add_argument(
        "--masks",
        help="comma-separated intended-divergence classes applied: the manifest's, and the built-in "
        "regulatory_feature_rows (--warnings mode takes it with the plugin warning classes below), "
        "order_only (DOMAINS, FREQS, ZYG and NearestExonJB compared as multisets, the "
        "5UTR_annotation groups by their pairs), individual_allele_order (an individual's "
        "alleles compared as a set), conservation_short_name_spelling (the reference's "
        "phlyoP<N>way key read as phyloP<N>way), pc_description_text (the reference's PC "
        "header line read with the candidate's description), motif_score_reverse_strand_mnv "
        "(MOTIF_SCORE_CHANGE on a reverse-strand motif entry whose allele is two or more bases), "
        "colocated_list_cut_at_cache_cap (a PUBMED or VAR_SYNONYMS list the reference writes with one last item "
        "more than the candidate, on a stored list of exactly 1,024 bytes that ends in no separator), "
        "the registry pair classes start_cooccurrence_swap and splice_family_swap (a key whose "
        "consequence sets differ by the shape the SNP/indel comparator's excluding rule admits, "
        "set aside whole), the field classes gxf_reverse_strand_utr_sequence (the fsTer or "
        "extTer count of an HGVSp on a reverse-strand transcript read from a GTF or GFF file), "
        "conservation_max_text_order (a Conservation MAX below the candidate's negative value), "
        "gwas_haplotype_cross_product (under --gwas-catalog, a GWAS_ field the reference alone fills "
        "or both fill differently, the empty fields of an annotated reference row included, or one "
        "the candidate alone fills where both rows are annotated under different studies, on a row "
        "at a location of a catalog row of two or more SNPs whose allele is another of its SNPs' risk "
        "alleles, or any of them where the row's SNPs and several locations differ in count, and "
        "whose reference GWAS_study names its study), "
        "gwas_interaction_first_snp_allele_whitespace (under --gwas-catalog, a GWAS_ field the "
        "candidate alone fills on a row at the location of an SNP written before x in a catalog row, "
        "its allele that SNP's risk allele and the candidate's GWAS_study naming the row's study, the "
        "empty fields of the candidate's annotated row included), "
        "alphamissense_two_digit_transcript_version (an am_ field the candidate "
        "alone fills), cadd_structural_allele_match (a CADD_ field the candidate alone fills on a "
        "structural allele's row; in --warnings mode the reference's Plugin 'CADD' went wrong "
        "lines), splicevault_sibling_allele_rows (a SpliceVault_ field the reference fills with the "
        "value its row of another allele of the record carries on the same feature, where the "
        "candidate leaves it empty or fills it differently), "
        "splicevault_json_later_rows (the JSON splicevault object the candidate alone fills), "
        "utrannotator_annotation_pair_order (a 5UTR_annotation whose groups carry the same pairs in "
        "different orders, the candidate's pairs in byte order of their keys), "
        "utrannotator_gained_start_window (a 5UTR_consequence one engine alone writes as the gained "
        "start, with the row's 5UTR_annotation: the candidate at a CapDistanceToStart inside an "
        "allele longer than the reference allele and past that allele's end, the reference inside a "
        "reference allele longer than the allele and past the allele's end, or the reference at a "
        "negative CapDistanceToStart; a sequence allele at a two-base Location under an "
        "identifier needs --input; an insertion's row is refused where the reference reports a start "
        "codon lost and the candidate's gained start stands the allele's length past it), "
        "utrannotator_consequence_precedence (under --utrannotator-sweep, a "
        "5UTR_consequence whose two values differ and are both terms those runs write on the row, "
        "with the row's 5UTR_annotation where each side's is one they write beside its term), "
        "utrannotator_structural_allele_text (the 5UTR_consequence and 5UTR_annotation the reference "
        "alone fills on a row whose Allele is a structural class word, its term one of the "
        "plugin's five and the candidate's row filling neither field; and those the candidate "
        "fills on a default or tab row of a symbolic deletion inside one exon of the 5' UTR, or "
        "of a plain symbolic insertion that carries the frameshift its stated length decides), "
        "utrannotator_lost_start_after_insertion (a 5UTR_consequence that is the lost-start term on "
        "the reference's row alone, with the row's 5UTR_annotation, where the annotation's one "
        "CapDistanceToStart is the lower cDNA_position of an insertion's row, or lies one or two below "
        "it on a row whose Allele, read on the transcript's strand, ends with the start codon's first "
        "base or its first two: the insertion only moves the start codon on, which stands whole), "
        "utrannotator_gained_start_remade_in_place (the two UTRAnnotator fields of an insertion's row on "
        "which the reference writes the gained start, itself or in a run of --utrannotator-sweep, at a "
        "CapDistanceToStart one or two below the lower cDNA_position, where the Allele, read on the "
        "transcript's strand, opens with the start codon's last two bases or its last and the candidate "
        "writes no gained start there: the inserted bases leave the reference's own start codon in place), "
        "utrannotator_frameshift_at_stop_boundary (the 5UTR_annotation of a row on which both write the "
        "frameshift and the candidate lists, beyond every group of the reference's, only uORFs whose "
        "alt_type_length differs from their ref_type_length and whose stop codon begins directly after the "
        "insertion point the row's HGVSc and HGVS_OFFSET give: the reference's test leaves out an upstream "
        "ORF an insertion written directly before its stop codon shifts; a row without an ins or dup "
        "notation between 5' UTR positions is refused), "
        "utrannotator_lost_stop_written_before_stop_codon (the two UTRAnnotator fields of an insertion's "
        "row on which the reference writes the frameshift and the candidate the lost stop, each group of "
        "the candidate's carrying the KozakContext and KozakStrength of exactly one uORF the reference "
        "lists and the alternative stop codon where that uORF's alt_type_length ends it: the inserted "
        "bases can be written inside the uORF's stop codon, which the reference reads only for a record "
        "written there), "
        "nearest_exon_jb_item_order (a NearestExonJB whose two values carry the same items in different orders), "
        "loeuf_zero_minimum (a LOEUF of 0 the candidate alone fills), maxentscan_overlapping_windows "
        "(on a row whose one-base record lies in the acceptor window of one intron of its transcript and the "
        "donor window of the next, read from --json-cache and nothing without it: a MaxEntScan_ field both "
        "fill with different numbers, and a window text of a verbose run where the reference writes nine "
        "bases and the candidate twenty-three), alphamissense_sibling_allele_rows "
        "(an am_ field the reference fills with the value its row of another allele of the record "
        "carries on the same feature beside the same Amino_acids), revel_sibling_allele_rows (a REVEL "
        "of the same shape), revel_other_reference_amino_acid_rows (a REVEL the reference fills with the "
        "score of the table row its walk reaches, a row whose reference amino acid is not the row's own, "
        "where the candidate writes nothing or the score of the first row of the row's allele that names "
        "both residues of Amino_acids; needs --revel-table, and --revel-no-match for a run that set the "
        "plugin's no_match), gnomadmt_info_value_second_equals (a gnomAD_ field both fill, the "
        "reference's value the candidate's up to its first =), loftee_anc_allele_multiallelic_record "
        "(on a one-base substitution's stop-gained or splice-site row of a multi-allelic record, a "
        "LoF_filter the candidate writes as the reference's plus ANC_ALLELE, with the row's LoF, LC "
        "beside the item and HC without it), loftee_first_alt_donor_filter (on a "
        "splice_donor_variant row of a multi-allelic record, a LoF_filter whose two values differ "
        "by GC_TO_GT_DONOR alone, either side carrying it, with the row's LoF read the same way, "
        "and not on the row of the record's first alternate; nothing where no entry states the "
        "record's allele order), "
        "loftee_gerp_zoom_level_summary (on a stop-gained or frameshift row whose BP_DIST is at "
        "least 79, a LoF_info whose two values differ in the GERP_DIST item alone, and, where "
        "DIST_FROM_LAST_EXON is at most 50 and the two totals straddle -58, a LoF_filter differing "
        "by END_TRUNC alone on the side of the lower total, with the row's LoF), "
        "nearest_exon_jb_stretch_priming (a NearestExonJB the candidate alone fills on a row whose "
        "DISTANCE is at most 12), nearest_exon_jb_shifted_span (under --shift_3prime 1 "
        "--shift_length, a NearestExonJB on a transcript row with a positive SHIFT_LENGTH that one "
        "side fills and the other leaves empty, or whose two values name different exons), "
        "nearest_exon_jb_id_less_exon_fold (a NearestExonJB both fill with items of an empty exon "
        "id, the reference's one item whose distance and side are those of one of the candidate's, "
        "on a row the candidate gives two or more items, or whose Consequence carries both "
        "splice-site terms and whose cDNA position names two positions), "
        "downstream_reverse_strand_rotation (under --hgvs, a DownstreamProtein or ProteinLengthChange "
        "both fill differently on a frameshift_variant row of a STRAND -1 transcript whose HGVS_OFFSET exceeds "
        "the number of bases its HGVSc inserts, duplicates or deletes and is not a multiple of it), "
        "splice_family_swap on the Downstream fields (a DownstreamProtein or ProteinLengthChange the "
        "reference alone fills on a key whose Consequence sets pair by that class's rule, the "
        "reference's row a frameshift row none of whose terms names a splice site), "
        "spliceai_chopped_multiallelic (a SpliceAI_ field the candidate alone fills on a substitution "
        "allele's row of a record whose rows carry an allele written -), "
        "spliceai_multi_entry_cell (a SpliceAI_ field on a row whose score record holds several comma-joined "
        "entries, the candidate's value the entry of the row's gene and the reference's the record's whole "
        "value on a row of the first entry's gene or none on a row of a later entry's; of record lines that "
        "stand more than once at an allele the reference's is the last the row's gene leads and the "
        "candidate's entry that of the last that names it; needs --spliceai-table, and --json-cache for a "
        "row that writes no SYMBOL), "
        "plugin_contig_mitochondrial_spelling (a field of a tabix-backed plugin other than gnomADMt, "
        "of those the plugins name themselves, that the candidate alone fills on a row of a record "
        "spelt MT, chrM or chrMT, and beside a --manifest whose assembly is GRCh37 an AlphaMissense "
        "field alone; the three sibling-allele classes, the two LoF_filter classes and "
        "spliceai_chopped_multiallelic read a record's alleles from the VCF and JSON entries and "
        "from a default or tab row whose Uploaded_variation spells them or, under --input, whose "
        "identifier names one input line: the one line that carries it, or of several the one "
        "that writes the row at its Location with its allele; they leave any other row under an "
        "identifier a difference), "
        "refseq_edited_model_hgvsc_offset "
        "(an HGVSc on a BAM_EDIT OK row differing in its coordinates alone), "
        "refseq_edited_model_shifted_insertion_bases (under --hgvs, two insertion notations of one "
        "transcript and coordinate pair on a BAM_EDIT OK insertion row both engines give one non-zero "
        "HGVS_OFFSET, the reference's bases the row's Allele on the transcript's strand by STRAND and "
        "the candidate's those bases turned 3' by the offset), "
        "refseq_edited_model_seam_insertion_hgvsc (under --hgvs, the HGVSc of an insertion standing on the "
        "bases an edited model holds and the genome lacks, on a row carrying BAM_EDIT OK or FAILED on both "
        "sides: the reference's an ins of the allele, the candidate's a notation that, applied to the row's "
        "reference codons, gives its alternate codons, or none where the allele is the model's bases), "
        "refseq_edited_model_unminimised_allele_hgvsc (under --hgvs, the HGVSc of a row carrying "
        "BAM_EDIT OK and one GIVEN_REF and one USED_REF on both sides: two substitutions, two ins or two "
        "delins of one sequence name, numbering and coordinate text, the reference's spelling the row's "
        "USED_REF and Allele whole on the transcript's strand and the candidate's spelling them without the "
        "bases the clip takes from GIVEN_REF and the Allele, which leave of GIVEN_REF the positions the "
        "coordinates cover), "
        "refseq_edited_model_hgvs_offset (an HGVS_OFFSET on a coding BAM_EDIT OK row whose HGVSc is a "
        "duplication, which the reference leaves empty and the candidate fills with a non-zero integer), "
        "hgvs_shifted_insertion_clamped_end (under --shift_3prime 1 --shift_length, the HGVSc of an "
        "insertion's row with a positive SHIFT_LENGTH, a cDNA position starting at it and an "
        "HGVS_OFFSET: the reference's at one position, a substitution, a delins or an insertion "
        "before it, the candidate's a duplication ending at the base before or an insertion at the "
        "same point, of as many bases as the allele), "
        "hgvs_intron_offset_abutting_insertion (the HGVSIntronOffset plugin's two fields the reference alone "
        "fills, 0 and 0, or 1 and 0 past a reverse-strand transcript's last base, or 0 and -1 before a "
        "forward-strand transcript's first base, on an insertion's row both engines annotate outside the "
        "transcript at a DISTANCE of 0, with no HGVSc and no HGVS_OFFSET), "
        "hgvs_shifted_indel_written_span_gate (under --hgvs, an HGVSc, HGVS_OFFSET or HGVSIntronOffset "
        "field the reference leaves empty on an insertion's or deletion's row annotated inside the "
        "transcript, where the candidate writes the shifted span's indel notation of the row's transcript "
        "beside a non-zero offset, the row's SHIFT_LENGTH positive where printed; beside a del "
        "notation the plugin's two fields are the intron offsets of its coordinates), "
        "hgvsp_shifted_indel_written_span_gate (under --hgvs, an HGVSp the reference leaves empty on an "
        "insertion's or deletion's row annotated in an intron, with no cDNA, CDS or protein position, "
        "where both engines write one shifted HGVSc placing the allele in the coding sequence beside one "
        "non-zero HGVS_OFFSET, in an exon whose first base is coding base 2 or later, and the candidate "
        "writes the protein notation), "
        "hgvs_shifted_deletion_clamped_span (under --hgvs, two del notations of one transcript and "
        "numbering ending at the same cDNA coordinate, the reference's starting later than the "
        "candidate's, on a deletion's row both give a non-zero HGVS_OFFSET, the candidate's span "
        "agreeing with the row's cDNA position where --shift_length prints the shifted one; on "
        "such a row an HGVS_IntronStartOffset each engine writes as the intron offset of its own "
        "notation's start coordinate), "
        "hgvsp_stop_codon_protein_comparison (under --hgvs, two HGVSp values naming one protein and "
        "the stop at one residue number: the reference's the stop replaced by a stop, by the first "
        "residues that follow the protein or by all of them and an inserted stop, residues inserted "
        "before the stop, or the stop unchanged on a stop_lost row; the candidate's the stop unchanged "
        "for the first and an extension named at the stop for the others; or the reference's a range "
        "from a residue to the stop replaced by residues on a stop_lost row, the candidate's a "
        "frameshift named at that residue with the first of them; or the reference's an extension "
        "of unknown count on a stop_lost row, the candidate's the same extension with a count; or "
        "the reference's the stop written as deleted with an extension count on a stop_lost row, the "
        "candidate's an extension with that count and the first new residue; or the reference's one "
        "residue and a stop inserted between two neighbouring residues, the residue the right flank's, "
        "on a stop_retained_variant row, the candidate's the stop unchanged at the residue after the "
        "right flank; the row's terms are the reference's; no "
        "row whose FLAGS name "
        "cds_end_NF), "
        "hgvsp_stop_named_on_incomplete_coding_sequence (under --hgvs, a reference HGVSp that names "
        "the stop alone, as unchanged, deleted, replaced by one residue or extended, where the "
        "candidate writes none, on a row whose FLAGS name cds_end_NF, whose terms hold "
        "incomplete_terminal_codon_variant, or whose description is the stop unchanged at residue 1 "
        "under FLAGS naming cds_start_NF), "
        "hgvs_clipped_reference_against_whole_alternate (under --hgvs, the HGVSc of a sequence "
        "allele's row whose cDNA position has an open side: the reference's a notation of the row's "
        "transcript, or empty; the candidate's a notation of the same transcript and numbering "
        "whose every coordinate lies within the row's span plus one of the reference's, empty "
        "against a reference ins or dup, empty against a reference del, delins or substitution on "
        "the row of an allele shorter than its span, or a del against the empty reference; with --json-cache, "
        "on a substitution's row whose cDNA position is open or empty, two delins notations whose "
        "bases are runs of the row's allele, the reference's holding the bases outside the "
        "transcript's exon spans that the clamp cut from the candidate's; and the HGVSIntronOffset "
        "plugin's two fields the reference alone fills beside an ins or dup it prints against an empty "
        "candidate, the candidate alone fills beside an ins of the bases a longer allele adds inside the "
        "span where the reference prints no notation, on a reverse-strand transcript's row whose cDNA position "
        "is open at its start and whose HGVSc stays a difference, or both fill "
        "beside two ins notations the class takes as a pair, each the intron part of its own side's "
        "notation's coordinate it is read from; and the plugin's end offset alone, on a forward-strand "
        "transcript's row whose cDNA position is open at its start, where the reference prints an ins "
        "from one base before the transcript's first base of the allele's bases less those inside, "
        "beside the pair 0 and -1, and the candidate the duplication of the bases the allele adds, "
        "beside the pair 0 and 0 and the bases it was carried in HGVS_OFFSET, the duplication ending "
        "as many bases into the transcript as the span's bases inside and that HGVS_OFFSET together; "
        "the HGVS_OFFSET of such a row stays a difference), "
        "hgvs_added_bases_left_at_the_span_end (under --hgvs, on a row whose allele holds more bases "
        "than its Location spans and whose reference prints an HGVSc and no HGVS_OFFSET: a reference "
        "ins one base past a coordinate beside the HGVSIntronOffset plugin's pair 1 and 0, against a "
        "candidate duplication of as many bases ending at that coordinate, or against a candidate "
        "that writes nothing where the added bases have no placement inside the transcript; and a "
        "reference ins between two exonic neighbours, or dup, against a candidate dup or ins the "
        "bases its HGVS_OFFSET states further 3', or against the dup of the bases that end at that "
        "ins beside no HGVS_OFFSET), "
        "hgvs_insertion_left_where_the_shift_stops (under --hgvs, on an insertion's row whose "
        "reference HGVS_OFFSET is in size the most steps the reference's shift walks an allele of "
        "that length, or is absent where that is none, beside an ins or dup: a candidate dup or ins beside a larger HGVS_OFFSET of the "
        "same sign, or a candidate that writes no HGVSc where the most 3' placement lies past the "
        "transcript's last base; with them the HGVSIntronOffset plugin's two fields, each row's the "
        "intron parts of its own notation), "
        "hgvs_added_bases_written_outside_the_transcript (under --hgvs, on an insertion's row both "
        "engines annotate upstream or downstream of the transcript: a candidate dup or ins of the "
        "row's transcript with the HGVSIntronOffset plugin's two fields of its coordinates, beside "
        "no HGVS_OFFSET on a downstream row and a non-zero one on an upstream row, larger than the "
        "DISTANCE on an insertion's row and no smaller on the row of an allele longer than its "
        "span, against a reference that prints no HGVSc and no HGVS_OFFSET and writes "
        "another pair or none), "
        "hgvsp_intron_end_insertion (under --hgvs, an HGVSp the reference fills and the candidate "
        "leaves empty on a row whose HGVSc, one notation of the row's transcript on both sides beside "
        "no HGVS_OFFSET, duplicates intron positions ending on the one before an exon's first coding "
        "base, or inserts bases between that position and the base, whether or not they end in a "
        "pair that ends an intron, the record written where an intron of more than 13 bases of that "
        "transcript ends before a coding exon base; a row with an HGVS_OFFSET and the end of an "
        "intron of 13 bases or fewer stay differences; a record that holds the insertion with a "
        "reference base it shares is placed by the gap its allele, clipped of that base to the "
        "bases the HGVSc inserts, lies in; a default or tab row that gives such a record's Location "
        "is placed by the gaps of that span, its allele read on either strand or, minimised, "
        "equal to the inserted bases or as long as the duplicated positions; a VCF or JSON entry is "
        "placed from its line of --input, in VCF or in Ensembl format, by an ALT that clips against "
        "the REF to the insertion the HGVSc describes, however the entry spells its allele; needs "
        "--json-cache), "
        "hgvsp_intron_start_insertion (under --hgvs, an HGVSp the reference fills and the candidate "
        "leaves empty on a row whose HGVSc, one notation of the row's transcript on both sides, "
        "inserts bases between an exon's last coding base and the first position of the intron "
        "after it, beside no HGVS_OFFSET or one non-zero integer on both sides, the record's gap, "
        "moved by that offset or as the Location of a run under --shift_genomic 1 names it, one "
        "where an exon of that transcript ends inside the coding region before an intron of any "
        "length; a record that holds the insertion with a reference base it shares carries no "
        "offset and is placed by the gap its clipped allele lies in; a default or tab row that "
        "gives such a record's Location is placed by the gaps of that span, its allele read on "
        "either strand or, minimised, equal to the inserted bases; on a row that carries the offset "
        "the allele is read as the record writes it or as a run under --shift_3prime 1 prints it, "
        "turned by the offset; a VCF or JSON entry is placed from its line of --input as under "
        "hgvsp_intron_end_insertion; needs --json-cache), "
        "shift_intron_term_unshifted_span (under --shift_3prime 1 --shift_length, a Consequence set "
        "the reference writes as the candidate's plus intron_variant, on a row with a positive "
        "SHIFT_LENGTH whose HGVSc places the shifted insertion or deletion at the acceptor's two "
        "bases or in the exon), "
        "refseq_edited_model_deleted_span_term (the Consequence of a BAM_EDIT OK transcript row of a "
        "sequence allele inside the bases the model's alignment deleted: intergenic_variant in the "
        "reference, exactly one of coding_sequence_variant, 5_prime_UTR_variant and 3_prime_UTR_variant "
        "in the candidate, no cDNA position and one IMPACT on both; the JSON record's "
        "most_severe_consequence follows), "
        "refseq_edited_model_tail_term (the Consequence of a BAM_EDIT OK transcript row of a sequence "
        "allele on one of the final bases of a model whose edits inserted bases: intergenic_variant in "
        "the reference, 3_prime_UTR_variant in the candidate, a cDNA position and one IMPACT on both; "
        "the JSON record's most_severe_consequence follows), "
        "mature_mirna_unshifted_span (under --shift_3prime 1 --shift_length, the Consequence of a "
        "transcript row with a positive SHIFT_LENGTH whose two sets are equal but for one engine's "
        "mature_miRNA_variant against the other's non_coding_transcript_exon_variant; the JSON "
        "record's most_severe_consequence follows), "
        "stop_codon_unshifted_span (under --shift_3prime 1 --shift_length, the Consequence of a "
        "transcript row with a positive SHIFT_LENGTH whose two sets are equal but for the reference's "
        "coding_sequence_variant against the candidate's stop_retained_variant; the JSON record's "
        "most_severe_consequence follows), "
        "stop_retained_without_stop_codon_test (the Consequence of an insertion's or deletion's "
        "transcript row whose two sets are equal but for the reference's stop_retained_variant, against "
        "nothing, coding_sequence_variant or frameshift_variant, both rows showing that the stop codon "
        "test is not made: an allele letter that names no base, a cDNA position missing at one end, or "
        "no CDS start; "
        "the cells term_pair_following_cell reads from the two sets of terms and the JSON record's "
        "most_severe_consequence follow), "
        "stop_codon_overlap_genomic_span (the Consequence of an insertion's or deletion's transcript "
        "row with both cDNA positions and a CDS start whose two sets are equal but for the reference's "
        "coding_sequence_variant against the candidate's stop_retained_variant, both holding "
        "splice_region_variant; "
        "the cells term_pair_following_cell reads from the two sets of terms and the JSON record's "
        "most_severe_consequence follow), "
        "stop_retained_three_base_stop_codon_test (the Consequence of an insertion's or deletion's "
        "transcript row whose stop codon test is made, both cDNA positions, a CDS start and an allele of "
        "bases on both sides: the reference's stop_retained_variant against coding_sequence_variant, or "
        "it, alone or beside an in-frame or protein-altering term, against frameshift_variant with or "
        "without stop_gained; "
        "the cells term_pair_following_cell reads from the two sets of terms and the JSON record's "
        "most_severe_consequence follow), "
        "stop_codon_removed_without_stop_codon_test (the Consequence of an insertion's or deletion's "
        "transcript row without a cDNA position at one end or a CDS start on both sides, whose two sets "
        "are equal but for the reference's stop_retained_variant against the candidate's stop_lost; "
        "the cells term_pair_following_cell reads from the two sets of terms and the JSON record's "
        "most_severe_consequence follow), "
        "stop_codon_overlap_insertion_range (the Consequence of an insertion's transcript row with both "
        "cDNA positions and a CDS start on both sides, whose two sets are equal but for the reference's "
        "coding_sequence_variant against the candidate's stop_retained_variant away from "
        "splice_region_variant, or the reference's frameshift_variant with or without stop_gained against "
        "the candidate's stop_retained_variant with or without inframe_insertion; "
        "the cells term_pair_following_cell reads from the two sets of terms and the JSON record's "
        "most_severe_consequence follow), "
        "splice_term_past_the_record_span (On a record whose alternate allele is longer than its "
        "reference allele, Ensembl VEP writes a splice donor, splice acceptor, donor fifth base, donor "
        "region, polypyrimidine tract or intron term from positions past the end of the record's span; "
        "where no base of the span lies in a window of that term, vep-rs does not write it. The class takes "
        "the Consequence of a transcript row whose cDNA and CDS positions are the same on both sides and "
        "whose two sets are equal but for the reference's terms among splice_donor_variant, "
        "splice_acceptor_variant, splice_donor_5th_base_variant, splice_donor_region_variant, "
        "splice_polypyrimidine_tract_variant and intron_variant, whose window is the intron, and the "
        "candidate's splice_donor_region_variant or splice_region_variant where a term that withholds it "
        "went. With the record's alleles (a default or tab row whose Uploaded_variation spells them, a VCF or "
        "JSON entry's own, or the records of --input for a row under an identifier) the record is no pure "
        "insertion, the reference bases are two or more and as many as the Location's positions, the allele "
        "is longer, and one intron stands beside the span, at a place the row's cDNA positions allow, at "
        "which the record's differing regions give the reference's terms and the span gives the candidate's. "
        "Without them the Location holds three or more positions, the allele more nucleotide letters than the "
        "span has bases and the two cDNA positions are numbers as far apart as the span's two ends, the "
        "candidate holding none of the six, with or without its splice_region_variant where one of the first "
        "four went. A row that differs by the terms of a stop codon class as well is taken where that "
        "class's rule takes the rest of the difference; "
        "the cells term_pair_following_cell reads from the two sets of terms and the JSON record's "
        "most_severe_consequence follow), "
        "refseq_edited_model_deleted_span_codons (the Codons of a BAM_EDIT OK transcript row whose "
        "USED_REF has more bases than its CDS span: the alternate codons equal, the reference codons "
        "of one length differing in an upper-case base alone, the reference's upper-case bases "
        "USED_REF's leading bases on the transcript's strand), "
        "frequency_filter_zero_value (a FREQS value the reference reads NA and the candidate 0, and "
        "the rows of a record the candidate alone keeps on a FREQS of zeros), "
        "merged_same_stable_id_collapse (a transcript row the candidate alone writes under --refseq "
        "or --merged, of a model the reference wrote no row for on a record it wrote a transcript row "
        "for, its SOURCE RefSeq or absent: Ensembl VEP keeps one transcript per stable identifier "
        "among the models a buffer loads), "
        "hgnc_id_copy (an HGNC_ID one side leaves empty and the other fills with an id, on a transcript "
        "row whose SYMBOL_SOURCE is not HGNC: Ensembl VEP lends a same-symbol gene's id per buffer), "
        "organelle_peptide_paths_translation_table (the HGVSp of a frameshift, of a lost stop's "
        "extension or of a duplication against an insertion, both values naming one protein, on a "
        "transcript the candidate JSON cache files under a codon table that reads some codon unlike "
        "the standard one, the row's consequence set, positions, amino acids, codons, HGVSc and "
        "HGVS_OFFSET equal on both sides; "
        "reads --json-cache, and without it applies from a manifest alone), "
        "the Downstream plugin's classes downstream_utr_frame_past_coding_end (the DownstreamProtein of "
        "a frameshift deletion's row that starts on a codon's third base and ends on the last base of "
        "the translateable sequence: the reference's the 3' UTR translated from its first base, the "
        "candidate's the two coding bases before the deletion translated with the UTR) and "
        "downstream_start_past_own_codon (the DownstreamProtein of any other frameshift deletion's "
        "row, with or without inserted bases, that starts on a codon's third base: the reference's "
        "translated from the codon after the first one the change touches, the candidate's from that "
        "codon), each on a row whose consequence set (or two sets that differ as start_cooccurrence_swap "
        "admits), CDS_position (or a later allele's two, as plugin_hgvs_shifted_cds_coordinates reads them) and ProteinLengthChange are equal on "
        "both sides, both values recomputed from the candidate cache's sequences under the transcript's "
        "codon table on the record's own CDS coordinates, and plugin_hgvs_shifted_cds_coordinates (a "
        "DownstreamProtein, a ProteinLengthChange or an NMD on a row whose reference carries a non-zero "
        "HGVS_OFFSET: under --hgvs Ensembl VEP's plugin reads the CDS coordinates and the allele the "
        "protein notation left shifted, the candidate the record's own; the reference's value is "
        "recomputed on the shifted coordinates, the candidate's on the own, and Ensembl VEP's rule gives "
        "another value on the own; on the row of a later alternate allele of its record whose reference "
        "prints the shifted coordinates as its CDS_position the own coordinates are the candidate's; "
        "a row the start rule changes as well needs the Downstream class of its "
        "own coordinates named beside this one); the three read --json-cache and accept nothing without it, "
        "plugin_value_follows_splice_term_past_the_record_span (a DownstreamProtein, a ProteinLengthChange, an "
        "NMD, a LoF, a LoF_info, a LoF_filter or a LoF_flags that one side alone writes, the side whose terms "
        "let the plugin run, on a row whose two Consequence sets differ as splice_term_past_the_record_span "
        "admits, with the record's alleles where the row carries them or --input gives them and without them "
        "where neither does), "
        "plugin_value_follows_stop_codon_term_pair (a DownstreamProtein, a ProteinLengthChange, an NMD, a LoF, "
        "a LoF_info, a LoF_filter or a LoF_flags that one side alone writes, the side whose terms let the plugin run, on a "
        "row whose two Consequence sets differ as one of the stop codon classes "
        "stop_retained_without_stop_codon_test, stop_codon_overlap_insertion_range and "
        "stop_retained_three_base_stop_codon_test admits), "
        "blosum62_hgvs_shift_shared_across_alleles (with --json-cache, a BLOSUM62 the reference alone writes on a later "
        "allele's row whose CDS and protein positions hgvs_shift_shared_across_alleles takes, where the reference's "
        "Amino_acids is one residue over one, the first the translation of the moved pair's codon, the candidate's own "
        "pair is none the plugin scores, and the value is the matrix cell of the reference's pair), "
        "plugin_value_on_reference_skipped_record (a pHaplo, a pTriplo or a count of existing upstream open "
        "reading frames that the candidate alone writes on a row of a record Ensembl VEP dropped; such a row has no "
        "reference entry, so this comparison sets it aside whole under reference_skipped_record and takes no "
        "pair for the class, which the plugin tuple scorer decides), "
        "the structural-variant comparator's classes reference_skipped_record (a row the candidate alone "
        "writes on a record Ensembl VEP drops before annotation, an unsupported type or a span over "
        "--max_sv_size, and wrote nothing for; needs --input in every format), "
        "breakend_mate_local_read (a bracket-mate row whose consequence sets are the local-coordinate "
        "and mate-coordinate derivations over the cache's transcript model, skipped whole; needs "
        "--json-cache) and bnd_own_chromosome_rows_lost (an own-end transcript row the candidate alone "
        "writes at a Location where the reference kept the mate rows of another chromosome and no own-end "
        "row, the transcript within 5,000 bases on the Location's chromosome; needs --json-cache), "
        "dosage_sensitivity_regulatory_row_warnings (--warnings mode: the reference's Plugin "
        "'DosageSensitivity' went wrong lines), fork_warning_pid_prefix (--warnings mode: the "
        "<pid> : Ensembl VEP's parent wrote before a forked child's line is removed from the "
        "reference's lines), the manifest-applied classes (each read from "
        "the manifest beside the input) transcript_annotator_sibling_allele_rows, "
        "transcript_annotator_insertion_window (the table's fields the candidate alone fills on "
        "an insertion's row, from a row at the insertion's start that Ensembl VEP's two-base "
        "window never reaches), "
        "hgvs_shift_shared_across_alleles (with --json-cache also decided without a manifest on the row of a "
        "second or later alternate allele of its record: its CDS_position, Protein_position and Amino_acids "
        "where the reference prints the record's own coordinates moved by the HGVS_OFFSET of the first "
        "allele's row and the candidate the record's own), hgvs_shift_sequence_edge_window, "
        "frequency_filter_allele_key_order, "
        "minimal_symbolic_alleles (the keys of a structural record with two or more symbolic "
        "ALTs under --minimal: a two-ALT record's reference rows at the trimmed span, a longer "
        "record's first two alleles' rows marked MINIMISED, and the candidate's rows at the "
        "record's span), shift_3prime_shared_predicates (a row whose terms Ensembl VEP reads at "
        "another allele's shifted span under --shift_3prime 1) and shift_start_codon_test (a "
        "boundary deletion's terms and HGVSp under --shift_3prime 1) "
        "and, in the recoder format, "
        "refseq_transcript_entries and lrg_transcript_entries (the RefSeq and LRG notations "
        "dropped on both sides), archive_dbsnp_synonyms (the Archive dbSNP synonym entries "
        "dropped) and identifier_indel_transcript_notations (hgvsc, hgvsp and mane_select "
        "dropped for an identifier input's insertion or deletion alleles)",
    )
    ap.add_argument(
        "--keep-only", choices=(REGULATORY_FEATURE_ROWS,),
        help="comparison mode: retain only RegulatoryFeature/MotifFeature rows and their records on both sides",
    )
    ap.add_argument(
        "--floors",
        help="COL=MIN,... minimum populated counts per default-format column, both sides; in JSON "
        "counted under the column's JSON key or, for a plugin's column, inside the plugin's object, a "
        "column with no row count there (record-level, Feature_type, Extra, a custom source's "
        "columns) reported as not applicable; in tab and VCF a column neither header declares is "
        "not applicable; the module docstring states the rules",
    )
    ap.add_argument("--field-pairs", action="append", help="COL=EXPECTED=>ACTUAL[;...] value pairs accepted as equal")
    ap.add_argument(
        "--nearest-sites", type=Path, help="the candidate JSON cache whose start sites settle NEAREST values"
    )
    ap.add_argument(
        "--json-cache",
        type=Path,
        help=(
            "the candidate JSON cache whose transcript models and spans the classes "
            "breakend_mate_local_read and bnd_own_chromosome_rows_lost read"
            "; hgvs_shift_sequence_edge_window reads its transcripts beside --reference-fasta"
            "; hgvs_shift_shared_across_alleles reads its mapper pairs, coding sequences and codon tables"
            "; organelle_peptide_paths_translation_table reads its codon tables"
            "; downstream_utr_frame_past_coding_end and downstream_start_past_own_codon read its coding sequences, "
            "3' UTRs and codon tables, as plugin_hgvs_shifted_cds_coordinates does with the transcripts' strands, "
            "exons, coding bounds and mapper pairs"
            "; hgvs_clipped_reference_against_whole_alternate counts its exon spans for a substitution's row "
            "whose two notations show the clamp (without it that arm takes no row)"
            "; hgvsp_intron_end_insertion and hgvsp_intron_start_insertion read its exons and coding bounds"
            "; maxentscan_overlapping_windows reads its exons and strands (without it the class takes no row)"
            "; spliceai_multi_entry_cell reads its transcripts' symbols for the gene of a row that writes no SYMBOL"
        ),
    )
    ap.add_argument(
        "--reference-fasta",
        type=Path,
        help="the reference FASTA both engines read, uncompressed beside its .fai: under --masks "
        "hgvs_shift_sequence_edge_window, with --json-cache, the class is decided from its sequence "
        "and the cache's transcripts where no manifest documents the row",
    )
    ap.add_argument(
        "--gwas-catalog",
        type=Path,
        help="the GWAS catalog associations file the run's GWAS plugin read, whose rows of two or more SNPs "
        "the classes gwas_haplotype_cross_product and gwas_interaction_first_snp_allele_whitespace read",
    )
    ap.add_argument(
        "--revel-table",
        type=Path,
        help="the REVEL table the run's REVEL plugin read, with its tabix index beside it, whose rows the "
        "class revel_other_reference_amino_acid_rows reads",
    )
    ap.add_argument(
        "--revel-no-match",
        action="store_true",
        help="the run set the REVEL plugin's no_match, so the class reads no row's transcript ids",
    )
    ap.add_argument(
        "--spliceai-table",
        type=Path,
        action="append",
        help="a score file the run's SpliceAI plugin read, with its tabix index beside it (repeatable: the snv= "
        "file and the indel= file), whose records the class spliceai_multi_entry_cell reads at a row's variant; "
        "the gene of a row that writes no SYMBOL is read from --json-cache",
    )
    ap.add_argument(
        "--utrannotator-sweep",
        type=Path,
        action="append",
        help="a default-format Ensembl VEP output of the same input under a fixed Perl hash seed (repeatable, one "
        "per seed), whose rows the classes utrannotator_consequence_precedence and "
        "utrannotator_gained_start_remade_in_place read",
    )
    ap.add_argument("--default", type=Path, help="self-parity: the default-format output")
    ap.add_argument("--tab", type=Path, help="self-parity: the --tab output")
    ap.add_argument("--vcf", type=Path, help="self-parity: the --vcf output")
    ap.add_argument("--json", type=Path, help="self-parity: the --json output")
    ap.add_argument("--parquet", type=Path, help="self-parity: the Parquet output (file or directory)")
    ap.add_argument("--expected-skipped", type=Path, help="warnings: the expected --skipped_variants_file output")
    ap.add_argument("--actual-skipped", type=Path, help="warnings: the actual --skipped_variants_file output")
    ap.add_argument(
        "--unordered",
        action="store_true",
        help="compare the shared entries without their order (a run whose reference rows Ensembl VEP "
        "emits in hash order: --pick_allele, --pick_allele_gene, --per_gene)",
    )
    ap.add_argument("--report", type=Path, help="where the JSON report is written")
    ap.add_argument("--examples", type=int, default=EXAMPLES_PER_GROUP, help="examples kept per difference group")
    return ap


def run(args: argparse.Namespace, front_end: FrontEnd | None = None) -> dict:
    if args.keep_only is not None and (args.warnings or args.self_parity):
        raise ComparatorError("--keep-only is available only in --format comparison mode")
    if args.warnings:
        if args.expected is None or args.actual is None:
            raise ComparatorError("--warnings needs --expected and --actual")
        return compare_warnings(args.expected, args.actual, args.expected_skipped, args.actual_skipped, parse_masks(args.masks))
    if args.self_parity:
        missing = [name for name in ("default", "tab", "vcf", "json") if getattr(args, name) is None]
        if missing:
            raise ComparatorError("--self-parity needs --" + ", --".join(missing))
        inputs = input_records(args.input) if args.input is not None else None
        parquet_text = parquet_as_tab(args.parquet, canonical=False) if args.parquet is not None else None
        report = SelfParity(
            read_text(args.default),
            read_text(args.tab),
            read_text(args.vcf),
            read_text(args.json),
            parquet_text,
            inputs,
            args.examples,
        ).run()
        report["files"] = {
            name: str(getattr(args, name)) for name in ("default", "tab", "vcf", "json", "parquet", "input")
            if getattr(args, name) is not None
        }
        return report
    if args.format is None or args.expected is None or args.actual is None:
        raise ComparatorError("--format, --expected and --actual are required")
    masks = parse_masks(args.masks)
    documented = Documented.from_manifest(
        load_manifest(args.manifest),
        masks,
        args.input.name if args.input is not None else None,
        run=args.run,
        reference_sequence=args.reference_sequence,
    )
    structural_variants_requested = args.format != "recoder" and any(k in masks for k in STRUCTURAL_VARIANT_CLASSES)
    # The records are read where a parser or class consumes them; a default or tab comparison
    # names its input for the manifest's flag set alone and reads none of its lines, but for
    # the records a class of ``RECORD_REFERENCE_CLASSES`` reads by position.
    reads_records = (
        args.format in ("vcf", "json") or structural_variants_requested or any(k in masks for k in ALLELE_READING_CLASSES)
    )
    inputs = input_records(args.input) if args.input is not None and reads_records else None
    floors = parse_floors(args.floors)
    field_pairs = parse_field_pairs(args.field_pairs)
    nearest_sites = NearestSites.from_cache(args.nearest_sites) if args.nearest_sites is not None else None
    structural_variants = None
    if structural_variants_requested:
        structural_variants = StructuralVariants(masks, args.input, args.json_cache)
        if inputs is not None and len(structural_variants.input_variants) != len(inputs):
            raise ComparatorError(
                f"{args.input}: {len(inputs)} data lines, {len(structural_variants.input_variants)} records the "
                "structural-variant comparator reads; a record needs CHROM, POS, ID, REF and ALT"
            )
    sequence_edges = None
    if args.reference_fasta is not None and HGVS_SHIFT_SEQUENCE_EDGE_WINDOW in masks:
        if args.json_cache is None:
            raise ComparatorError(
                f"--masks {HGVS_SHIFT_SEQUENCE_EDGE_WINDOW} with --reference-fasta needs --json-cache, the candidate JSON cache"
            )
        sequence_edges = SequenceEdges(args.reference_fasta, args.json_cache)
    codon_tables = None
    if ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE in masks and args.json_cache is not None and args.format != "recoder":
        codon_tables = OrganelleCodonTables(args.json_cache)
    downstream_sequences = None
    transcript_classes = (*TRANSCRIPT_READING_CLASSES, HGVS_SHIFT_SHARED_ACROSS_ALLELES)
    if args.json_cache is not None and args.format != "recoder" and any(k in masks for k in transcript_classes):
        downstream_sequences = DownstreamSequences(args.json_cache)
    transcript_edges = None
    if args.json_cache is not None and HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE in documented.masks:
        transcript_edges = TranscriptEdges(args.json_cache, inputs)
    intron_ends = None
    intron_gap_classes = [k for k in INTRON_GAP_CLASSES if k in masks]
    if intron_gap_classes:
        if args.json_cache is None:
            raise ComparatorError(f"--masks {','.join(intron_gap_classes)} needs --json-cache, the candidate JSON cache")
        intron_ends = IntronEnds(args.json_cache, inputs)
    splice_windows = None
    if args.json_cache is not None and MAXENTSCAN_OVERLAPPING_WINDOWS in masks and args.format != "recoder":
        splice_windows = SpliceWindows(args.json_cache, inputs)
    return compare_files(
        args.format,
        args.expected,
        args.actual,
        inputs,
        documented,
        masks,
        floors,
        field_pairs,
        args.examples,
        nearest_sites,
        unordered=args.unordered,
        keep_only=args.keep_only,
        structural_variants=structural_variants,
        gwas_pairs=GwasCatalogPairs.from_catalog(args.gwas_catalog) if args.gwas_catalog is not None else None,
        utrannotator_sweep=UtrannotatorSeedSweep.from_outputs(args.utrannotator_sweep) if args.utrannotator_sweep else None,
        sequence_edges=sequence_edges,
        codon_tables=codon_tables,
        downstream_sequences=downstream_sequences,
        transcript_edges=transcript_edges,
        intron_ends=intron_ends,
        front_end=front_end,
        input_path=args.input if any(k in masks for k in RECORD_REFERENCE_CLASSES) else None,
        revel_table=RevelTable(args.revel_table, args.revel_no_match) if args.revel_table is not None else None,
        splice_windows=splice_windows,
        spliceai_table=SpliceAiTable(args.spliceai_table, args.json_cache) if args.spliceai_table else None,
    )


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        front_end = FrontEnd.named(os.environ.get(FRONT_END_VARIABLE))
        report = run(args, front_end)
        write_report(report, args.report)
    except ComparatorError as exc:
        print(f"ERROR: [{COMPONENT}] {exc}", file=sys.stderr)
        return 2
    if front_end is not None:
        print(f"[{COMPONENT}] {front_end.note}", file=sys.stderr)
    label = report.get("format") or report["mode"]
    print(f"{label}: {report['verdict']}")
    for reason in report["reasons"]:
        print(f"  {reason}")
    return 0 if report["verdict"] == "green" else 1


if __name__ == "__main__":
    sys.exit(main())
