#!/usr/bin/env python3
"""Build a golden corpus: exemplar input records for every consequence-set combination a
VEP reference output contains, keyed by Ensembl release, species and assembly.

The corpus is the fixture behind `crates/vep-cli/tests/golden.rs` and
`crates/vep-cli/tests/format_parity.rs`: a small VCF whose records, annotated by VEP and by
vep-rs against the same transcript cache, must yield the same rows in every output format.

Subcommands:

  select    Stream one or more VEP default-format reference outputs, enumerate the distinct
            consequence sets they contain, choose up to K exemplar records per set
            (stratified by variant class and transcript strand, and by the family key of
            `--family-columns` when given; when `--cache` names the JSON cache the fixture
            is pruned from, the records whose transcripts cost the fixture cache the fewest
            bytes come first, then the records that exercise transcript flags, coordinate
            ranges and unmappable ends; every variant class and every family stratum seen
            keeps at least one exemplar), resolve each exemplar back to its input VCF
            record, and write `variants.vcf` plus `manifest.json`. A row's variant class is
            read from its Location and Allele and, where those print a two-base substitution
            and a two-base insertion alike, from its input record's REF and ALT
            (`variant_class`). The manifest's `covers`
            list names the parity-matrix rows the corpus exercises: `--covers <item>`
            (repeatable) writes it, each name a row of `--matrix` when one is given, a corpus
            regenerated over an existing manifest keeps that manifest's list, its `masks`
            (each one of `MASK_CLASSES`, else the run is refused before anything is
            written), its `reference_warnings` (the texts Ensembl VEP wrote about the
            run rather than about one record), its `perl_mode` (how the manifest's own
            reference was generated, one of `PERL_MODES`; `offline` when absent, and
            the mode a flag set naming none inherits), its `plugin_files` (the data files an
            Ensembl VEP plugin reads, per plugin, paths inside the corpus that the
            reference generator checks and digests) and its `flag_sets` (each set's run
            definition and `covers`; the divergence keys `classify --flag-set` wrote
            describe the old records and are dropped, as the manifest's own are), and a
            run with neither
            writes nothing. `--include <sid>=<name>` (repeatable) names a record of suite
            `sid`'s input, by its ID or by VEP's name for an ID-less record, that the
            corpus carries whatever the sampling chooses: its rows in the suite's
            reference join the chosen exemplars after the sampling, so a record chosen
            for what it co-locates with rather than for its consequence set stays in
            even when another record of the same stratum is sampled; the names are
            recorded as `included_records`. `--keep-records` keeps every record of the
            manifest being regenerated the same way, each resolved again from its
            suite's input and rowed again from its suite's reference, so a regeneration
            loses no record (the count is recorded as `kept_records`, and the manifest's
            `included_records` are carried into the new map as named records). A named or
            kept record is an exemplar of every combination its reference rows show, so
            `variant_classes` and `family_strata` count each of those rows, and it survives
            `--max-records`; one with no row in its suite's reference is an error, and
            nothing is written, so a corpus whose manifest carries a record the reference
            skipped (`reference_warning`) is regenerated without `--keep-records`, which
            keeps no record by name: that record returns only when the sampling picks it.
            `--keep-samples N` keeps `FORMAT` and the first N sample
            columns of every record and of the `#CHROM` header, and adds records until
            every genotype class the inputs show among those samples has
            `--genotype-exemplars` records (`--require-genotype-class` names a class the
            corpus must carry whether or not the inputs show it). `--record-classes
            regulatory_class` classes every record by the regulatory rows of its reference
            (`regulatory_classes`: which of the two feature types it reaches, whether it
            reaches no transcript, and the shape of each motif allele) and adds records until
            every class the references show has `--record-class-exemplars` records.
  classify  Compare a VEP default-format output for `variants.vcf` with a vep-rs output for
            the same file and record, per exemplar row VEP emits that vep-rs does not, the
            documented divergence class that explains it (or `unexplained_residual`), so the
            golden test asserts the documented difference instead of equality. A record no
            VEP row names is one VEP dropped before annotation: it must carry
            `reference_warning` (VEP's reason, from its warnings file, recorded when the
            corpus is built), it gets `reference_rows` 0, and every key only vep-rs emits
            for it is classed `reference_skipped_record`; any other key only vep-rs emits is
            `unexplained_residual`. A zero-row record without the warning, or a warned record
            with rows, is an error, unless the record carries no row by its nature: a
            non-variant record (ALT `.`), one the run's expected JSON names as a feature
            with no allele to annotate (`most_severe_consequence` `?`: every ALT a `*` or
            `<DEL:*>`, equal to the reference or a repeat) or, under a row filter of the run
            (`--coding_only`, `--no_intergenic`) or its frequency filter (`--check_frequency`,
            `--filter_common`), one vep-rs writes no row for either. Under the frequency
            filter a record Ensembl VEP dropped and vep-rs kept because an allele's frequency
            is exactly 0 (every vep-rs row of it lists a FREQS item of value 0) has its keys
            classed `frequency_filter_zero_value`. The
            records are the run's input's, numbered as the golden tests number them: the
            manifest's for `variants.vcf`, and for a flag set over another `input` that
            file's lines, a line the manifest records (a twin of `variants.vcf`) keeping
            its entry and a file in another format or delimiter, `variants.vcf` re-encoded
            record for record, taking the manifest's records in order. A row under an
            identifier is read with the alleles of the record that writes it, found among
            those records by position, where a class reads them (`record_reference`).
            A run whose `masks` list
            `regulatory_feature_rows` has its VEP rows of `Feature_type` RegulatoryFeature
            and MotifFeature left out first, the rows the golden test sets aside for a run
            under `--regulatory`; one whose `masks` list `order_only` has its `DOMAINS`,
            `FREQS` and `ZYG` values compared as multisets of their items under `--field`,
            as the golden test compares them. A field divergence two causes change (a
            `DownstreamProtein` Ensembl VEP reads at its HGVS-shifted span on a row whose own
            span opens on a codon's third base) names the second cause's class under
            `second_expected_divergence`, and the format comparator applies the entry only
            where both classes are named. A `divergences` entry may carry, written by hand,
            `vep_rs_consequence_sets_with_reference_sequence`: vep-rs's sets on the key in a
            run that has a reference sequence, where `vep_rs_consequence_sets` are those of a
            run without one. classify never writes the key; it carries a key the target run
            holds onto the entry it writes for the same Location, Allele, Feature, records
            and Ensembl VEP set, and refuses the run when it writes no such entry or the key
            is of another shape. `covers` and `flag_sets` are kept; with `--flag-set
            <name>` the divergence keys are written into that flag set instead of the
            manifest's top level, so each flag set documents its own run, under its own
            `masks` when its entry carries the key and the manifest's otherwise.
  check     Recompute the observed combination space from the reference outputs and fail
            when a combination has no exemplar in the manifest; with `--matrix` fail when a
            `covers` name (the manifest's or a flag set's) is not a matrix row; with
            `--family-columns` fail when a family stratum the references show has no
            exemplar among the corpus's records; on a corpus built with `--keep-samples`
            fail when `variants.vcf` does not carry `FORMAT` and the manifest's samples on
            its header and every record, or when a genotype class the inputs show among
            those samples (or one the manifest or `--require-genotype-class` requires) has
            no record in the corpus; on a corpus built with `--record-classes` fail when a
            class the references show has no record in the corpus.
  covers    List (`--list`) or extend (`--add <item>`, repeatable) the `covers` of a corpus,
            or of one of its flag sets (`--flag-set <name>`); a name that is not a row of
            `--matrix` is refused and nothing is written.
  flag-set  Write one `flag_sets` entry: a run of the corpus under its own `flags` beside
            the manifest's default run, with its own expected outputs. An entry carries
            `name`, `flags` (the whole flag list of that run), `input` (a file in the corpus,
            `variants.vcf` by default), `stdin` (the input piped on standard input, with
            `-i STDIN`; its `flags` must name `--format`, which Ensembl VEP needs on
            STDIN), `unordered` (the reference's rows within a record come out in
            Ensembl VEP's hash order, so the golden tests compare the rows and not
            their order), `formats` (the output formats compared, all four by default),
            `expected_dir` (`expected_<name>` by default), `perl_mode` (written as
            `database` when the reference was generated against the Ensembl database
            rather than offline, `gtf` or `gff` when Ensembl VEP read its transcripts from
            the file given to `--gtf` or `--gff`; the manifest's own `perl_mode` when the
            option is not given) and `covers`; `input` and `expected_dir` are paths inside
            the corpus directory. A set that exists has those keys rewritten from the
            options and keeps every other key of its entry, its hand-written `masks` (the
            classes of its run alone, each one of `MASK_CLASSES`, as the manifest's) and
            `reference_warnings` (the texts Ensembl VEP wrote about its run alone) among them;
            its divergence keys describe the run before the rewrite and go, so
            `classify --flag-set` runs after it.

Row keys follow `compare_vep_outputs.py`: (Location, Allele, Feature, Feature_type,
consequence set), consequence terms sorted and comma-joined. A record is identified by
VEP's (Uploaded_variation, Location) pair, which every row of one input record shares.

The family key of `--family-columns` is one class per named Extra column plus the variant
class: a bare name (`HGVSc`) classes a row by the column's presence (`HGVSc+`, `HGVSc-`), a
name ending in `=` (`IMPACT=`) by its value (`IMPACT=HIGH`); `check --family-columns`
asserts every stratum the references show is covered, which is the per-family coverage
claim behind a flag family's corpus. `custom_match=<exact>/<overlap>` names the short names
of one `--custom` file read as a `type=exact` and as a `type=overlap` source and classes a
row by how they answered (`custom_match_class`): `exact`, `multi_allelic` (the exact source
matched one allele of a record whose fields list several), `overlap` (the overlap source
alone matched), `sv_span` (the overlap source alone matched a structural variant whose span
covers the record) or `miss`. `<a>/<b>` names two columns and classes a row by their
agreement (`agreement_class`): `same`, `differ` or `-` when either is absent, so
`GIVEN_REF/USED_REF` keeps a row whose transcript sequence differs from the input's reference
allele beside the rows where the two agree. A list-valued column stratified by value
(`REFSEQ_MATCH=`) makes one stratum per distinct list, so every code the lists carry has an
exemplar.

The genotype classes of `--keep-samples` are read from the kept samples' GT values as
Ensembl VEP reads a genotype for `--individual` (ensembl-io BaseVCF4.pm 941-951 and
Parser/VCF.pm 627-651 at release 116): per sample `het`, `hom_alt`, `hom_ref`, `missing`,
`partially_missing`, `haploid`, `star` and `multi_allelic`; per record `all_missing`,
`no_carrier` and `one_hom_ref_rest_missing` (`genotype_labels` and
`record_genotype_classes` define each; a structural record shows none, since Ensembl VEP
annotates it without genotypes). A record shows every class any kept sample shows, the
manifest lists each record's classes and the count of records per class, and `check`
asserts every class the inputs show is covered.

The record classes of `--record-classes regulatory_class` are read from a record's reference
rows under `--regulatory` (`regulatory_classes` defines each): `rf_only`, `mf_only` or
`rf_and_mf` for the regulatory feature types its rows name, `regulatory_only` for a record
with such a row and no Transcript row, and per MotifFeature row `motif_snv`, `motif_mnv` or
`motif_indel` by the shape of its allele, read with the input record's REF and ALT
(`allele_shape`: a two-base substitution and a two-base insertion print the same row), and
`motif_reverse` for a motif on the reverse strand. A record shows every class any of its rows
shows; the manifest lists each record's
classes and the count of records per class, `select` adds records until every class the
references show has `--record-class-exemplars` records, and `check` asserts every class the
references show is covered.

Selection is deterministic under `--seed`: reservoir sampling bounds the candidates kept per
combination and the same seed reproduces the same corpus from the same references.
"""
from __future__ import annotations

import argparse
import gzip
import json
import random
import re
import sys
import zlib
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Collection, Iterable, Iterator

GENERATOR_VERSION = "2"
RESERVOIR = 64
FAMILY_RESERVOIR = 8
SEQUENCE_CLASSES = ("snv", "deletion", "insertion", "mnv_or_complex")
OUTPUT_FORMATS = ("default", "tab", "vcf", "json")
# How a reference was generated: Ensembl VEP on its cache, connected to its database, or
# reading its transcripts from the file given to `--gtf` or `--gff` (the golden harness then
# sets the source's `SOURCE` column and header lines aside, since vep-rs reads the cache
# vep-cache-builder built from the same file, which carries no source name).
PERL_MODES = ("offline", "database", "gtf", "gff")
DEFAULT_INPUT = "variants.vcf"
# The species of a manifest that names none; `--species` writes another, in Ensembl's
# production spelling, the cache directory's name.
DEFAULT_SPECIES = "homo_sapiens"
SPECIES_RE = re.compile(r"[a-z][a-z0-9_]*")
# The manifest `masks` class under which the reference's regulatory rows are left out, and
# the `Feature_type` values of those rows.
REGULATORY_FEATURE_ROWS = "regulatory_feature_rows"
REGULATORY_FEATURE_TYPES = frozenset({"RegulatoryFeature", "MotifFeature"})
# The `masks` class under which the set-valued columns agree as a multiset of their items
# in any order (Ensembl VEP prints DOMAINS in the cache's protein_features order,
# BaseTranscriptVariation.pm:634-643, FREQS in hash order, OutputFactory.pm:1152-1153, and
# the ZYG list of --individual_zyg in hash order, OutputFactory.pm:983), and those columns;
# the default format joins the items with `,` (Utils.pm:185-198).
ORDER_ONLY = "order_only"
ORDER_ONLY_COLUMNS = ("DOMAINS", "FREQS", "ZYG")
# The classes the golden harness alone gives a meaning (its warnings, header and column
# comparisons, which `classify` does not make): the reference's warnings about a plugin module
# vep-rs cannot load, the `phlyoP<N>way` key Ensembl VEP's --phyloP set misspells
# (Config.pm:555), its two descriptions of a custom source's PC field
# (AnnotationSourceAdaptor.pm:395-397), and the SOURCE column and source header lines of a
# reference read from a GTF or GFF, which the harness names `file_cache_writes_no_source_field`
# and sets aside before the comparison (the cache vep-cache-builder builds from the same file
# carries no source name). The harness's `MASK_CLASSES` names these four after
# `REGULATORY_FEATURE_ROWS` and `ORDER_ONLY`, and nothing else.
HARNESS_ONLY_MASK_CLASSES = (
    "unloaded_module_plugin", "conservation_short_name_spelling", "pc_description_text",
    "file_cache_writes_no_source_field",
)
# The classes of the behaviours vep-rs corrects where Ensembl VEP 116.2 is wrong (each a page
# or section of docs/intended-divergences.md): a corpus documents their rows under
# `expected_divergence` in its three lists and names the class in `masks`; the format
# comparator's built-in classes carry the same names, and the golden harness's
# `CORRECTED_BEHAVIOUR_CLASSES` names every class a committed corpus names in `masks`.
CORRECTED_BEHAVIOUR_CLASSES = (
    "start_cooccurrence_swap", "splice_family_swap", "gxf_reverse_strand_utr_sequence",
    "conservation_max_text_order", "gwas_haplotype_cross_product",
    "gwas_interaction_first_snp_allele_whitespace",
    "alphamissense_two_digit_transcript_version", "cadd_structural_allele_match",
    "splicevault_sibling_allele_rows", "splicevault_json_later_rows",
    "transcript_annotator_sibling_allele_rows", "transcript_annotator_insertion_window",
    "utrannotator_annotation_pair_order", "nearest_exon_jb_item_order",
    "utrannotator_gained_start_window",
    "utrannotator_consequence_precedence", "utrannotator_structural_allele_text",
    "utrannotator_lost_start_after_insertion",
    "utrannotator_gained_start_remade_in_place", "utrannotator_frameshift_at_stop_boundary",
    "utrannotator_lost_stop_written_before_stop_codon",
    "loeuf_zero_minimum", "maxentscan_overlapping_windows",
    "alphamissense_sibling_allele_rows", "revel_sibling_allele_rows",
    "revel_other_reference_amino_acid_rows",
    "gnomadmt_info_value_second_equals", "loftee_anc_allele_multiallelic_record",
    "loftee_first_alt_donor_filter", "loftee_gerp_zoom_level_summary",
    "nearest_exon_jb_stretch_priming", "nearest_exon_jb_shifted_span",
    "nearest_exon_jb_id_less_exon_fold", "downstream_reverse_strand_rotation",
    "spliceai_chopped_multiallelic", "spliceai_multi_entry_cell", "plugin_contig_mitochondrial_spelling",
    "dosage_sensitivity_regulatory_row_warnings", "hgvs_shift_shared_across_alleles",
    "hgvs_shift_sequence_edge_window", "refseq_edited_model_hgvsc_offset",
    "refseq_edited_model_shifted_insertion_bases",
    "refseq_edited_model_seam_insertion_hgvsc",
    "refseq_edited_model_unminimised_allele_hgvsc",
    "refseq_edited_model_hgvs_offset", "frequency_filter_zero_value",
    "frequency_filter_allele_key_order", "organelle_peptide_paths_translation_table",
    "minimal_symbolic_alleles", "shift_3prime_shared_predicates", "shift_start_codon_test",
    "hgvs_shifted_insertion_clamped_end", "shift_intron_term_unshifted_span",
    "hgvs_intron_offset_abutting_insertion",
    "hgvs_shifted_indel_written_span_gate", "hgvs_shifted_deletion_clamped_span",
    "hgvs_clipped_reference_against_whole_alternate",
    "hgvs_added_bases_left_at_the_span_end",
    "hgvs_insertion_left_where_the_shift_stops",
    "hgvs_added_bases_written_outside_the_transcript",
    "hgvsp_stop_codon_protein_comparison",
    "hgvsp_stop_named_on_incomplete_coding_sequence",
    "hgvsp_shifted_indel_written_span_gate",
    "hgvsp_intron_end_insertion",
    "hgvsp_intron_start_insertion",
    "refseq_edited_model_deleted_span_codons", "hgnc_id_copy",
    "downstream_utr_frame_past_coding_end", "downstream_start_past_own_codon",
    "plugin_hgvs_shifted_cds_coordinates", "plugin_value_follows_stop_codon_term_pair",
    "plugin_value_follows_splice_term_past_the_record_span",
    "plugin_value_on_reference_skipped_record", "blosum62_hgvs_shift_shared_across_alleles",
    "stop_retained_without_stop_codon_test",
    "stop_codon_overlap_genomic_span",
    "stop_retained_three_base_stop_codon_test",
    "stop_codon_removed_without_stop_codon_test",
    "stop_codon_overlap_insertion_range",
    "splice_term_past_the_record_span",
)
# Every class a `masks` list, the manifest's or a flag set's, may name.
MASK_CLASSES = (REGULATORY_FEATURE_ROWS, ORDER_ONLY, *HARNESS_ONLY_MASK_CLASSES, *CORRECTED_BEHAVIOUR_CLASSES)
CLASSIFY_KEYS = ("divergences", "vep_rs_only_tuples", "divergence_summary", "vep_rs_only_summary",
                 "fields_compared", "field_divergences", "field_divergence_summary")
DEFAULT_COLUMNS = (
    "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
    "Consequence", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "Extra",
)
FIXED_COLUMNS = ("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")


def opener(path: Path):
    """`gzip.open` for a `.gz` or `.bgz` (BGZF is gzip), else `open`."""
    return gzip.open if str(path).endswith((".gz", ".bgz")) else open


def normalize_consequence_set(raw: str) -> str:
    return ",".join(sorted({t for t in raw.strip().split(",") if t}))


def parse_extra(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not raw or raw == "-":
        return out
    for part in raw.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        elif part:
            out[part] = ""
    return out


def variant_class(location: str, allele: str, substitution: bool = False) -> str:
    """Infer VEP's variant shape from the Location and Allele columns of one row: one of
    `SEQUENCE_CLASSES` for a sequence allele; `breakend`, or `breakend_insertion` when the
    bracket allele carries inserted sequence (more than the anchor base outside the brackets);
    `single_breakend`; and `symbolic_<word>` for the class word VEP prints for a symbolic
    allele (`symbolic_deletion`, `symbolic_Alu_insertion`). Each is annotated by a path of its
    own, so each is a stratum of its own. A two-base substitution and a two-base insertion
    print the same Location and Allele (Ensembl VEP writes an insertion between its two
    flanking bases), so a span of two adjacent positions is `mnv_or_complex` when
    `substitution` says the input record holds the allele as a substitution
    (`holds_allele`), and `insertion` otherwise."""
    if "[" in allele or "]" in allele:
        bases = re.sub(r"[\[\]][^\[\]]*[\[\]]", "", allele)
        return "breakend_insertion" if len(bases) > 1 else "breakend"
    if allele.startswith(".") or allele.endswith("."):
        return "single_breakend"
    if allele and not re.fullmatch(r"[ACGTNacgtn*-]+", allele):
        return f"symbolic_{allele}"
    span = location.split(":", 1)[1] if ":" in location else location
    if allele == "-":
        return "deletion"
    if "-" in span:
        a, b = span.split("-", 1)
        if a.isdigit() and b.isdigit() and int(b) == int(a) + 1:
            return "mnv_or_complex" if substitution else "insertion"
        return "mnv_or_complex"
    return "snv" if len(allele) == 1 else "mnv_or_complex"


CUSTOM_MATCH = "custom_match"
CUSTOM_MATCH_SPEC = re.compile(CUSTOM_MATCH + r"=([A-Za-z0-9_.-]+)/([A-Za-z0-9_.-]+)")
CUSTOM_MATCH_CLASSES = ("exact", "multi_allelic", "overlap", "sv_span", "miss")
# `<a>/<b>`: two columns classed by whether they agree (`agreement_class`).
AGREEMENT_SPEC = re.compile(r"([A-Za-z0-9_.-]+)/([A-Za-z0-9_.-]+)")
AGREEMENT_CLASSES = ("same", "differ", "-")


def parse_family_columns(arg: str | None) -> list[str]:
    """`--family-columns HGVSc,IMPACT=,custom_match=CE/CO,GIVEN_REF/USED_REF`: a column name
    (presence), a name ending in `=` (value), `custom_match=<exact>/<overlap>`
    (`custom_match_class`) or `<a>/<b>`, two columns classed by their agreement
    (`agreement_class`)."""
    specs = [s.strip() for s in (arg or "").split(",") if s.strip()]
    for spec in specs:
        if (
            not re.fullmatch(r"[A-Za-z0-9_.-]+=?", spec)
            and not CUSTOM_MATCH_SPEC.fullmatch(spec)
            and not AGREEMENT_SPEC.fullmatch(spec)
        ):
            raise SystemExit(
                f"ERROR: [build_golden_corpus] --family-columns wants NAME, NAME=, custom_match=EXACT/OVERLAP or A/B, got {spec!r}"
            )
    return specs


def _present(extra: dict[str, str], name: str) -> bool:
    return extra.get(name, "") not in ("", "-")


def agreement_class(extra: dict[str, str], a: str, b: str) -> str:
    """How two columns of one row relate: `same` when both are present with one value,
    `differ` when both are present with two (`GIVEN_REF` and `USED_REF` on a row whose
    transcript sequence differs from the input's reference allele), `-` when either is
    absent."""
    if not (_present(extra, a) and _present(extra, b)):
        return "-"
    return "same" if extra[a] == extra[b] else "differ"


def custom_match_class(extra: dict[str, str], exact: str, overlap: str, location: str, allele: str) -> str:
    """How one `--custom` file answered for a row through two sources over it: `exact` and
    `overlap` are the short names of its `type=exact` and `type=overlap` sources, the columns
    their record names are written under (`<name>_<field>` for a requested field, one value per
    matched record, comma-joined). `exact`: the exact source matched. `multi_allelic`: it matched
    and a field of the overlap source lists more values than the source lists records while the
    exact source lists fewer, so a record carries one value per allele and the exact match kept
    the matching allele's (`FILTER` lists a record's filters, never alleles, and is passed over).
    `overlap`: the overlap source alone matched a sequence variant. `sv_span`: it alone matched a
    structural variant (a symbolic or breakend allele), whose span covers the record. `miss`:
    neither matched."""
    if _present(extra, exact):
        records = extra[overlap].count(",") + 1 if _present(extra, overlap) else 1
        prefix = f"{overlap}_"
        for key, value in extra.items():
            if not key.startswith(prefix):
                continue
            field = key[len(prefix):]
            if field == "FILTER":
                continue
            overlap_values = value.count(",") + 1
            exact_values = extra.get(f"{exact}_{field}", "").count(",") + 1
            if overlap_values > records and exact_values < overlap_values:
                return "multi_allelic"
        return "exact"
    if _present(extra, overlap):
        return "overlap" if variant_class(location, allele) in SEQUENCE_CLASSES else "sv_span"
    return "miss"


def family_classes(extra: dict[str, str], specs: list[str], location: str = "", allele: str = "") -> tuple[str, ...]:
    """One class per spec: `NAME+` or `NAME-` for presence, `NAME=<value>` for a value (`-` when
    absent), `custom_match=<class>` for a `custom_match=<exact>/<overlap>` spec, `<a>/<b>=<class>`
    for an agreement spec."""
    classes = []
    for spec in specs:
        if match := CUSTOM_MATCH_SPEC.fullmatch(spec):
            classes.append(f"{CUSTOM_MATCH}={custom_match_class(extra, match.group(1), match.group(2), location, allele)}")
        elif match := AGREEMENT_SPEC.fullmatch(spec):
            classes.append(f"{spec}={agreement_class(extra, match.group(1), match.group(2))}")
        elif spec.endswith("="):
            name = spec[:-1]
            classes.append(f"{name}={extra.get(name) or '-'}")
        else:
            classes.append(f"{spec}{'+' if _present(extra, spec) else '-'}")
    return tuple(classes)


SAMPLE_CLASSES = ("het", "hom_alt", "hom_ref", "missing", "partially_missing", "haploid", "star", "multi_allelic")
RECORD_CLASSES = ("all_missing", "no_carrier", "one_hom_ref_rest_missing")
GENOTYPE_CLASSES = SAMPLE_CLASSES + RECORD_CLASSES
STRUCTURAL_ALLELE = re.compile(r"[<\[\]][^*]+[>\]\[]|^\.\w+|\w+\.$")


def genotype_labels(gt: str, alts: list[str], record: str) -> set[str]:
    """The classes one sample's GT shows at a record, read as Ensembl VEP reads a genotype
    for `--individual` (ensembl-io BaseVCF4.pm 941-951 and Parser/VCF.pm 627-651, release
    116): the bits split on `|` when the GT has one, else on `/`; `.` bits dropped before
    anything else is read; an allele containing `*` outside the non-reference set. `missing`
    (no called bit) stands alone; otherwise `partially_missing` when a bit was dropped,
    `haploid` for a one-bit genotype, `star` for a `*` allele, and one of `multi_allelic`
    (two or more distinct non-reference alleles), `hom_ref` (none), `hom_alt` (one, every
    called bit the same) or `het`. A bit that is not an index of the record's alleles is
    refused, naming `record`."""
    bits = gt.split("|") if "|" in gt else gt.split("/")
    called_bits = [b for b in bits if b != "."]
    if not called_bits:
        return {"missing"}
    for b in called_bits:
        # str.isdigit alone accepts a superscript such as "²", which int() rejects.
        if not (b.isascii() and b.isdigit()) or int(b) > len(alts):
            raise failure(f"record {record} has a GT {gt!r} outside its alleles")
    labels: set[str] = set()
    if len(called_bits) < len(bits):
        labels.add("partially_missing")
    if len(bits) == 1:
        labels.add("haploid")
    called = {int(b) for b in called_bits}
    star = {i + 1 for i, a in enumerate(alts) if "*" in a}
    if star & called:
        labels.add("star")
    non_ref = called - {0} - star
    if len(non_ref) > 1:
        labels.add("multi_allelic")
    elif not non_ref:
        labels.add("hom_ref")
    elif len(called) == 1:
        labels.add("hom_alt")
    else:
        labels.add("het")
    return labels


def is_structural_record(ref: str, alts: list[str], info: str) -> bool:
    """Whether Ensembl VEP reads a record as a structural variant, which it annotates before
    the per-individual step and so without any individual's genotype (Parser/VCF.pm 236-244):
    its alleles are not all plain bases and it carries an INFO `SVTYPE` whose value Perl reads
    as true, or a symbolic, breakend or single-breakend allele (`<*>` is none of these). The
    value is read as ensembl-io's `get_info` stores it (BaseVCF4.pm 576-589): the piece after
    the first `=`, so a bare `SVTYPE` flag and `SVTYPE=` hold none, `SVTYPE=0` holds Perl's
    false `0`, and of a repeated key the last one counts. The allele test runs on the
    comma-joined ALT string as Ensembl VEP's does, so a single-breakend allele counts only at
    the string's start or end."""
    if re.fullmatch(r"[ACGT]+", "".join([ref, *alts])):
        return False
    svtype = ""
    for part in info.split(";"):
        pieces = part.split("=")
        if pieces[0] == "SVTYPE":
            svtype = pieces[1] if len(pieces) > 1 else ""
    if svtype not in ("", "0"):
        return True
    return STRUCTURAL_ALLELE.search(",".join(alts)) is not None


def record_genotype_classes(fields: list[str], keep: int) -> set[str]:
    """The genotype classes one record shows among its first `keep` samples: every kept
    sample's `genotype_labels`, plus `all_missing` (no kept sample has a called bit),
    `no_carrier` (a called bit somewhere but no non-reference allele in any kept sample)
    and `one_hom_ref_rest_missing` (one kept sample `hom_ref`, every other one `missing`).
    A structural record shows no class (`is_structural_record`)."""
    record = " ".join(fields[:5])
    if len(fields) < 9 + keep:
        raise failure(f"record {record} has {max(len(fields) - 9, 0)} sample columns, {keep} are kept")
    alts = fields[4].split(",")
    if is_structural_record(fields[3], alts, fields[7]):
        return set()
    keys = fields[8].split(":")
    if "GT" not in keys:
        raise failure(f"record {record} has no GT in its FORMAT")
    gt_index = keys.index("GT")
    classes: set[str] = set()
    missing = hom_ref = 0
    for sample in fields[9:9 + keep]:
        values = sample.split(":")
        labels = genotype_labels(values[gt_index] if gt_index < len(values) else ".", alts, record)
        classes |= labels
        missing += "missing" in labels
        hom_ref += "hom_ref" in labels
    if missing == keep:
        classes.add("all_missing")
    elif missing + hom_ref == keep:
        classes.add("no_carrier")
        if hom_ref == 1:
            classes.add("one_hom_ref_rest_missing")
    return classes


def genotype_class_names(names: list[str] | None) -> list[str]:
    """The `--require-genotype-class` names of a command line, each a genotype class, in order
    and without repeats."""
    names = list(dict.fromkeys(names or []))
    for n in names:
        if n not in GENOTYPE_CLASSES:
            raise failure(f"--require-genotype-class wants one of {', '.join(GENOTYPE_CLASSES)}, got {n!r}")
    return names


# The record-class key of `--record-classes`: a record's regulatory rows under `--regulatory`
# (`regulatory_classes` names its eight classes).
REGULATORY_CLASS = "regulatory_class"
RECORD_CLASS_RESERVOIR = 16


def record_class_key(name: str | None) -> str | None:
    """The `--record-classes` key of a command line, one of `RECORD_CLASS_KEYS`; `None` when
    the option is not given."""
    if name is not None and name not in RECORD_CLASS_KEYS:
        raise failure(f"--record-classes wants one of {', '.join(RECORD_CLASS_KEYS)}, got {name!r}")
    return name


def vcf_lines(path: Path, columns: int | None = None) -> Iterator[tuple[bool, list[str]]]:
    """Every `#CHROM` header line (`True`, every column) and every record of at least eight
    fields (`False`) of a VCF, split on tabs; `##` lines are skipped. With `columns` a record
    is split into that many fields and one more holding the rest of the line, so a source of
    thousands of samples is not split beyond the columns read. A file that cannot be opened
    or read (a `.gz` that is not gzip fails on its first line, not at the open) is refused
    with its path."""
    try:
        with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    yield True, line.rstrip("\n").split("\t")
                    continue
                fields = line.rstrip("\n").split("\t", -1 if columns is None else columns)
                if len(fields) >= 8:
                    yield False, fields
    except (OSError, EOFError, zlib.error) as exc:
        raise failure(f"{path} is not a readable text file ({exc})") from None


def kept_samples(header: list[str], keep: int, path: Path, samples: list[str], source: str = "an earlier input") -> list[str]:
    """The first `keep` sample names of a `#CHROM` header line, refused when the line lacks
    `FORMAT` or that many samples, or when `samples` (the names `source` keeps) names others."""
    if len(header) < 9 or header[8] != "FORMAT":
        raise failure(f"{path} has no FORMAT column, so no sample columns")
    if len(header) < 9 + keep:
        raise failure(f"{path} has {len(header) - 9} sample columns, {keep} are kept")
    kept = header[9:9 + keep]
    if samples and kept != samples:
        raise failure(f"{path} keeps samples {','.join(kept)} where {source} keeps {','.join(samples)}: "
                      "every input must keep the same samples")
    return kept


@dataclass
class Row:
    suite: str
    uploaded_variation: str
    location: str
    allele: str
    feature: str
    feature_type: str
    consequence_set: str
    strand: str
    flags: str
    positions: tuple[str, str, str]
    family: tuple[str, ...] = ()
    # Whether the input record holds the row's allele as a substitution (`iter_rows`): what
    # tells a two-base substitution from a two-base insertion.
    substitution: bool = False

    @property
    def record_key(self) -> tuple[str, str]:
        return (self.uploaded_variation, self.location)

    @property
    def klass(self) -> str:
        return variant_class(self.location, self.allele, self.substitution)

    @property
    def family_stratum(self) -> str:
        """The family classes and the variant class, `;`-joined: the key `--family-columns` covers."""
        return ";".join((*self.family, self.klass))

    def preference(self) -> int:
        score = 0
        if self.flags:
            score += 4
        if any("?" in p for p in self.positions):
            score += 3
        if any("-" in p and p != "-" for p in self.positions):
            score += 2
        return score

    @property
    def span(self) -> int:
        """Genomic width of the Location, in bases (0 for a single position)."""
        coords = self.location.split(":", 1)[-1]
        if "-" not in coords:
            return 0
        a, b = coords.split("-", 1)
        try:
            return abs(int(b) - int(a))
        except ValueError:
            return 0


def iter_rows(
    path: Path,
    suite: str,
    family_specs: list[str] | None = None,
    substitutions: dict[tuple[str, str], list[tuple[int, str]]] | None = None,
) -> Iterator[Row]:
    """Every row of at least seven columns of a default-format reference; a file that cannot
    be opened or read is refused with its path, as `vcf_lines` refuses a VCF. With
    `substitutions` (the reference's input, as `substitution_alts` reads it) a row whose
    input record holds its allele as a substitution is marked so (`holds_allele`), which its
    variant class is read with."""
    try:
        with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 7:
                    continue
                while len(f) < 14:
                    f.append("-")
                extra = parse_extra(f[13])
                chrom, _, span = f[1].partition(":")
                alts = substitutions.get((chrom, f[0])) if substitutions else None
                yield Row(
                    suite=suite,
                    uploaded_variation=f[0],
                    location=f[1],
                    allele=f[2],
                    feature=f[4],
                    feature_type=f[5],
                    consequence_set=normalize_consequence_set(f[6]),
                    strand=extra.get("STRAND", ""),
                    flags=extra.get("FLAGS", ""),
                    positions=(f[7], f[8], f[9]),
                    family=family_classes(extra, family_specs, f[1], f[2]) if family_specs else (),
                    substitution=bool(alts) and holds_allele(span, f[2], alts),
                )
    except (OSError, EOFError, zlib.error) as exc:
        raise failure(f"{path} is not a readable text file ({exc})") from None


def iter_records(
    path: Path, suite: str, substitutions: dict[tuple[str, str], list[tuple[int, str]]] | None = None
) -> Iterator[list[Row]]:
    """The rows of one input record at a time: the consecutive rows of a reference sharing a
    `record_key`, which Ensembl VEP writes together (every row of a variant comes out of one
    output-factory call); `substitutions` is `iter_rows`'s."""
    group: list[Row] = []
    for row in iter_rows(path, suite, substitutions=substitutions):
        if group and row.record_key != group[0].record_key:
            yield group
            group = []
        group.append(row)
    if group:
        yield group


def allele_shape(row: Row, alleles: list[tuple[str, str]]) -> str:
    """The shape of a row's allele, `snv`, `mnv` or `indel`, from the input record's REF and
    ALT pair the allele came from: a pair of one length is a substitution, of two lengths an
    indel. A two-base substitution and a two-base insertion print the same Location and Allele
    (Ensembl VEP writes an insertion between its two flanking bases), so the row alone cannot
    tell them apart. The pairs are read as Ensembl VEP's parser leaves them (the rule of
    `record_names`): a bi-allelic indel loses its anchor base and the bases shared at either
    end, a multi-allelic record with an indel loses the first base of every allele when all
    share it and nothing more. A row whose allele matches no pair, or a record with no pairs,
    falls back to `variant_class`, which reads that shape as an insertion."""
    allele = row.allele.upper()
    pairs = [(r.upper(), a.upper()) for r, a in alleles if not (a.startswith("<") or "[" in a or "]" in a)]
    multi_indel = len(pairs) > 1 and any(len(r) != len(a) for r, a in pairs)
    shared_first = len({x[:1] for r, a in pairs for x in (r, a) if "*" not in x}) == 1
    for ref, alt in pairs:
        r, a = ref, alt
        if multi_indel and shared_first:
            # A `*` allele keeps its one character, as `record_names` keeps it.
            r, a = r[1:], (a if "*" in a else a[1:])
        elif not multi_indel and len(r) != len(a):
            if r[:1] == a[:1]:
                r, a = r[1:], a[1:]
            while r and a and r[0] == a[0]:
                r, a = r[1:], a[1:]
            while r and a and r[-1] == a[-1]:
                r, a = r[:-1], a[:-1]
        if (a or "-") != allele:
            continue
        if len(ref) != len(alt):
            return "indel"
        return "snv" if len(alt) == 1 else "mnv"
    klass = variant_class(row.location, row.allele)
    if klass == "snv":
        return "snv"
    if klass == "mnv_or_complex":
        return "mnv"
    return "indel"


def regulatory_classes(rows: list[Row], alleles: list[tuple[str, str]]) -> set[str]:
    """The `regulatory_class` classes a record's reference rows show, eight in all:
    `rf_only`, `mf_only` or `rf_and_mf` for the regulatory feature types among the rows'
    `Feature_type`; `regulatory_only` when the record has one of them and no Transcript row
    (its remaining row is the intergenic one Ensembl VEP writes for a variant no transcript
    overlaps); per MotifFeature row `motif_snv`, `motif_mnv` or `motif_indel` by the shape of
    its allele (`allele_shape`, over the record's REF and ALT pairs), and `motif_reverse` for a
    motif on the reverse strand. A record with no regulatory row shows no class."""
    types = {r.feature_type for r in rows}
    has_rf, has_mf = "RegulatoryFeature" in types, "MotifFeature" in types
    if not (has_rf or has_mf):
        return set()
    if has_rf and has_mf:
        classes = {"rf_and_mf"}
    elif has_rf:
        classes = {"rf_only"}
    else:
        classes = {"mf_only"}
    if "Transcript" not in types:
        classes.add("regulatory_only")
    for r in rows:
        if r.feature_type != "MotifFeature":
            continue
        if r.klass in SEQUENCE_CLASSES:
            classes.add(f"motif_{allele_shape(r, alleles)}")
        if r.strand == "-1":
            classes.add("motif_reverse")
    return classes


RECORD_CLASS_FUNCTIONS: dict[str, Callable[[list[Row], list[tuple[str, str]]], set[str]]] = {
    REGULATORY_CLASS: regulatory_classes,
}
RECORD_CLASS_KEYS = tuple(RECORD_CLASS_FUNCTIONS)


def input_alleles(path: Path) -> dict[tuple[str, str], list[tuple[str, str]]]:
    """Every (REF, ALT) pair of an input VCF's records by (chromosome, name), under every
    Uploaded_variation name the record can carry (`record_names`)."""
    out: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for is_header, f in vcf_lines(path):
        if is_header:
            continue
        pairs = [(f[3], alt) for alt in f[4].split(",")]
        for name in record_names("\t".join(f)):
            out.setdefault((f[0], name), []).extend(pairs)
    return out


def substitution_alts(path: Path) -> dict[tuple[str, str], list[tuple[int, str]]]:
    """The POS and bases of every ALT of an input VCF that is as long as its record's REF of
    two bases or more, a substitution of that many bases, by (chromosome, name), under every
    Uploaded_variation name the record can carry (`record_names`)."""
    out: dict[tuple[str, str], list[tuple[int, str]]] = {}
    for is_header, f in vcf_lines(path, 8):
        if is_header or len(f[3]) < 2:
            continue
        alts = [(int(f[1]), alt) for alt in f[4].split(",") if len(alt) == len(f[3]) and re.fullmatch(r"[A-Za-z]+", alt)]
        if alts:
            for name in record_names("\t".join(f)):
                out.setdefault((f[0], name), []).extend(alts)
    return out


def holds_allele(span: str, allele: str, alts: Iterable[tuple[int, str]]) -> bool:
    """Whether a substitution of `alts` (`substitution_alts`) holds `allele` at the positions
    of `span`, a row's Location after its chromosome: the row is then that substitution's,
    whole or as Ensembl VEP prints it once the bases its alleles share are dropped (the anchor
    base of a multi-allelic record with an indel, either end under `--minimal`). The record's
    position is part of the test, so a row never reads a record that only shares its name."""
    start, _, end = span.partition("-")
    try:
        a, b = int(start), int(end or start)
    except ValueError:
        return False
    return any(pos <= a and len(allele) == b - a + 1 and alt[a - pos:b - pos + 1].upper() == allele.upper() for pos, alt in alts)


def record_class_pass(
    suites: list[tuple[str, Path, Path]], seed: int, key: str
) -> tuple[dict[tuple[str, str, str], set[str]], dict[str, list[list[Row]]], Counter[str]]:
    """One streaming pass over the references for the record classes of `key`, each record
    classed from its rows and its input record's alleles: every classed record's classes by
    (suite, Uploaded_variation, Location), per class a bounded uniform sample of the records
    showing it (each as its rows), and the count of records per class."""
    rng = random.Random(seed)
    classes_of: dict[tuple[str, str, str], set[str]] = {}
    samples: dict[str, list[list[Row]]] = defaultdict(list)
    seen: Counter[str] = Counter()
    for sid, ref, inp in suites:
        alleles = input_alleles(inp)
        for rows in iter_records(ref, sid, substitution_alts(inp)):
            chrom = rows[0].location.split(":", 1)[0]
            classes = RECORD_CLASS_FUNCTIONS[key](rows, alleles.get((chrom, rows[0].uploaded_variation), []))
            if not classes:
                continue
            classes_of[(sid, *rows[0].record_key)] = classes
            for c in sorted(classes):
                seen[c] += 1
                reservoir_step(rng, samples[c], rows, seen[c], RECORD_CLASS_RESERVOIR)
    return classes_of, samples, seen


def parse_suite_arg(arg: str) -> tuple[str, Path, Path]:
    """`sid=reference.txt=input.vcf[.gz]`."""
    parts = arg.split("=")
    if len(parts) != 3:
        raise SystemExit(f"ERROR: [build_golden_corpus] --suite wants sid=reference=input, got {arg!r}")
    return parts[0], Path(parts[1]), Path(parts[2])


BND_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
SHARD_NAME = re.compile(r"^(\d+)-(\d+)\.json(\.gz)?$")


class CacheCost:
    """Compressed bytes a row adds to the fixture cache `prune_json_cache.py` cuts for it:
    every transcript of the source cache whose span, extended by `flank`, overlaps the row's
    location (and a bracket breakend's mate position), costed at its own deflate size. A
    shard is read on first touch; the source stores a transcript in every shard it reaches,
    so only the shards overlapping the flanked location are read."""

    def __init__(self, cache: Path, flank: int) -> None:
        self.flank = flank
        self.shards: dict[str, list[tuple[int, int, Path]]] = {}
        for chrom_dir in (cache / "transcripts").iterdir():
            self.shards[chrom_dir.name] = sorted(
                (int(m.group(1)), int(m.group(2)), shard)
                for shard in chrom_dir.iterdir()
                if (m := SHARD_NAME.match(shard.name))
            )
        self.loaded: dict[Path, list[tuple[int, int, int]]] = {}
        self.memo: dict[tuple[str, int, int], int] = {}

    def transcripts(self, shard: Path) -> list[tuple[int, int, int]]:
        if shard not in self.loaded:
            with opener(shard)(shard, "rt", encoding="utf-8") as handle:
                records = json.load(handle)
            self.loaded[shard] = [
                (int(t["start"]), int(t["end"]),
                 len(zlib.compress(json.dumps(t, separators=(",", ":"), ensure_ascii=False).encode("utf-8"))))
                for t in records
            ]
        return self.loaded[shard]

    def span_cost(self, chrom: str, start: int, end: int) -> int:
        key = (chrom, start, end)
        if key not in self.memo:
            lo, hi = start - self.flank, end + self.flank
            total = 0
            for shard_start, shard_end, shard in self.shards.get(chrom, ()):
                if shard_start <= hi and shard_end >= lo:
                    total += sum(b for t_start, t_end, b in self.transcripts(shard) if t_start <= hi and t_end >= lo)
            self.memo[key] = total
        return self.memo[key]

    def row_cost(self, row: Row) -> int:
        chrom, _, coords = row.location.partition(":")
        a, _, b = coords.partition("-")
        try:
            start, end = int(a), int(b or a)
        except ValueError:
            return 0
        start, end = min(start, end), max(start, end)
        total = self.span_cost(chrom, start, end)
        for m in BND_MATE.finditer(row.allele):
            mate_chrom, mate_pos = m.group(1).removeprefix("chr"), int(m.group(2))
            if mate_chrom == chrom:
                # A same-chromosome mate is annotated over the whole breakend-pair span.
                total = self.span_cost(chrom, min(start, mate_pos), max(end, mate_pos))
            else:
                total += self.span_cost(mate_chrom, mate_pos, mate_pos)
        return total


# ----------------------------------------------------------------------------- select
def reservoir_step(rng: random.Random, bucket: list, item, offered: int, capacity: int) -> None:
    """One reservoir step: `bucket` stays a uniform sample of at most `capacity` of the
    `offered` items seen so far for its key."""
    if len(bucket) < capacity:
        bucket.append(item)
    else:
        j = rng.randrange(offered)
        if j < capacity:
            bucket[j] = item


def reservoir_pass(
    suites: list[tuple[str, Path, Path]],
    seed: int,
    family_specs: list[str] | None = None,
    named: dict[str, set[str]] | None = None,
) -> tuple[dict[str, list[Row]], Counter[str], dict[str, Counter[str]], dict[str, list[Row]], dict[tuple[str, str], list[Row]]]:
    """One streaming pass: per consequence set, a bounded uniform sample of rows plus counts,
    per family stratum (when `family_specs` is given) a smaller sample of its own, so a
    stratum too rare to reach any consequence set's sample still has candidates, and every
    row of the records `named` (per suite, their Uploaded_variation names), keyed by suite
    and name, which the samples' bounds never touch."""
    rng = random.Random(seed)
    reservoir: dict[str, list[Row]] = defaultdict(list)
    seen: Counter[str] = Counter()
    per_suite: dict[str, Counter[str]] = defaultdict(Counter)
    families: dict[str, list[Row]] = defaultdict(list)
    family_seen: Counter[str] = Counter()
    named_rows: dict[tuple[str, str], list[Row]] = defaultdict(list)
    for sid, ref, inp in suites:
        wanted_names = (named or {}).get(sid, set())
        for row in iter_rows(ref, sid, family_specs, substitution_alts(inp)):
            c = row.consequence_set
            seen[c] += 1
            per_suite[sid][c] += 1
            reservoir_step(rng, reservoir[c], row, seen[c], RESERVOIR)
            if family_specs:
                stratum = row.family_stratum
                family_seen[stratum] += 1
                reservoir_step(rng, families[stratum], row, family_seen[stratum], FAMILY_RESERVOIR)
            if row.uploaded_variation in wanted_names:
                named_rows[(sid, row.uploaded_variation)].append(row)
    return reservoir, seen, per_suite, families, named_rows


def exemplar_order(row: Row, max_span: int, cost: Callable[[Row], int] | None) -> tuple[int, int, int]:
    """Sort key of a candidate row: no wider than `max_span` (when > 0) first, then the fewest
    cache bytes under `cost`, then the highest preference."""
    return (1 if max_span and row.span > max_span else 0, cost(row) if cost else 0, -row.preference())


def ensure_class_exemplars(
    chosen: dict[str, list[Row]],
    reservoir: dict[str, list[Row]],
    max_span: int,
    cost: Callable[[Row], int] | None,
) -> dict[str, list[Row]]:
    """Every variant class the references show keeps an exemplar: a class no chosen row has
    takes, from the whole reservoir, the row that sorts first in `exemplar_order`, added to the
    combination it was sampled for."""
    covered = {r.klass for rows in chosen.values() for r in rows}
    best: dict[str, tuple[tuple[int, int, int], str, Row]] = {}
    for combo, rows in reservoir.items():
        for r in rows:
            if r.klass in covered:
                continue
            key = exemplar_order(r, max_span, cost)
            if r.klass not in best or key < best[r.klass][0]:
                best[r.klass] = (key, combo, r)
    for _key, combo, r in best.values():
        chosen.setdefault(combo, []).append(r)
    return chosen


def ensure_family_exemplars(
    chosen: dict[str, list[Row]],
    families: dict[str, list[Row]],
    max_span: int,
    cost: Callable[[Row], int] | None,
) -> dict[str, list[Row]]:
    """Every family stratum the references show keeps an exemplar: a stratum no chosen row
    has takes, from its own sample, the row that sorts first in `exemplar_order`, added to
    the combination it was sampled for."""
    covered = {r.family_stratum for rows in chosen.values() for r in rows}
    for stratum, rows in sorted(families.items()):
        if stratum in covered:
            continue
        r = min(rows, key=lambda row: exemplar_order(row, max_span, cost))
        chosen.setdefault(r.consequence_set, []).append(r)
    return chosen


def choose_exemplars(
    reservoir: dict[str, list[Row]],
    k: int,
    rng: random.Random,
    max_span: int = 0,
    focus_terms: Iterable[str] = (),
    focus_k: int | None = None,
    cost: Callable[[Row], int] | None = None,
) -> dict[str, list[Row]]:
    """Up to K rows per combination, distinct (class, strand, family classes) first, in
    `exemplar_order`.

    Rows wider than `max_span` bases (when > 0) come after every narrower row, so a
    combination takes a wide structural variant only when nothing narrower shows it:
    each wide record pulls every transcript of its span into the fixture cache. With
    `cost` (bytes the row adds to the fixture cache, `CacheCost.row_cost`) the cheapest
    rows come next, then the higher preference; without it the preference alone orders
    the rows and ties fall in random order.

    A combination containing any of `focus_terms` takes up to `focus_k` rows instead of
    `k`, so a corpus can carry many exemplars of the consequence shapes it exists to test
    beside one exemplar of every other shape.
    """
    focus = set(focus_terms)
    chosen: dict[str, list[Row]] = {}
    for combo in sorted(reservoir):
        want = focus_k if focus_k is not None and focus & set(combo.split(",")) else k
        cands = list(reservoir[combo])
        rng.shuffle(cands)
        cands.sort(key=lambda r: exemplar_order(r, max_span, cost))
        picked: list[Row] = []
        strata_seen: set[tuple[str, str, tuple[str, ...]]] = set()
        for r in cands:
            stratum = (r.klass, r.strand, r.family)
            if stratum in strata_seen:
                continue
            picked.append(r)
            strata_seen.add(stratum)
            if len(picked) == want:
                break
        for r in cands:
            if len(picked) == want:
                break
            if r not in picked:
                picked.append(r)
        chosen[combo] = picked
    return chosen


def vep_auto_names(chrom: str, pos: int, ref: str, alts: list[str]) -> list[str]:
    """Reproduce VEP's auto-generated Uploaded_variation for an ID-less VCF record, in the
    two forms releases 116 and 115 write.

    Release 116 (ensembl-vep PR 1941) names the record from the raw line,
    `CHR_POS_REF/ALT1/ALT2`, and that name comes first. Release 115's names
    follow: the alleles are the RAW REF and ALTs joined by `/` (VEP's
    `nontrimmed_allele_string`) and the start is the variant's after parsing: unchanged
    when no allele differs in length from REF; for a bi-allelic indel the first shared
    base is chopped and the start then advances over the shared leading bases
    (`trim_sequences`, whose trailing trim moves only the end); for a multi-allelic
    indel the first base is chopped when every allele (ignoring `*`) shares it. Symbolic
    and breakend alleles: `chr_start_<ALT...>` with start POS or POS+1 (both are
    returned; which one VEP used depends on whether the REF base was padding).
    """
    alleles = [ref] + alts
    raw = "/".join(alleles)
    names = [f"{chrom}_{pos}_{raw}"]
    if any(a.startswith("<") or "[" in a or "]" in a for a in alts):
        joined = "/".join(alts)
        names += [f"{chrom}_{pos}_{joined}", f"{chrom}_{pos + 1}_{joined}"]
    elif len(alts) == 1:
        r, a, start = ref, alts[0], pos
        if len(r) != len(a):
            if r[:1] == a[:1]:
                r, a, start = r[1:] or "-", a[1:] or "-", start + 1
            while r and a and r[0] == a[0]:
                r, a, start = r[1:], a[1:], start + 1
            names.append(f"{chrom}_{start}_{raw}")
    elif any(len(a) != len(ref) for a in alts):
        firsts = {x[:1] for x in alleles if "*" not in x}
        if len(firsts) == 1:
            names.append(f"{chrom}_{pos + 1}_{raw}")
    return list(dict.fromkeys(names))


MULTI_ALLELIC_STRATA = ("multi_allelic_snv", "multi_allelic_indel", "multi_allelic_indel_trimmable")


def multi_allelic_stratum(ref: str, alts: list[str]) -> str:
    """Which multi-allelic stratum a record belongs to: SNV-only alleles; alleles of
    differing length; or such alleles that VEP's bi-allelic minimisation WOULD trim
    further after the shared first base is chopped (a shared prefix or suffix left
    between REF and some ALT), the shape on which a per-allele minimiser and VEP's
    untouched multi-allelic record disagree on positions and distances."""
    if not any(len(a) != len(ref) for a in alts):
        return "multi_allelic_snv"
    alleles = [ref] + alts
    if len({x[:1] for x in alleles}) == 1:
        alleles = [x[1:] or "-" for x in alleles]
    r = alleles[0]
    for a in alleles[1:]:
        if r == "-" or a == "-" or len(r) == len(a):
            continue
        if r[0] == a[0] or r[-1] == a[-1]:
            return "multi_allelic_indel_trimmable"
    return "multi_allelic_indel"


@dataclass
class ResolvedInputs:
    """What one scan of the input VCFs yields: per suite the record lines whose
    Uploaded_variation (ID, or VEP's auto-name) was wanted; the `#CHROM` header line the
    corpus writes (eight columns, or `FORMAT` and the kept samples after them); per suite
    and multi-allelic stratum a uniform sample of records; and, with kept samples, the names
    kept, how many records show each genotype class and a uniform sample of records per
    class."""
    found: dict[str, dict[tuple[str, str], list[str]]]
    header: str
    multi: dict[str, dict[str, list[str]]]
    samples: list[str]
    genotype_seen: Counter[str]
    genotype_pool: dict[str, list[tuple[str, str]]]


def resolve_records(
    inputs: list[tuple[str, Path]],
    wanted: dict[str, set[str]],
    multi_per_suite: int = 0,
    seed: int = 1,
    keep_samples: int = 0,
    genotype_per_class: int = 0,
) -> ResolvedInputs:
    """Scan each input VCF once. A record line is its first eight columns, or with
    `keep_samples` its first eight, `FORMAT` and the first `keep_samples` sample columns.
    The multi-allelic sample (`multi_per_suite` sequence records per suite and stratum:
    SNV-only alleles; alleles of differing length) and the genotype-class sample
    (`genotype_per_class` records per class the kept samples show) are records no
    consequence set can select for: the first exercises VEP's treatment of a multi-allelic
    record as one variant, the second the genotype shapes of `--individual`."""
    rng = random.Random(seed)
    found: dict[str, dict[tuple[str, str], list[str]]] = defaultdict(lambda: defaultdict(list))
    multi: dict[str, dict[str, list[str]]] = defaultdict(lambda: {k: [] for k in MULTI_ALLELIC_STRATA})
    multi_seen: Counter[tuple[str, str]] = Counter()
    header: str | None = None
    samples: list[str] = []
    genotype_seen: Counter[str] = Counter()
    genotype_pool: dict[str, list[tuple[str, str]]] = defaultdict(list)
    kept_columns = 9 + keep_samples if keep_samples else 8
    for sid, path in inputs:
        want = wanted.get(sid, set())
        if not want and not multi_per_suite and not keep_samples:
            continue
        named = False
        for is_header, f in vcf_lines(path, kept_columns):
            if is_header:
                if keep_samples:
                    samples = kept_samples(f, keep_samples, path, samples)
                    named = True
                if header is None:
                    header = "\t".join(f[:kept_columns])
                continue
            if keep_samples and not named:
                raise failure(f"{path} has no #CHROM header line naming its samples")
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            alts = alt.split(",")
            record = "\t".join(f[:kept_columns])
            if keep_samples:
                # Sorted so the draws fall on the classes in one order in every interpreter.
                for c in sorted(record_genotype_classes(f, keep_samples)):
                    genotype_seen[c] += 1
                    if genotype_per_class:
                        reservoir_step(rng, genotype_pool[c], (sid, record), genotype_seen[c], genotype_per_class)
            if multi_per_suite and len(alts) > 1 and not any(a.startswith("<") or "[" in a or "]" in a or a == "*" for a in alts):
                stratum = multi_allelic_stratum(ref, alts)
                multi_seen[(sid, stratum)] += 1
                reservoir_step(rng, multi[sid][stratum], record, multi_seen[(sid, stratum)], multi_per_suite)
            if not want:
                continue
            names = list(dict.fromkeys([*record_names(record), *vep_auto_names(chrom, pos, ref, alts)]))
            for n in names:
                if n in want:
                    found[sid][(n, chrom)].append(record)
                    break
    if keep_samples and not samples:
        raise failure("no input has a #CHROM header line naming its samples")
    return ResolvedInputs(found, header or "\t".join(FIXED_COLUMNS), multi, samples, genotype_seen, genotype_pool)


def run_options(args: argparse.Namespace) -> dict:
    """The manifest keys the golden harness passes to both engines: `flags` (extra
    command-line flags such as `--hgvs`) and `fasta` (a gzipped reference beside the
    manifest, decompressed and passed to `--fasta`)."""
    options: dict = {}
    if getattr(args, "flag", None):
        options["flags"] = list(args.flag)
    if getattr(args, "fasta_name", None):
        options["fasta"] = args.fasta_name
    return options


def failure(message: str) -> SystemExit:
    return SystemExit(f"ERROR: [build_golden_corpus] {message}")


def is_name_list(value: object) -> bool:
    return isinstance(value, list) and all(isinstance(n, str) and n.strip() for n in value)


def masks_shape_error(masks: object, owner: str) -> str | None:
    """The reason a `masks` list (`owner` names it: the manifest's own or a flag set's)
    cannot be applied: not a list of names, or a class other than `MASK_CLASSES`."""
    if not is_name_list(masks):
        return f"{owner} is not a list of names"
    unknown = sorted(set(masks) - set(MASK_CLASSES))
    if unknown:
        return f"{owner} names unknown classes {unknown} (the classes are {', '.join(MASK_CLASSES)})"
    return None


FLAG_SET_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]*")


def plugin_files_shape_error(value: object) -> str | None:
    """The reason a `plugin_files` map cannot be read: it must map each plugin name to a
    list of paths inside the corpus (relative, no `..` segment, no whitespace)."""
    if not isinstance(value, dict) or not all(isinstance(k, str) and k.strip() for k in value):
        return "plugin_files is not a map of plugin names"
    for plugin, files in value.items():
        if not is_name_list(files):
            return f"plugin_files for {plugin} is not a list of paths"
        bad = [f for f in files if not is_corpus_path(f) or any(c.isspace() for c in f)]
        if bad:
            return f"plugin_files for {plugin} names paths outside the corpus {bad}"
    return None


def is_flag_set_name(value: object) -> bool:
    return isinstance(value, str) and value != "default" and FLAG_SET_NAME.fullmatch(value) is not None


def is_corpus_path(value: object) -> bool:
    """A non-empty path relative to the corpus directory with no `..` segment."""
    if not isinstance(value, str) or not value:
        return False
    path = Path(value)
    return not path.is_absolute() and ".." not in path.parts


def names_format(flags: list[str]) -> bool:
    """Whether a flag list carries `--format` (as its own token or as `--format=<format>`)."""
    return any(f == "--format" or f.startswith("--format=") for f in flags)


def flag_set_shape_error(fs: object) -> str | None:
    """The reason one `flag_sets` entry cannot be run: the shapes `flag-set` writes, so a
    hand-written entry is held to the same rules; its hand-written `masks` (classes of
    `MASK_CLASSES`) and `reference_warnings` (non-empty texts) are held to the rules the
    manifest's own are, as the golden harness holds them."""
    name = fs.get("name") if isinstance(fs, dict) else None
    if not is_flag_set_name(name):
        return "has a flag_sets entry without a name (letters, digits, _ and -, not default)"
    if "flags" in fs and not (isinstance(fs["flags"], list) and all(isinstance(f, str) for f in fs["flags"])):
        return f"flag set {name} flags is not a list of strings"
    for key in ("input", "expected_dir"):
        if key in fs and not is_corpus_path(fs[key]):
            return f"flag set {name} {key} is not a path inside the corpus directory"
    for key in ("stdin", "unordered"):
        if key in fs and not isinstance(fs[key], bool):
            return f"flag set {name} {key} is not true or false"
    if "formats" in fs and not (isinstance(fs["formats"], list) and fs["formats"]
                                and all(f in OUTPUT_FORMATS for f in fs["formats"])):
        return f"flag set {name} formats is not a non-empty subset of {','.join(OUTPUT_FORMATS)}"
    if fs.get("stdin") and not names_format(fs.get("flags") or []):
        return f"flag set {name} runs on stdin without --format among its flags"
    if "perl_mode" in fs and fs["perl_mode"] not in PERL_MODES:
        return f"flag set {name} perl_mode is not one of {', '.join(PERL_MODES)}"
    for key, kind in (("covers", "names"), ("reference_warnings", "texts")):
        if key in fs and not is_name_list(fs[key]):
            return f"flag set {name} {key} is not a list of {kind}"
    if "masks" in fs:
        return masks_shape_error(fs["masks"], f"flag set {name} masks")
    return None


def with_reference_sets_error(entry: dict) -> str | None:
    """The reason a `divergences` entry's `WITH_REFERENCE_SETS` cannot be read, `None` for an entry
    without the key or with one of its shape: a non-empty list of distinct consequence sets in sorted
    order, each a sorted, comma-joined set of terms, other than the entry's `vep_rs_consequence_sets`
    as a whole; a list that holds Ensembl VEP's set (`vep_consequence_set`) holds nothing else. The
    format comparator and the golden test refuse the same shapes."""
    if WITH_REFERENCE_SETS not in entry:
        return None
    sets = entry[WITH_REFERENCE_SETS]
    if not isinstance(sets, list) or not sets:
        return "is not a non-empty list"
    if not all(isinstance(s, str) and s and normalize_consequence_set(s) == s for s in sets):
        return "holds an item that is not a sorted, comma-joined set of terms"
    if sets != sorted(set(sets)):
        return "is not a sorted list of distinct sets"
    if sets == entry.get("vep_rs_consequence_sets"):
        return "repeats vep_rs_consequence_sets and states no second outcome"
    vep_set = entry.get("vep_consequence_set")
    if len(sets) > 1 and isinstance(vep_set, str) and normalize_consequence_set(vep_set) in sets:
        return "names Ensembl VEP's set beside another"
    return None


def entry_key_text(entry: dict) -> str:
    """How an error names a classify entry: its Location, Allele and Feature."""
    return " ".join(str(entry.get(name, "")) for name in ("location", "allele", "feature"))


def classify_lists_shape_error(run: dict, owner: str) -> str | None:
    """The reason the classify lists of one run (the manifest's own or a flag set's, named `owner`)
    cannot be read: a `divergences` entry whose `WITH_REFERENCE_SETS` is of another shape
    (`with_reference_sets_error`), or the key on an entry of `vep_rs_only_tuples` or
    `field_divergences`, which state no Ensembl VEP set to hold it against."""
    for section in ("divergences", "vep_rs_only_tuples", "field_divergences"):
        entries = run.get(section)
        for entry in entries if isinstance(entries, list) else []:
            if not isinstance(entry, dict) or WITH_REFERENCE_SETS not in entry:
                continue
            reason = with_reference_sets_error(entry) if section == "divergences" else "belongs on a divergences entry"
            if reason:
                return f"{owner} {section} entry on {entry_key_text(entry)}: {WITH_REFERENCE_SETS} {reason}"
    return None


def manifest_shape_error(manifest: dict) -> str | None:
    """The reason a manifest's `covers`, `masks`, `perl_mode`, `reference_warnings`, `plugin_files`,
    `flag_sets`, `cache_layout`, `fasta_in_cache`, `included_records`, `record_class_key` (one of
    `RECORD_CLASS_KEYS` when present) or sample keys cannot be read:
    `covers` and `reference_warnings` must be lists of non-empty strings, `masks` a list of
    `MASK_CLASSES` (`masks_shape_error`), `perl_mode` one of `PERL_MODES`, `plugin_files` a map
    of plugin names to paths inside the corpus
    (`plugin_files_shape_error`), `flag_sets` a
    list of entries of the shape `flag_set_shape_error` holds them to, with distinct names,
    `cache_layout`, when present, `ensembl` (the cache at `cache_root/<species>/<release>_<assembly>`,
    so `release` and `assembly` must be named) rather than the default `json_cache/`, `species`,
    when present, an Ensembl species name (`homo_sapiens` when absent), `fasta_in_cache`,
    when present, a boolean that is true only beside a `fasta` (the golden run then finds the reference
    inside its cache directory instead of receiving `--fasta`), `keep_samples` a count,
    `samples` the names of that many kept samples when the count is above zero, and
    `included_records`, when present, a map of non-empty suite names to non-empty lists of
    record names. A `vep_rs_consequence_sets_with_reference_sequence` key in a classify list of the
    manifest or of a flag set has the shape `classify_lists_shape_error` holds it to."""
    if "covers" in manifest and not is_name_list(manifest["covers"]):
        return "covers is not a list of names"
    if "perl_mode" in manifest and manifest["perl_mode"] not in PERL_MODES:
        return f"perl_mode is not one of {', '.join(PERL_MODES)}"
    if "masks" in manifest:
        reason = masks_shape_error(manifest["masks"], "masks")
        if reason:
            return reason
    if "reference_warnings" in manifest and not is_name_list(manifest["reference_warnings"]):
        return "reference_warnings is not a list of texts"
    if "plugin_files" in manifest:
        reason = plugin_files_shape_error(manifest["plugin_files"])
        if reason:
            return reason
    keep = manifest.get("keep_samples", 0)
    if not isinstance(keep, int) or isinstance(keep, bool) or keep < 0:
        return "keep_samples is not a count"
    if keep and not (is_name_list(manifest.get("samples")) and len(manifest["samples"]) == keep):
        return f"samples does not name {keep} kept samples"
    required = manifest.get("required_genotype_classes", [])
    if not (isinstance(required, list) and all(c in GENOTYPE_CLASSES for c in required)):
        return "required_genotype_classes is not a list of genotype classes"
    if "record_class_key" in manifest and manifest["record_class_key"] not in RECORD_CLASS_KEYS:
        return f"record_class_key is not one of {', '.join(RECORD_CLASS_KEYS)}"
    flag_sets = manifest.get("flag_sets", [])
    if not isinstance(flag_sets, list):
        return "flag_sets is not a list"
    names: set[str] = set()
    for fs in flag_sets:
        reason = flag_set_shape_error(fs)
        if reason:
            return reason
        if fs["name"] in names:
            return f"names flag set {fs['name']} twice"
        names.add(fs["name"])
    if "species" in manifest and not (isinstance(manifest["species"], str) and SPECIES_RE.fullmatch(manifest["species"])):
        return f"species is {json.dumps(manifest['species'])}, not an Ensembl species name"
    if "cache_layout" in manifest and manifest["cache_layout"] != "ensembl":
        return f'cache_layout is {json.dumps(manifest["cache_layout"])}, not "ensembl"'
    if manifest.get("cache_layout") == "ensembl" and not all(isinstance(manifest.get(k), str) and manifest[k] for k in ("release", "assembly")):
        return "cache_layout ensembl wants release and assembly"
    included = manifest.get("included_records", {})
    if not (isinstance(included, dict)
            and all(isinstance(sid, str) and sid.strip() and entry and is_name_list(entry) for sid, entry in included.items())):
        return "included_records is not a map of suite to record names"
    if "fasta_in_cache" in manifest and not isinstance(manifest["fasta_in_cache"], bool):
        return f'fasta_in_cache is {json.dumps(manifest["fasta_in_cache"])}, not a boolean'
    if manifest.get("fasta_in_cache") and not (isinstance(manifest.get("fasta"), str) and manifest["fasta"]):
        return "fasta_in_cache wants fasta"
    for owner, run in (("the manifest's", manifest), *((f"flag set {fs['name']}", fs) for fs in flag_sets)):
        reason = classify_lists_shape_error(run, owner)
        if reason:
            return reason
    return None


def read_manifest(corpus: Path) -> dict:
    """The corpus manifest, refused unless it is a JSON object whose `covers`, `flag_sets` and
    `cache_layout` have the shapes `manifest_shape_error` describes."""
    path = corpus / "manifest.json"
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise failure(f"{path} is not readable JSON ({exc})") from None
    if not isinstance(manifest, dict):
        raise failure(f"{path} is not a JSON object")
    reason = manifest_shape_error(manifest)
    if reason:
        raise failure(f"{path} {reason}")
    return manifest


def write_manifest(corpus: Path, manifest: dict) -> None:
    (corpus / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")


def matrix_items(path: Path) -> set[str]:
    """The `item` column of a parity matrix (tab-separated, `item` first)."""
    if not path.is_file():
        raise failure(f"{path} is not a file: name the parity matrix with --matrix")
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise failure(f"{path} is not a readable text file ({exc})") from None
    if not lines or lines[0].split("\t", 1)[0] != "item":
        raise failure(f"{path} is not a parity matrix (first column must be item)")
    return {line.split("\t", 1)[0] for line in lines[1:] if line.strip()}


def covers_owner(flag_set: str | None) -> str:
    """How an error names a `covers` list: the manifest's own, or that of the flag set `flag_set`."""
    return f"flag set {flag_set} covers" if flag_set is not None else "covers"


def covers_lists(manifest: dict) -> list[tuple[str, list]]:
    """Every `covers` list of a manifest with its owner: the manifest's own, then each flag set's."""
    lists = [(covers_owner(None), manifest.get("covers") or [])]
    for fs in manifest.get("flag_sets") or []:
        lists.append((covers_owner(fs["name"]), fs.get("covers") or []))
    return lists


def names_outside_matrix(names: Iterable[str], items: set[str], matrix: Path, owner: str) -> list[str]:
    """One error line per name that is not a row of the matrix, naming the list it came from."""
    return [f"{owner} names {n}, which is not a row of {matrix}" for n in names if n not in items]


def covers_names(args_names: list[str] | None) -> list[str]:
    """The `--covers` names of a command line, each non-empty, in order and without repeats."""
    names = list(dict.fromkeys(args_names or []))
    if not is_name_list(names):
        raise failure("--covers wants a non-empty row name")
    return names


def flag_set_named(manifest: dict, name: str) -> dict:
    for fs in manifest.get("flag_sets") or []:
        if fs["name"] == name:
            return fs
    raise failure(f"the manifest has no flag set {name!r}")


def path_inside_corpus(value: str, option: str) -> str:
    """A file or directory name relative to the corpus directory: no absolute path, no `..` segment."""
    if not is_corpus_path(value):
        raise failure(f"{option} wants a path inside the corpus directory, got {value!r}")
    return value


@dataclass
class NamedRecord:
    """A record the corpus carries outside the sampling: its suite, the label an error names it
    by (the `--include` argument, a kept record's first five columns, or a carried
    `included_records` name) and the Uploaded_variation names its reference rows can carry."""
    suite: str
    label: str
    names: tuple[str, ...]


def named_records(args: argparse.Namespace, existing: dict, suite_ids: set[str]) -> list[NamedRecord]:
    """The records `--include` names and, under `--keep-records`, every record of the manifest
    being regenerated and every name of its `included_records`, so a carried name is held to the
    checks a `--include` gets; a record of a suite that is not among `--suite` is refused."""
    named: list[NamedRecord] = []
    for arg in args.include or []:
        sid, sep, name = arg.partition("=")
        if not (sep and sid and name):
            raise failure(f"--include wants sid=name, got {arg!r}")
        named.append(NamedRecord(sid, f"--include {arg}", (name,)))
    if args.keep_records:
        if not existing:
            raise failure(f"--keep-records wants a manifest to regenerate in {args.out}")
        for rec in existing.get("records") or []:
            label = "kept record " + " ".join(rec["vcf"].split("\t")[:5])
            named.append(NamedRecord(rec["suite"], label, tuple(record_names(rec["vcf"]))))
        for sid, names in (existing.get("included_records") or {}).items():
            named.extend(NamedRecord(sid, f"carried include {sid}={name}", (name,)) for name in names)
    for n in named:
        if n.suite not in suite_ids:
            raise failure(f"{n.label} names suite {n.suite!r}, which is not among --suite")
    return named


def included_by_suite(includes: list[str], carried: dict[str, list[str]] | None = None) -> dict[str, list[str]]:
    """The `--include` arguments as suite -> sorted names, merged with the `included_records` a
    regeneration under `--keep-records` carries from the manifest it keeps: the new manifest's
    `included_records`."""
    by_suite: dict[str, set[str]] = defaultdict(set)
    for sid, names in (carried or {}).items():
        by_suite[sid].update(names)
    for arg in includes:
        sid, _sep, name = arg.partition("=")
        by_suite[sid].add(name)
    return {sid: sorted(names) for sid, names in sorted(by_suite.items())}


def cmd_select(args: argparse.Namespace) -> int:
    suites = [parse_suite_arg(s) for s in args.suite]
    out = Path(args.out)
    existing = read_manifest(out) if (out / "manifest.json").is_file() else {}
    species = args.species if args.species is not None else (existing.get("species") or DEFAULT_SPECIES)
    if not SPECIES_RE.fullmatch(species):
        print(f"ERROR: [build_golden_corpus] --species {species!r} is not an Ensembl species name "
              "(lower-case letters, digits and underscores)", file=sys.stderr)
        return 1
    named = named_records(args, existing, {sid for sid, _ref, _inp in suites})
    covers = list(dict.fromkeys([*(existing.get("covers") or []), *covers_names(args.covers)]))
    masks = existing.get("masks") or []
    perl_mode = existing.get("perl_mode")
    reference_warnings = existing.get("reference_warnings") or []
    plugin_files = existing.get("plugin_files") or {}
    # A set's divergence keys name the old records by ordinal; `classify --flag-set` rebuilds them.
    flag_sets = [{k: v for k, v in fs.items() if k not in CLASSIFY_KEYS} for fs in existing.get("flag_sets") or []]
    if not covers:
        print("ERROR: [build_golden_corpus] select writes no manifest without --covers: "
              "name the parity-matrix rows the corpus exercises", file=sys.stderr)
        return 1
    if args.matrix is not None:
        matrix = Path(args.matrix)
        errors = names_outside_matrix(covers, matrix_items(matrix), matrix, covers_owner(None))
        for e in errors:
            print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
        if errors:
            return 1
    family_specs = parse_family_columns(args.family_columns)
    keep_samples = args.keep_samples
    required = genotype_class_names(args.require_genotype_class)
    if keep_samples < 0:
        raise failure(f"--keep-samples wants a count of sample columns, got {keep_samples}")
    if args.genotype_exemplars < 0:
        raise failure(f"--genotype-exemplars wants a count of records per class, got {args.genotype_exemplars}")
    if required and not keep_samples:
        raise failure("--require-genotype-class wants --keep-samples: a corpus without sample columns shows no genotype")
    class_key = record_class_key(args.record_classes)
    if args.record_class_exemplars < 0:
        raise failure(f"--record-class-exemplars wants a count of records per class, got {args.record_class_exemplars}")
    if args.cache and not (Path(args.cache) / "transcripts").is_dir():
        raise failure(f"{args.cache} has no transcripts/ directory")
    names_by_suite: dict[str, set[str]] = defaultdict(set)
    for n in named:
        names_by_suite[n.suite].update(n.names)
    rng = random.Random(args.seed)
    reservoir, seen, per_suite, families, named_rows = reservoir_pass(suites, args.seed, family_specs, names_by_suite)
    unrowed = [n for n in named if not any((n.suite, name) in named_rows for name in n.names)]
    for n in unrowed:
        print(f"ERROR: [build_golden_corpus] {n.label} has no row in the {n.suite} reference", file=sys.stderr)
    if unrowed:
        return 1
    cost = CacheCost(Path(args.cache), args.flank).row_cost if args.cache else None
    chosen = choose_exemplars(reservoir, args.k, rng, args.max_span, args.focus_term or (), args.focus_k, cost)
    chosen = ensure_class_exemplars(chosen, reservoir, args.max_span, cost)
    if family_specs:
        chosen = ensure_family_exemplars(chosen, families, args.max_span, cost)
    for rows in named_rows.values():
        for r in rows:
            bucket = chosen.setdefault(r.consequence_set, [])
            if r not in bucket:
                bucket.append(r)
    classes_of: dict[tuple[str, str, str], set[str]] = {}
    class_seen: Counter[str] = Counter()
    if class_key:
        # Every record class the references show gets its exemplars: a class short of them
        # takes, from its own sample, the records that sort first in `exemplar_order`, each
        # counted towards every class it shows.
        classes_of, class_samples, class_seen = record_class_pass(suites, args.seed, class_key)
        chosen_records = {(r.suite, *r.record_key) for rows in chosen.values() for r in rows}
        showing = Counter(c for record in chosen_records for c in classes_of.get(record, ()))
        for c in sorted(class_seen):
            candidates = sorted(class_samples[c], key=lambda rows: min(exemplar_order(r, args.max_span, cost) for r in rows))
            for rows in candidates:
                if showing[c] >= args.record_class_exemplars:
                    break
                record = (rows[0].suite, *rows[0].record_key)
                if record in chosen_records:
                    continue
                chosen_records.add(record)
                showing.update(classes_of[record])
                for r in rows:
                    chosen.setdefault(r.consequence_set, []).append(r)
    wanted: dict[str, set[str]] = defaultdict(set)
    for combo, rows in chosen.items():
        for r in rows:
            wanted[r.suite].add(r.uploaded_variation)
    resolved = resolve_records(
        [(sid, inp) for sid, _ref, inp in suites], wanted, args.multi_allelic_per_suite, args.seed,
        keep_samples, args.genotype_exemplars,
    )
    records: dict[str, dict] = {}
    missing: list[str] = []
    for combo, rows in chosen.items():
        for r in rows:
            chrom = r.location.split(":", 1)[0]
            lines = resolved.found.get(r.suite, {}).get((r.uploaded_variation, chrom), [])
            if not lines:
                missing.append(f"{r.suite}\t{r.uploaded_variation}\t{r.location}")
                continue
            for line in lines:
                key = f"{r.suite}\t{line}"
                rec = records.setdefault(key, {"suite": r.suite, "vcf": line, "combinations": [], "exemplar_rows": []})
                if combo not in rec["combinations"]:
                    rec["combinations"].append(combo)
                rec["exemplar_rows"].append({
                    "uploaded_variation": r.uploaded_variation, "location": r.location, "allele": r.allele,
                    "feature": r.feature, "variant_class": r.klass, "strand": r.strand,
                    **({"family_stratum": r.family_stratum} if family_specs else {}),
                })
    for sid in sorted(resolved.multi):
        for stratum, lines in sorted(resolved.multi[sid].items()):
            for line in lines:
                key = f"{sid}\t{line}"
                rec = records.setdefault(key, {"suite": sid, "vcf": line, "combinations": [], "exemplar_rows": []})
                rec.setdefault("strata", []).append(stratum)
    if keep_samples:
        for rec in records.values():
            rec["genotype_classes"] = sorted(record_genotype_classes(rec["vcf"].split("\t"), keep_samples))
        showing = Counter(c for rec in records.values() for c in rec["genotype_classes"])
        for c in sorted(resolved.genotype_seen):
            for sid, line in resolved.genotype_pool[c]:
                if showing[c] >= args.genotype_exemplars:
                    break
                key = f"{sid}\t{line}"
                if key not in records:
                    classes = sorted(record_genotype_classes(line.split("\t"), keep_samples))
                    records[key] = {"suite": sid, "vcf": line, "combinations": [], "exemplar_rows": [], "genotype_classes": classes}
                    showing.update(classes)
    if class_key:
        # By the record's names rather than its exemplar rows: a record the multi-allelic or the
        # genotype sampling carried has none, and shows its classes all the same.
        classes_by_name: dict[tuple[str, str, str], set[str]] = defaultdict(set)
        for (sid, name, location), classes in classes_of.items():
            chrom = location.split(":", 1)[0]
            classes_by_name[(sid, name, chrom)] |= classes
        for rec in records.values():
            chrom = rec["vcf"].split("\t", 1)[0]
            rec["record_classes"] = sorted({
                c for name in record_names(rec["vcf"]) for c in classes_by_name.get((rec["suite"], name, chrom), ())
            })
    if len(records) > args.max_records:
        # The first record of every combination, family stratum, genotype class and record
        # class and every named record survive the cap, which those alone may exceed; the rest
        # fill in sorted record order, so the cut is deterministic.
        by_group: dict[str, list[str]] = defaultdict(list)
        for key, rec in records.items():
            for c in rec["combinations"]:
                by_group[f"combination {c}"].append(key)
            for row in rec["exemplar_rows"]:
                if "family_stratum" in row:
                    by_group[f"stratum {row['family_stratum']}"].append(key)
            for c in rec.get("genotype_classes", ()):
                by_group[f"genotype class {c}"].append(key)
            for c in rec.get("record_classes", ()):
                by_group[f"record class {c}"].append(key)
        keep = {keys[0] for keys in by_group.values()}
        for key, rec in records.items():
            suite_names = names_by_suite.get(rec["suite"], ())
            if any(row["uploaded_variation"] in suite_names for row in rec["exemplar_rows"]):
                keep.add(key)
        for key in sorted(records):
            if len(keep) >= args.max_records:
                break
            keep.add(key)
        records = {k: v for k, v in records.items() if k in keep}

    def sort_key(item):
        f = item[1]["vcf"].split("\t")
        chrom = f[0]
        chrom_rank = (0, int(chrom)) if chrom.isdigit() else (1, {"X": 1, "Y": 2, "MT": 3}.get(chrom, 9))
        return (chrom_rank, int(f[1]), f[3], f[4])

    ordered = sorted(records.items(), key=sort_key)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "variants.vcf").open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(resolved.header + "\n")
        for _key, rec in ordered:
            fh.write(rec["vcf"] + "\n")
    covered = {c for _k, rec in ordered for c in rec["combinations"]}
    genotype_classes = Counter(c for _k, rec in ordered for c in rec.get("genotype_classes", ()))
    uncovered_genotype_classes = sorted((set(resolved.genotype_seen) | set(required)) - set(genotype_classes))
    record_classes = Counter(c for _k, rec in ordered for c in rec.get("record_classes", ()))
    uncovered_record_classes = sorted(set(class_seen) - set(record_classes))
    carried = existing.get("included_records") if args.keep_records else None
    included = included_by_suite(args.include or [], carried)
    manifest = {
        "release": args.release,
        "species": species,
        "assembly": args.assembly,
        "generator_version": GENERATOR_VERSION,
        "seed": args.seed,
        "exemplars_per_combination": args.k,
        **({"focus_terms": sorted(args.focus_term), "focus_exemplars_per_combination": args.focus_k} if args.focus_term else {}),
        **({"included_records": included} if included else {}),
        **({"kept_records": len(existing.get("records") or [])} if args.keep_records else {}),
        "max_span": args.max_span,
        **({"cache_cost_flank": args.flank} if args.cache else {}),
        **({"keep_samples": keep_samples, "samples": resolved.samples, "genotype_exemplars_per_class": args.genotype_exemplars,
            "required_genotype_classes": required} if keep_samples else {}),
        **({"record_class_key": class_key, "record_class_exemplars": args.record_class_exemplars} if class_key else {}),
        **run_options(args),
        "covers": covers,
        **({"masks": masks} if masks else {}),
        **({"perl_mode": perl_mode} if perl_mode else {}),
        **({"reference_warnings": reference_warnings} if reference_warnings else {}),
        **({"plugin_files": plugin_files} if plugin_files else {}),
        **({"flag_sets": flag_sets} if flag_sets else {}),
        # File names only: the directory each suite was staged in belongs to the build
        # host, not to the corpus, and the tests read neither field.
        "suites": {sid: {"reference": Path(ref).name, "input": Path(inp).name, "rows": sum(per_suite[sid].values()), "combinations": len(per_suite[sid])} for sid, ref, inp in suites},
        "combinations": {c: {"observed_rows": seen[c], "exemplar_records": sum(1 for _k, rec in ordered if c in rec["combinations"])} for c in sorted(seen)},
        "variant_classes": dict(sorted(Counter(r.klass for rows in chosen.values() for r in rows).items())),
        **({"family_columns": family_specs,
            "family_strata": dict(sorted(Counter(r.family_stratum for rows in chosen.values() for r in rows).items()))}
           if family_specs else {}),
        **({"genotype_classes": dict(sorted(genotype_classes.items()))} if keep_samples else {}),
        **({"record_classes": dict(sorted(record_classes.items()))} if class_key else {}),
        "records": [
            {"suite": rec["suite"], "vcf": rec["vcf"], "combinations": sorted(rec["combinations"]),
             **({"strata": rec["strata"]} if rec.get("strata") else {}),
             **({"genotype_classes": rec["genotype_classes"]} if keep_samples else {}),
             **({"record_classes": rec["record_classes"]} if class_key else {})}
            for _k, rec in ordered
        ],
        "multi_allelic_per_suite": args.multi_allelic_per_suite,
        "uncovered_combinations": sorted(set(seen) - covered),
        **({"uncovered_genotype_classes": uncovered_genotype_classes} if keep_samples else {}),
        **({"uncovered_record_classes": uncovered_record_classes} if class_key else {}),
        "unresolved_exemplars": missing,
    }
    write_manifest(out, manifest)
    classes = manifest["variant_classes"]
    class_counts = ", ".join(f"{k} {v}" for k, v in classes.items())
    genotype_note = (f", genotype classes observed {len(resolved.genotype_seen)}, covered {len(genotype_classes)}"
                     if keep_samples else "")
    record_class_note = (f", record classes observed {len(class_seen)}, covered {len(record_classes)}" if class_key else "")
    # A record both included and kept is two entries of `named` and one key of `named_rows`.
    named_note = f", named records {len(named_rows)}" if named else ""
    print(f"combinations observed {len(seen)}, covered {len(covered)}, records {len(ordered)}, unresolved exemplars {len(missing)}, "
          f"variant classes {len(classes)} ({class_counts}){genotype_note}{record_class_note}{named_note}")
    return 0 if not (missing or manifest["uncovered_combinations"] or uncovered_genotype_classes or uncovered_record_classes) else 1


# --------------------------------------------------------------------------- classify
def run_masks(manifest: dict, target: dict) -> list[str]:
    """The mask classes of a run, as the golden harness reads them: the flag set's own `masks`
    when its entry (`target`) carries the key, else the manifest's; `target` is the manifest
    itself for the default run."""
    if "masks" in target:
        return list(target["masks"])
    return list(manifest.get("masks") or [])


def masked_feature_types(masks: Iterable[str]) -> frozenset[str]:
    """The `Feature_type` values a run's `masks` leave out of the reference: the regulatory
    rows' under `regulatory_feature_rows`; `order_only` leaves out no row (`fields_agree`
    gives it its meaning)."""
    return REGULATORY_FEATURE_TYPES if REGULATORY_FEATURE_ROWS in masks else frozenset()


def fields_agree(name: str, vep: str, vep_rs: str, order_only: bool) -> bool:
    """Whether two values of the Extra field `name` agree: exactly, or, for a column of
    `ORDER_ONLY_COLUMNS` under `order_only`, as a multiset of their comma-separated items;
    an absent value (`-`) agrees with itself alone."""
    if order_only and name in ORDER_ONLY_COLUMNS and "-" not in (vep, vep_rs):
        return Counter(vep.split(",")) == Counter(vep_rs.split(","))
    return vep == vep_rs


def load_tuples(path: Path, masked: frozenset[str] = frozenset()) -> tuple[dict[tuple[str, str, str, str], set[str]], dict[tuple[str, str, str, str], set[str]]]:
    """(Location, Allele, Feature, Feature_type) -> consequence-set strings, and -> Uploaded_variation
    names, leaving out the rows whose Feature_type is in `masked`."""
    out: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    names: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for r in iter_rows(path, "x"):
        if r.feature_type in masked:
            continue
        key = (r.location, r.allele, r.feature, r.feature_type)
        out[key].add(r.consequence_set)
        names[key].add(r.uploaded_variation)
    return out, names


def load_fields(path: Path, fields: list[str], masked: frozenset[str] = frozenset()) -> dict[tuple[str, str, str, str, str], dict[str, str]]:
    """(Location, Allele, Feature, Feature_type, consequence set) -> the named Extra fields ('-' when absent)
    and the row's `cDNA_position` and `CDS_position` columns, leaving out the rows whose Feature_type is
    in `masked`."""
    out: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            while len(f) < 14:
                f.append("-")
            if f[5] in masked:
                continue
            extra = parse_extra(f[13])
            key = (f[1], f[2], f[4], f[5], normalize_consequence_set(f[6]))
            row = {name: extra.get(name, "-") for name in fields}
            row.setdefault("cDNA_position", f[7])
            row.setdefault("CDS_position", f[8])
            out.setdefault(key, row)
    return out



def record_names(vcf_line: str) -> list[str]:
    """The Uploaded_variation names a VCF record line can carry: the first `;`-joined element
    of its ID, else VEP's auto-names when that element is `.`, empty or `0` (Parser/VCF.pm:347
    takes the first element; Parser.pm:462 replaces a name Perl reads as false,
    OutputFactory.pm:900 a `.`). An empty or `0` element also offers the name Parser.pm:462
    builds, `CHR_START_ALLELES` from the alleles the parser left (the anchor chop of
    Parser/VCF.pm:300-336, then `trim_sequences` on a bi-allelic indel) in the case written,
    which the output keeps when upper-casing changed the string (OutputFactory.pm:900)."""
    f = vcf_line.split("\t")
    chrom, pos, vid, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
    first_id = vid.split(";")[0]
    if first_id not in (".", "", "0"):
        return [first_id]
    alts = alt.split(",")
    names = vep_auto_names(chrom, pos, ref, alts)
    if first_id != "." and not any(a.startswith("<") or "[" in a or "]" in a for a in alts):
        start = pos
        if len(alts) == 1 and len(ref) != len(alts[0]):
            a = alts[0]
            if ref[:1] == a[:1]:
                ref, a, start = ref[1:] or "-", a[1:] or "-", start + 1
            if "-" not in (ref, a):
                while ref and a and ref[0] == a[0]:
                    ref, a, start = ref[1:], a[1:], start + 1
                while ref and a and ref[-1] == a[-1]:
                    ref, a = ref[:-1], a[:-1]
                ref, a = ref or "-", a or "-"
            alts = [a]
        elif len(alts) > 1 and any(len(a) != len(ref) for a in alts):
            if len({a[:1] for a in [ref, *alts] if "*" not in a}) == 1:
                ref = ref[1:] or "-"
                alts = [a if "*" in a else a[1:] or "-" for a in alts]
                start += 1
        names.append(f"{chrom}_{start}_{ref}/{'/'.join(alts)}")
    return list(dict.fromkeys(names))


def record_name_index(records: list[dict]) -> dict[str, list[int]]:
    """Uploaded_variation (the ID, else VEP's auto-name) -> ordinals of the corpus records it names."""
    index: dict[str, list[int]] = defaultdict(list)
    for i, rec in enumerate(records):
        for n in record_names(rec["vcf"]):
            index[n].append(i)
    return index


REFERENCE_SKIPPED = "reference_skipped_record"
UNEXPLAINED = "unexplained_residual"
# The classes a field divergence can be explained by (`classify_field_divergence`), each a page
# of the documented divergences.
MOTIF_SCORE_REVERSE_STRAND_MNV = "motif_score_reverse_strand_mnv"
CELL_TYPE_COMMA_NAMES = "cell_type_comma_names"
REFSEQ_EDITED_MODEL_HGVSC_OFFSET = "refseq_edited_model_hgvsc_offset"
REFSEQ_EDITED_MODEL_HGVS_OFFSET = "refseq_edited_model_hgvs_offset"
ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE = "organelle_peptide_paths_translation_table"
HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE = "hgvs_shifted_indel_written_span_gate"
HGVS_SHIFTED_DELETION_CLAMPED_SPAN = "hgvs_shifted_deletion_clamped_span"
HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE = "hgvs_clipped_reference_against_whole_alternate"
HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT = "hgvs_added_bases_written_outside_the_transcript"
HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE = "hgvsp_stop_named_on_incomplete_coding_sequence"
# The class of a SpliceAI field on a row whose score record holds several comma-joined entries:
# Ensembl VEP files the record under its first entry's symbol with every entry as one value
# (SpliceAI.pm:326-332, :305), where vep-rs reads each entry under its own symbol. The format
# comparator's predicate decides it from the score files `classify --spliceai-table` names.
SPLICEAI_MULTI_ENTRY_CELL = "spliceai_multi_entry_cell"
# The fields of a row that predicate reads, which ride along with the compared fields.
SPLICEAI_ROW_FIELDS = (
    "SpliceAI_pred", "SpliceAI_pred_SYMBOL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL",
    "SpliceAI_cutoff", "SYMBOL",
)
# The nine fields SpliceAI.pm:276-287 cut from a value for the VCF format, whatever the run's
# options, and the key of the object the JSON format writes them under.
SPLICEAI_MEMBER_FIELDS = SPLICEAI_ROW_FIELDS[1:10]
SPLICEAI_JSON_FIELD = "spliceai"
SPLICEAI_NUMBER = re.compile(r"-?\d+(\.\d+)?")
# A description that names the stop alone: unchanged, deleted, replaced by one residue or extended.
STOP_NAMED_ALONE = re.compile(r"Ter([1-9]\d*)(=|del|[A-Z][a-z]{2}(?:extTer(?:[1-9]\d*|\?))?)")
HGVSP_STOP_CODON_PROTEIN_COMPARISON = "hgvsp_stop_codon_protein_comparison"
# The eight descriptions the class reads: the stop unchanged, an extension named at it, the
# stop alone replaced, residues inserted between the last residue and the stop, a range that
# ends at the stop replaced, a frameshift, the stop deleted before an extension's count, and
# one residue and a stop inserted between two residues.
STOP_UNCHANGED = re.compile(r"Ter([1-9]\d*)=")
STOP_EXTENSION = re.compile(r"Ter([1-9]\d*)([A-Z][a-z]{2})extTer([1-9]\d*|\?)")
STOP_REPLACED = re.compile(r"Ter([1-9]\d*)delins((?:[A-Z][a-z]{2})+)")
STOP_INSERTION = re.compile(r"[A-Z][a-z]{2}([1-9]\d*)_Ter([1-9]\d*)ins((?:[A-Z][a-z]{2})+)")
STOP_RANGE_REPLACED = re.compile(r"([A-Z][a-z]{2})([1-9]\d*)_Ter([1-9]\d*)delins((?:[A-Z][a-z]{2})+)")
STOP_FRAMESHIFT = re.compile(r"([A-Z][a-z]{2})([1-9]\d*)([A-Z][a-z]{2})fsTer([1-9]\d*|\?)")
STOP_DELETED_EXTENSION = re.compile(r"Ter([1-9]\d*)delextTer([1-9]\d*|\?)")
RESIDUE_AND_STOP_INSERTED = re.compile(r"([A-Z][a-z]{2})([1-9]\d*)_([A-Z][a-z]{2})([1-9]\d*)ins([A-Z][a-z]{2})Ter")
# An HGVSc body of any kind but a repeat count: the numbering, the one or two cDNA coordinates,
# then a substitution (`-47A>C`; the reference base is `N` over a masked reference window), a
# `delins` or `ins` with its bases, or a plain inversion, deletion or duplication.
HGVSC_CLIPPED_BODY = re.compile(r"(?P<numbering>[cn])\.(?P<start>\*?-?\d+(?:[+-]\d+)?)(?:_(?P<end>\*?-?\d+(?:[+-]\d+)?))?"
                                r"(?:[ACGTN]>(?P<alternate_base>[ACGT])|(?P<kind>delins|ins)(?P<bases>[ACGT]+)|(?P<plain>inv|del|dup))")
# The fields the HGVSIntronOffset plugin writes beside an HGVSc, and every field an HGVSc the
# reference lacks brings with it.
HGVS_INTRON_OFFSET_FIELDS = ("HGVS_IntronStartOffset", "HGVS_IntronEndOffset")
# The exon coordinate and the intron offset `hgvs_transcript` reads off a cDNA coordinate
# (TranscriptVariationAllele.pm:1530-1531): the first integer and the one after it.
COORDINATE_PARTS_AS_READ = re.compile(r"(-?[0-9]+)\+?(-?[0-9]+)?")
# A coordinate one base before an exonic one, as `_get_cDNA_position` writes the position
# before a forward-strand transcript's first base.
ONE_BEFORE_COORDINATE = re.compile(r"(\*?-?\d+)-1")
SHIFTED_INDEL_GATE_FIELDS = ("HGVSc", "HGVS_OFFSET", *HGVS_INTRON_OFFSET_FIELDS)
# The terms of a row annotated outside the transcript, on which Ensembl VEP prints no HGVSc
# (OutputFactory.pm:1721).
OUTSIDE_TRANSCRIPT_TERMS = frozenset(("upstream_gene_variant", "downstream_gene_variant"))
# A non-zero HGVS_OFFSET, the only value Ensembl VEP prints for the field.
NON_ZERO_INTEGER = re.compile(r"-?[1-9]\d*")
# An HGVSc body describing an insertion or deletion: the numbering, the one or two cDNA
# coordinates and the kind (`c.520_541del`, `n.14dup`, `c.-258_-257insCGAGGGA`).
HGVSC_INDEL_BODY = re.compile(r"([cn])\.(\*?-?\d+(?:[+-]\d+)?)(?:_(\*?-?\d+(?:[+-]\d+)?))?(del|dup|ins[ACGT]+)")
# A cDNA coordinate of an HGVSc body (`-71`, `250+601`, `*12-3`); with each replaced by `#` two
# bodies that differ in their coordinates alone compare equal.
HGVSC_COORDINATE = re.compile(r"\*?-?\d+(?:[+-]\d+)?")
# The signed integers of an HGVSc body: `250+601` gives 250 and 601, `*12-3` gives 12 and -3.
HGVSC_INTEGER = re.compile(r"-?\d+")
# The class of an `HGVSp` Ensembl VEP writes and vep-rs leaves empty on an insertion the input
# writes between an intron's last base and the first coding base of the next exon
# (`intron_end_insertion`), and the four forms of HGVSc body such a row carries, `N` the exon's
# first coding base: a copy of the intron's last base (`c.N-1dup`), a duplication of two or more
# intron positions ending on it (`c.N-K_N-1dup`, or `c.M+J_N-1dup` from the intron's first half),
# an insertion after it (`c.N-1_Nins` and the inserted bases) of one base or of bases that do not
# end in a pair that ends an intron, and an insertion there of two or more bases that end in such
# a pair, so that every `c.N-1_Nins` insertion is of a form.
HGVSP_INTRON_END_INSERTION = "hgvsp_intron_end_insertion"
INTRON_END_COPY_BODY = re.compile(r"c\.([1-9]\d*)-1dup")
INTRON_END_DUPLICATION_BODY = re.compile(r"c\.([1-9]\d*)([+-])([1-9]\d*)_([1-9]\d*)-1dup")
INTRON_END_INSERTION_BODY = re.compile(r"c\.([1-9]\d*)-1_([1-9]\d*)ins([A-Z]+)")
# The two pairs of bases that end an intron, of the major kind and of the minor kind, as
# character codes, on the transcript's strand as an HGVSc spells inserted bases.
INTRON_END_PAIRS = frozenset(bytes(pair).decode() for pair in ((65, 71), (65, 67)))
# The longest intron at whose end vep-rs keeps the notation: Ensembl VEP reads an intron of this
# many bases or fewer with the exons beside it.
SHORT_INTRON_BASES = 13
# The class of an `HGVSp` Ensembl VEP writes and vep-rs leaves empty on an insertion the 3' rule
# leaves between a coding exon's last base and the first base of the intron after it
# (`intron_start_insertion`), and the form of HGVSc body such a row carries, `N` the exon's last
# coding base in both coordinates (`c.N_N+1ins` and the inserted bases).
HGVSP_INTRON_START_INSERTION = "hgvsp_intron_start_insertion"
INTRON_START_INSERTION_BODY = re.compile(r"c\.([1-9]\d*)_([1-9]\d*)\+1ins([A-Z]+)")
# The base each base pairs with on the other strand, as character codes: a record that holds an
# insertion with a reference base it shares spells its allele on the forward strand, and the
# HGVSc of either class spells the inserted bases on the transcript's (`held_insertion_gaps`).
BASE_COMPLEMENT = {65: 84, 67: 71, 71: 67, 84: 65}
# The flags that drop a record's rows by what they overlap (OutputFactory.pm:482-500 and
# :540-558), so a record of a run under one may have no row on either side.
ROW_FILTERS = frozenset({"--coding_only", "--no_intergenic"})
# The flags under which Ensembl VEP's variation source drops whole records by the frequency of
# their known alleles, with no warning (BaseCacheVariation.pm:308-375), so a record of such a
# run may have no reference row and rows on the vep-rs side.
FREQUENCY_FILTER_FLAGS = frozenset({"--check_frequency", "--filter_common"})
# The class of a record Ensembl VEP 116.2 dropped by frequency and vep-rs keeps because one of
# its alleles has a population frequency of exactly 0: 116.2 read the zero as no frequency
# (`|| 'NA'`, BaseCacheVariation.pm:419), which fails the threshold test, where vep-rs judges
# the zero as Ensembl VEP 116.3 does (`// 'NA'`). Its vep-rs rows carry a FREQS item of value 0.
FREQUENCY_FILTER_ZERO_VALUE = "frequency_filter_zero_value"
# The class of a plugin field on a row of a multi-allelic record `--minimal` split, where Ensembl
# VEP's output pass hands the plugin the rejoined record's span and alleles (OutputFactory.pm:2291-2293,
# :2314, :2226) and vep-rs the allele's own; `rejoined_record_plugin_field` names the fields of the
# covered plugins whose names are fixed (dbNSFP's and TranscriptAnnotator's are their tables' column
# names, so a dataset reaching them extends the predicate).
MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD = "minimal_plugin_sees_the_rejoined_record"
REJOINED_RECORD_PLUGIN_FIELDS = frozenset((
    "Conservation", "NearestExonJB", "TSSDistance", "ClinPred", "LOEUF", "PrimateAI", "REVEL", "ada_score",
    "rf_score", "pHaplo", "pTriplo", "Existing_InFrame_oORFs", "Existing_OutOfFrame_oORFs", "Existing_uORFs",
))
REJOINED_RECORD_PLUGIN_PREFIXES = (
    "gnomAD_", "GWAS_", "AVI_", "OpenTargets_", "am_", "BayesDel_", "CADD_", "MaxEntScan_", "SpliceAI_pred",
    "SpliceVault_", "5UTR_",
)
# `maxentscan_overlapping_windows`: on a row whose one-base record lies in the acceptor window of
# one intron of the row's transcript and the donor window of the next, Ensembl VEP's MaxEntScan
# scores the window of whichever intron its interval tree returns first (MaxEntScan.pm:270) and
# vep-rs the upstream intron's acceptor, so the numbers differ and the window texts of a `verbose`
# run are a donor's nine bases against an acceptor's twenty-three.
MAXENTSCAN_OVERLAPPING_WINDOWS = "maxentscan_overlapping_windows"
MAXENTSCAN_FIELD_PREFIX = "MaxEntScan_"
MAXENTSCAN_WINDOW_TEXT_FIELDS = ("MaxEntScan_ref_seq", "MaxEntScan_alt_seq")
# The class of a `DownstreamProtein` on a row of a deletion, with or without inserted bases, whose
# first base is a codon's third: Ensembl VEP's Downstream plugin rounds that base's coordinate down
# to a codon boundary, the end of the codon the deletion changes, and translates from the codon
# after it (Downstream.pm:112-115), where vep-rs translates from the changed codon.
DOWNSTREAM_START_PAST_OWN_CODON = "downstream_start_past_own_codon"
# The class of such a row whose deletion inserts no base and runs to the last coding base: the
# changed coding sequence ends one base before that boundary, so the plugin takes nothing of it
# (Downstream.pm:115) and translates the 3' UTR from its first base (Downstream.pm:116-119), one
# base out of the changed sequence's frame, where vep-rs translates the two bases the changed
# codon keeps with the UTR, in frame.
DOWNSTREAM_UTR_FRAME_PAST_CODING_END = "downstream_utr_frame_past_coding_end"
# The class of a `DownstreamProtein` or `ProteinLengthChange` on a row whose HGVS notation is
# shifted 3': under `--hgvs` Ensembl VEP's `hgvs_protein` reads the row's CDS pair again at the
# shifted span and leaves it there (TranscriptVariationAllele.pm:844-858, :2406-2410) with the
# allele's bases rotated for that span (:1330-1345, :1778), and the Downstream plugin, run after
# it (OutputFactory.pm:398), computes on both (Downstream.pm:108-110), so its two values change
# with a flag that asks for a notation; vep-rs computes on the record's own span and the allele
# as written, with the flag and without it.
PLUGIN_HGVS_SHIFTED_CDS_COORDINATES = "plugin_hgvs_shifted_cds_coordinates"
# The key under which a `field_divergences` entry names the class of a second cause that changes
# the same value: the format comparator applies such an entry only where both classes are named.
SECOND_EXPECTED_DIVERGENCE = "second_expected_divergence"
# The key under which a `divergences` entry states vep-rs's consequence sets in a run that has a
# reference sequence, `vep_rs_consequence_sets` being those of a run without one; a set equal to
# `vep_consequence_set` states that the two engines agree there. Written by hand, carried by classify.
WITH_REFERENCE_SETS = "vep_rs_consequence_sets_with_reference_sequence"
# The fields the Downstream plugin writes, in the order `downstream_plugin_values` and
# `downstream_own_values` give them.
DOWNSTREAM_FIELDS = ("DownstreamProtein", "ProteinLengthChange")
# The keys of a cached mapper pair, in the order `CodingModel.pairs` holds their values: the
# pair's first and last cDNA base, its first and last genomic base and its orientation.
MAPPER_PAIR_KEYS = ("from_start", "from_end", "to_start", "to_end", "ori")
# A `CDS_position` column naming one coding base or a span of them, neither end a gap (`?`).
CDS_POSITION_SPAN = re.compile(r"(\d+)(?:-(\d+))?")
# The genetic codes Bio::Tools::CodonTable translates with, by table id, as `CODES` of
# crates/vep-core/src/codon.rs lists them: each codon's amino acid, the codons in the order TTT TTC
# TTA TTG TCT ... GGG. Ids 7, 8 and 17 to 20 hold no table there and ids from 25 on are unknown to
# it, so a transcript under one of them has no code here.
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
# Each table's amino acid by codon.
GENETIC_CODES = {
    table: {first + second + third: amino_acids[16 * i + 4 * j + k]
            for i, first in enumerate("TCAG") for j, second in enumerate("TCAG") for k, third in enumerate("TCAG")}
    for table, amino_acids in CODON_TABLE_AMINO_ACIDS.items()
}
# The bases an IUPAC nucleotide letter stands for, `U` read as `T`.
IUPAC_BASES = {
    "T": "T", "U": "T", "C": "C", "A": "A", "G": "G", "Y": "TC", "W": "TA", "K": "TG", "M": "CA", "S": "CG",
    "R": "AG", "H": "TCA", "B": "TCG", "D": "TAG", "V": "CAG", "N": "TCAG", "X": "TCAG",
}
# Each base's complement, for an allele read on a reverse-strand transcript.
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def rejoined_record_plugin_field(name: str) -> bool:
    """Whether `name` is a field of a plugin the class `minimal_plugin_sees_the_rejoined_record` covers."""
    return name in REJOINED_RECORD_PLUGIN_FIELDS or name.startswith(REJOINED_RECORD_PLUGIN_PREFIXES)


def indel_edited_models(cache_dir: Path) -> dict[str, frozenset[int]]:
    """The cache's transcripts carrying an `_rna_edit` attribute that inserts or deletes bases
    against the genome, by stable id, each with the distances its edits put between the edited
    and the unedited cDNA numbering: the running sums of the edits' length changes in start
    order, as absolute values, 0 left out. The value is read as Ensembl VEP's
    `get_misalignment_offset` reads it (TranscriptVariationAllele.pm:1605-1618 at release/116):
    `<start> <end> <bases>` split on single spaces with trailing empty fields dropped, a
    two-field value a deletion of `end - start + 1` bases, a three-field value whose bases span
    exactly that many a substitution and any other an insertion of its bases; a value without
    two leading integers, which no cache writes, is skipped. Empty when the cache is absent."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc

    edited: dict[str, frozenset[int]] = {}
    for _chrom, tx in bmc.iter_cache_transcripts(cache_dir):
        edits: list[tuple[int, int]] = []
        for attribute in tx.get("attributes") or []:
            if not isinstance(attribute, dict) or not str(attribute.get("code", "")).startswith("_rna_edit"):
                continue
            fields = str(attribute.get("value", "")).split(" ")
            while fields and fields[-1] == "":
                fields.pop()
            if len(fields) < 2 or not all(re.fullmatch(r"[+-]?[0-9]+", f) for f in fields[:2]):
                continue
            start, end = int(fields[0]), int(fields[1])
            if len(fields) == 3:
                if end - start + 1 == len(fields[2]):
                    continue
                edits.append((start, len(fields[2])))
            else:
                edits.append((start, -1 - (end - start)))
        if not edits:
            continue
        offsets: set[int] = set()
        running = 0
        for _start, length_change in sorted(edits):
            running += length_change
            offsets.add(abs(running))
        offsets.discard(0)
        edited[str(tx.get("stable_id", ""))] = frozenset(offsets)
    edited.pop("", None)
    return edited


def cdna_genomic_position(exons: list[tuple[int, int]], forward: bool, cdna: int) -> int | None:
    """The genomic position of the cDNA position `cdna`, counted along `exons` (in genomic order)
    in the transcript's direction; `None` for a position outside them."""
    remaining = cdna
    for start, end in exons if forward else reversed(exons):
        length = end - start + 1
        if 1 <= remaining <= length:
            return start + remaining - 1 if forward else end - remaining + 1
        remaining -= length
    return None


def coding_transcript_exons(tx: dict) -> tuple[bool, list[tuple[int, int]], int, int] | None:
    """Whether a cache transcript lies on the forward strand, its exons in genomic order, and the
    lower and the higher genomic bound of its coding region, stop codon included: the genomic
    positions of its `cdna_coding_start` and `cdna_coding_end` (`cdna_genomic_position`). `None`
    for a transcript without a coding region and for one with a coordinate that is no integer."""
    try:
        forward = int(tx.get("strand", 1)) >= 0
        exons = sorted((int(exon["start"]), int(exon["end"])) for exon in tx.get("exons") or [])
        coding = [int(tx[key]) for key in ("cdna_coding_start", "cdna_coding_end")]
    except (KeyError, TypeError, ValueError):
        return None
    bounds = [cdna_genomic_position(exons, forward, cdna) for cdna in coding]
    if None in bounds:
        return None
    return forward, exons, min(bounds), max(bounds)


def transcript_intron_end_gaps(tx: dict) -> set[int]:
    """The gaps of one cache transcript in which an intron of more than `SHORT_INTRON_BASES` ends
    before an exon base inside the coding region, each the lower of the two genomic positions
    either side of the gap: on the forward strand the intron's last base, the exon's first base
    following it; on the reverse strand the exon's first base in the transcript's direction, its
    highest genomic base, the intron's last base following it. Empty where
    `coding_transcript_exons` reads no coding region."""
    coding = coding_transcript_exons(tx)
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


def transcript_intron_start_gaps(tx: dict) -> set[int]:
    """The gaps of one cache transcript in which an exon ends inside the coding region before an
    intron of any length, each the lower of the two genomic positions either side of the gap: on
    the forward strand the exon's last base, the intron's first base following it; on the reverse
    strand the intron's first base in the transcript's direction, its highest genomic base, the
    exon's last base in that direction following it. Empty where `coding_transcript_exons` reads
    no coding region."""
    coding = coding_transcript_exons(tx)
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


def cache_transcript_gaps(cache_dir: Path, gaps_of: Callable[[dict], set[int]]) -> dict[tuple[str, str], frozenset[int]]:
    """The gaps `gaps_of` reads from each transcript of the cache, by chromosome and stable id; a
    transcript without a gap has no entry. Empty when the cache is absent."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc

    gaps: dict[tuple[str, str], set[int]] = defaultdict(set)
    for chrom, tx in bmc.iter_cache_transcripts(cache_dir):
        found = gaps_of(tx)
        if found:
            gaps[(chrom, str(tx["stable_id"]))].update(found)
    return {key: frozenset(found) for key, found in gaps.items()}


def intron_end_gaps(cache_dir: Path) -> dict[tuple[str, str], frozenset[int]]:
    """The gaps of the cache's transcripts in which an intron ends before a coding exon base
    (`transcript_intron_end_gaps`), by chromosome and stable id (`cache_transcript_gaps`)."""
    return cache_transcript_gaps(cache_dir, transcript_intron_end_gaps)


def intron_start_gaps(cache_dir: Path) -> dict[tuple[str, str], frozenset[int]]:
    """The gaps of the cache's transcripts in which a coding exon ends before an intron
    (`transcript_intron_start_gaps`), by chromosome and stable id (`cache_transcript_gaps`)."""
    return cache_transcript_gaps(cache_dir, transcript_intron_start_gaps)


def transcript_two_window_positions(tx: dict) -> set[int]:
    """The genomic positions of one cache transcript that lie in the acceptor window of one intron
    and the donor window of the next intron in the transcript's order, the windows as
    MaxEntScan.pm:282-291 places them: the donor's from three bases before an intron to its sixth
    base, the acceptor's from an intron's twentieth base from its end to three bases after it.
    Empty for a transcript with a coordinate that is no integer."""
    try:
        forward = int(tx.get("strand", 1)) >= 0
        exons = sorted((int(exon["start"]), int(exon["end"])) for exon in tx.get("exons") or [])
    except (KeyError, TypeError, ValueError):
        return set()
    introns = [(exons[index][1] + 1, exons[index + 1][0] - 1) for index in range(len(exons) - 1)]
    positions: set[int] = set()
    for index in range(len(introns) - 1):
        lower_end, upper_start = introns[index][1], introns[index + 1][0]
        if forward:
            acceptor, donor = (lower_end - 19, lower_end + 3), (upper_start - 3, upper_start + 5)
        else:
            acceptor, donor = (upper_start - 3, upper_start + 19), (lower_end - 5, lower_end + 3)
        positions.update(range(max(acceptor[0], donor[0]), min(acceptor[1], donor[1]) + 1))
    return positions


def two_window_positions(cache_dir: Path) -> dict[tuple[str, str], frozenset[int]]:
    """The positions of the cache's transcripts that lie in two canonical splice windows
    (`transcript_two_window_positions`), by chromosome and stable id (`cache_transcript_gaps`)."""
    return cache_transcript_gaps(cache_dir, transcript_two_window_positions)


def finite_number(value: str) -> bool:
    """Whether a field's text is a finite number."""
    try:
        number = float(value)
    except ValueError:
        return False
    return float("-inf") < number < float("inf")


def maxentscan_overlapping_windows(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                   two_windows: dict[tuple[str, str], frozenset[int]] | None) -> bool:
    """Whether a differing `MaxEntScan_` pair stands on a row whose Location names one position
    that lies in two canonical windows of the row's transcript in `two_windows`
    (`two_window_positions`), under the Location's chromosome as written or without a leading
    `chr` and the row's Feature as written or without a trailing version, and has the class's
    shape: two finite numbers, or for a window text of a `verbose` run Ensembl VEP's nine bases
    against vep-rs's twenty-three. Without the positions nothing is named."""
    location, _allele, feature, _feature_type, _consequences = key
    if not name.startswith(MAXENTSCAN_FIELD_PREFIX) or not two_windows:
        return False
    chrom, _, position = location.rpartition(":")
    if not position.isdigit():
        return False
    names = (feature, re.sub(r"\.\d+$", "", feature))
    if not any(int(position) in two_windows.get((contig, stable_id), ())
               for contig in (chrom, chrom.removeprefix("chr")) for stable_id in names):
        return False
    if name in MAXENTSCAN_WINDOW_TEXT_FIELDS:
        return re.fullmatch(r"[ACGT]{9}", vep) is not None and re.fullmatch(r"[ACGT]{23}", vep_rs) is not None
    return finite_number(vep) and finite_number(vep_rs)


def peptide_description_shape(description: str) -> str:
    """`dup`, `ins` or `other` for the text after `p.`; a `delins` is `other`."""
    if description.endswith("dup"):
        return "dup"
    if re.search(r"(?<!del)ins", description):
        return "ins"
    return "other"


def organelle_peptide_path_divergence(consequences: str, vep: str, vep_rs: str) -> bool:
    """Whether a differing HGVSp pair of a row with the given consequence set comes out of one of
    the four paths of TranscriptVariationAllele.pm that translate a sequence bare (the alternate
    sequence in `_get_fs_peptides`, `_stop_loss_extra_AA` and `_get_del_peptides`, the reference
    sequence in `_check_for_peptide_duplication`): `_get_fs_peptides` and `_stop_loss_extra_AA`
    write the whole description of a frameshift (its residues, a `Ter` at an immediate stop, the
    stop distance and the lost-start shapes); `_stop_loss_extra_AA` the count of an extension
    after a lost stop; `_check_for_peptide_duplication` the choice between a duplication and an
    insertion. `_get_del_peptides` writes a deletion whose reference window is empty: its rows
    re-type to an insertion or a synonymous description and carry no `ext`, no corpus row reaches
    it, and the rule does not name it, so such a row stays unexplained. A substitution, a `delins`
    and a plain deletion read their residues through `peptide`, which already follows the
    transcript's table."""
    terms = set(consequences.split(","))
    if "frameshift_variant" in terms:
        return True
    if "stop_lost" in terms:
        return "ext" in vep and "ext" in vep_rs
    shapes = {peptide_description_shape(vep), peptide_description_shape(vep_rs)}
    return bool(terms & {"inframe_insertion", "protein_altering_variant"}) and shapes == {"dup", "ins"}


def organelle_codon_tables(corpus: Path) -> dict[str, int]:
    """Each transcript of the corpus's cache whose `codon_table` is not the standard code (1),
    by stable id, with that table: the transcripts whose peptides Ensembl VEP translates
    with two codes."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc

    tables: dict[str, int] = {}
    for _chrom, tx in bmc.iter_cache_transcripts(corpus_cache_dir(corpus, read_manifest(corpus))):
        vefc = tx.get("variation_effect_feature_cache") or {}
        # A converted cache carries the table as the dump's string (`"5"`), a native one as a number.
        table = str(vefc.get("codon_table", ""))
        if table.isdigit() and int(table) not in (0, 1):
            tables[tx["stable_id"]] = int(table)
    return tables


def downstream_sequences(corpus: Path) -> dict[str, tuple[str, str, int]]:
    """Each transcript of the corpus's cache that carries a `translateable_seq`, by stable id, with
    the sequences Ensembl VEP's Downstream plugin translates and the table it translates them
    under: the `translateable_seq` (Downstream.pm:104-106), the `three_prime_utr`, the empty
    sequence when the cache spells none (the plugin reads the transcript's own, :116), and the
    `codon_table`, 1 when the cache names none or 0 (:130-133). A transcript whose table is no
    number is left out."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc

    sequences: dict[str, tuple[str, str, int]] = {}
    for _chrom, tx in bmc.iter_cache_transcripts(corpus_cache_dir(corpus, read_manifest(corpus))):
        vefc = tx.get("variation_effect_feature_cache") or {}
        coding, utr = vefc.get("translateable_seq"), vefc.get("three_prime_utr")
        # A converted cache carries the table as the dump's string (`"1"`), a native one as a number.
        table = str(vefc.get("codon_table") or 1)
        if isinstance(coding, str) and coding and table.isdecimal():
            sequences[tx["stable_id"]] = (coding, utr if isinstance(utr, str) else "", int(table) or 1)
    return sequences


@dataclass(frozen=True)
class CodingModel:
    """What places a genomic base on a cached transcript's coding sequence: the mapper pairs,
    each the values of `MAPPER_PAIR_KEYS`, the first and last cDNA base of the coding sequence,
    and the length of the cached peptide, `None` where the cache holds none."""
    pairs: tuple[tuple[int, ...], ...]
    coding_start: int
    coding_end: int
    peptide_length: int | None

    def cdna(self, genomic: int) -> int | None:
        """The cDNA coordinate of a genomic base, or `None` for a base no pair holds."""
        for from_start, _from_end, to_start, to_end, ori in self.pairs:
            if to_start <= genomic <= to_end:
                return from_start + (to_end - genomic if ori == -1 else genomic - to_start)
        return None

    def cds_pair(self, start: int, end: int, reverse: bool, padding: int) -> tuple[int, int] | None:
        """The CDS start and end Ensembl VEP gives the genomic span `start` to `end` on the
        transcript, `reverse` for one on the reverse strand, counted from 1 on a
        `translateable_seq` that opens with `padding` bases of start-exon phase
        (BaseTranscriptVariation.pm:252-264), or `None` where either is no coding base. A span
        of bases maps by its two ends. An insertion, whose start lies one past its end, takes
        the bases between its two flanks' coordinates (Mapper.pm:504-558): none between
        neighbours, which leaves the start one past the end, and the bases the model carries
        where two pairs touch on the genome; where no pair holds one flank it lies against the
        other. Each flank's coordinate is read alone, so a flank outside the coding sequence
        leaves no pair even where the bases between reach into it (TranscriptMapper.pm:441-443)."""
        if start == end + 1:
            before, after = (self.cdna(base) for base in ((start, end) if reverse else (end, start)))
            if before is None and after is None:
                return None
            if (before is not None and before < self.coding_start) or (after is not None and after > self.coding_end):
                return None
            first = after if before is None else before + 1
            last = before if after is None else after - 1
        else:
            first, last = (self.cdna(base) for base in ((end, start) if reverse else (start, end)))
            if first is None or last is None:
                return None
        if not (self.coding_start <= first <= self.coding_end and self.coding_start <= last <= self.coding_end):
            return None
        to_cds = padding + 1 - self.coding_start
        return first + to_cds, last + to_cds


def downstream_coding_models(corpus: Path) -> dict[tuple[str, str], CodingModel]:
    """Each transcript of the corpus's cache that carries a `cdna_coding_start`, a
    `cdna_coding_end` and mapper pairs, by chromosome and stable id, as a `CodingModel`; a
    converted cache carries the pairs' numbers as the dump's strings. The chromosome is part of
    the key because a cache files a pseudoautosomal transcript under X and under Y, with its
    own pairs on each. A transcript lacking one of the three, or holding a value that is no
    number, is left out."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc

    models: dict[tuple[str, str], CodingModel] = {}
    for chrom, tx in bmc.iter_cache_transcripts(corpus_cache_dir(corpus, read_manifest(corpus))):
        vefc = tx.get("variation_effect_feature_cache") or {}
        peptide = vefc.get("peptide")
        try:
            pairs = tuple(tuple(int(pair[key]) for key in MAPPER_PAIR_KEYS) for pair in vefc["mapper"]["pairs"])
            models[chrom, tx["stable_id"]] = CodingModel(pairs, int(tx["cdna_coding_start"]), int(tx["cdna_coding_end"]),
                                                         len(peptide) if isinstance(peptide, str) else None)
        except (KeyError, TypeError, ValueError):
            continue
    return models


def record_span(vcf_line: str) -> tuple[int, int]:
    """The genomic start and end Ensembl VEP's VCF parser gives a record (Parser/VCF.pm:254,
    :295-336): POS and the last REF base, the start one base later where an ALT differs from REF
    in length and REF and the ALTs open with one base, an ALT naming `*` left out of several, so
    an insertion's start lies one past its end."""
    f = vcf_line.split("\t")
    start, ref, alts = int(f[1]), f[3], f[4].split(",")
    end = start + len(ref) - 1
    if any(len(alt) != len(ref) for alt in alts):
        first_bases = {allele[:1] for allele in (ref, *alts) if len(alts) == 1 or "*" not in allele}
        if len(first_bases) == 1:
            start += 1
    return start, end


def row_record_span(location: str, vcf_lines: Iterable[str]) -> tuple[int, int] | None:
    """The one span among the records `vcf_lines` (`record_span`) that a row's Location prints,
    an insertion's two flanks lower first, or `None` where none or more than one does."""
    printed = location.rsplit(":", 1)[-1]
    spans: set[tuple[int, int]] = set()
    for vcf_line in vcf_lines:
        span = record_span(vcf_line)
        low, high = sorted(span)
        if printed == (str(low) if low == high else f"{low}-{high}"):
            spans.add(span)
    return spans.pop() if len(spans) == 1 else None


def ambiguous_codon_amino_acid(codon: str, code: dict[str, str]) -> str:
    """`Bio::Tools::CodonTable::_translate_ambiguous_codon` under the genetic code `code`: the
    amino acid of an upper-case codon that is no unambiguous codon, from the distinct amino acids
    its readings code: one, that one; aspartate and asparagine `B`; glutamate and glutamine `Z`;
    anything else `X`, as is a codon with a character that is no nucleotide letter."""
    readings = [""]
    for letter in codon:
        readings = [reading + base for reading in readings for base in IUPAC_BASES.get(letter, "")]
    amino_acids = {code[reading] for reading in readings}
    if len(amino_acids) == 1:
        return amino_acids.pop()
    if amino_acids == {"D", "N"}:
        return "B"
    if amino_acids == {"E", "Q"}:
        return "Z"
    return "X"


def translate_codons(sequence: str, table: int) -> str:
    """`Bio::Tools::CodonTable::translate` under the genetic code of id `table`, one of
    `GENETIC_CODES`, as `Bio::Seq::translate` calls it for a sequence with no frame and no
    partial-codon reading and as `translate_sequence` of crates/vep-core/src/codon.rs restates
    it: the sequence's whole codons in order, in either case, a trailing one or two bases
    dropped; a codon of unambiguous bases by the table, `---` as `-`, any other by
    `ambiguous_codon_amino_acid`. A stop stays `*`."""
    code = GENETIC_CODES[table]
    sequence = sequence.upper()
    peptide: list[str] = []
    for at in range(0, len(sequence) - 2, 3):
        codon = sequence[at:at + 3]
        if codon in code:
            peptide.append(code[codon])
        elif codon == "---":
            peptide.append("-")
        else:
            peptide.append(ambiguous_codon_amino_acid(codon, code))
    return "".join(peptide)


def downstream_plugin_values(coding: str, utr: str, table: int, peptide_length: int,
                             start: int, end: int, bases: str) -> tuple[str, str, bool] | None:
    """Ensembl VEP's `DownstreamProtein` and `ProteinLengthChange` (Downstream.pm:104-149) under
    the codon table `table`, one of `GENETIC_CODES`, for a change that puts `bases`, the allele
    on the transcript's strand, in place of the coding bases `start` to `end` (counted from 1 on
    `coding`, the transcript's `translateable_seq`; an insertion's start lies one past its end),
    and whether the protein's translation starts past the changed coding sequence; `None` where
    `start` lies more than one base past `coding`, which ends the plugin's call (:110). The
    protein is the translation of the changed sequence from a codon boundary, the 3' UTR
    appended, up to its first stop (:115-140). The boundary is the lower of `start` and `end`
    rounded down to a multiple of three (:112-113), the end of the codon holding that base when
    it is the codon's third, and a boundary past the changed sequence's end takes none of it
    (:115), which leaves the UTR translated from its first base. The length change is the
    length of the whole changed sequence's translation with the UTR, up to its first stop,
    less `peptide_length`, the cached peptide's (:124-148)."""
    if start - 1 > len(coding):
        return None
    changed = coding[:start - 1] + bases + coding[max(end, start - 1):]
    low = min(start, end)
    boundary = low - low % 3
    protein = translate_codons(changed[boundary:] + utr, table).partition("*")[0]
    whole = translate_codons(changed + utr, table).partition("*")[0]
    return protein, str(len(whole) - peptide_length), boundary > len(changed)


def downstream_own_values(coding: str, utr: str, table: int, peptide_length: int,
                          start: int, end: int, bases: str) -> tuple[str, str]:
    """vep-rs's `DownstreamProtein` and `ProteinLengthChange` for the change
    `downstream_plugin_values` describes (`translate_change` of
    crates/vep-builtins/src/plugins/downstream.rs): the span is clipped to `coding`, the protein
    is read from the first base of the codon holding the first changed base, and the length
    change is computed as Ensembl VEP's."""
    first = min(start - 1, len(coding))
    changed = coding[:first] + bases + coding[min(max(end, start - 1), len(coding)):]
    own_codon = first - first % 3
    protein = translate_codons(changed[own_codon:] + utr, table).partition("*")[0]
    whole = translate_codons(changed + utr, table).partition("*")[0]
    return protein, str(len(whole) - peptide_length)


def hgvs_shifted_allele(allele: str, offset: int, reverse: bool) -> str:
    """The bases `shift_feature_seqs` leaves on a row whose notation moved `offset` bases
    (TranscriptVariationAllele.pm:1330-1345): the allele as written, its first base moved to
    its end once per base of the shift on a forward-strand transcript and `len(allele)` less
    that many times on a reverse-strand one (`reverse`), not at all where that count is not
    positive."""
    turns = len(allele) - abs(offset) if reverse else abs(offset)
    if allele and turns > 0:
        turns %= len(allele)
        allele = allele[turns:] + allele[:turns]
    return allele


def downstream_field_class(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str, strand: str,
                           perl_row: dict[str, str] | None, rust_row: dict[str, str] | None,
                           coding_sequences: dict[str, tuple[str, str, int]] | None,
                           coding_models: dict[tuple[str, str], CodingModel] | None,
                           span: tuple[int, int] | None) -> tuple[str, str | None] | None:
    """The class of a differing `DownstreamProtein` or `ProteinLengthChange` pair (`name`) with
    the class of a second cause that changes the value (`None` for a row one cause changes), or
    `None`. The row carries `frameshift_variant` and no term naming a splice site, the rows the
    plugin annotates (Downstream.pm:95-99); both engines wrote a value and print one
    `CDS_position` with no gap; its transcript is one of `coding_sequences`
    (`downstream_sequences`), and the values are read under the transcript's codon table, so a
    transcript under a table `GENETIC_CODES` lacks gets no class.

    The row has two readings. Its own is the `CDS_position` column with the row's allele on the
    transcript's strand, none for `-`: one base or `start-end` in order, and for an insertion
    (`span`, the record's genomic start and end of `row_record_span`, has its start one past
    its end) the two flanking bases the column prints lower first, read as a start one past its
    end unless the transcript's model (`coding_models`, `downstream_coding_models`) carries
    exactly those bases at the insertion's point. The shifted one is what the plugin reads
    under `--hgvs` on a row whose `HGVS_OFFSET` in `perl_row` is a non-zero integer: the CDS
    pair of `span` moved by that offset along the genome (`CodingModel.cds_pair`) with the
    allele of `hgvs_shifted_allele`, where both ends of the moved span are coding bases
    (TranscriptVariationAllele.pm:1730-1738 returns before the move otherwise); on any other
    row, and without `span` or a model, it is the own reading.

    vep-rs's value must be `downstream_own_values` on the own reading and Ensembl VEP's
    `downstream_plugin_values` on the shifted one. Where Ensembl VEP's rule gives another value
    on the own reading, the flag changes its value and the class is
    `plugin_hgvs_shifted_cds_coordinates`, for either field; a `ProteinLengthChange` needs the
    cached peptide's length and has no other class. Otherwise the two rules differ on the own
    reading itself, a span on the coding sequence whose first base is a codon's third:
    `downstream_utr_frame_past_coding_end` where Ensembl VEP's translation starts past the
    changed coding sequence, else `downstream_start_past_own_codon`. A `DownstreamProtein` the
    flag changes whose vep-rs value is not Ensembl VEP's on the own reading either is changed
    by both causes: its second class is the one its own reading has."""
    location, allele, feature, _feature_type, consequences = key
    terms = consequences.split(",")
    if "frameshift_variant" not in terms or any("splice" in term for term in terms):
        return None
    coding, utr, codon_table = (coding_sequences or {}).get(feature, ("", "", 1))
    perl_row = perl_row or {}
    column = perl_row.get("CDS_position", "")
    printed = CDS_POSITION_SPAN.fullmatch(column)
    if not coding or codon_table not in GENETIC_CODES or printed is None or "-" in (vep, vep_rs):
        return None
    model = (coding_models or {}).get((location.rsplit(":", 1)[0], feature))
    peptide_length = model.peptide_length if model is not None else None
    index = DOWNSTREAM_FIELDS.index(name)
    if (rust_row or {}).get("CDS_position", column) != column or (index and peptide_length is None):
        return None
    reverse = strand == "-1"
    padding = len(coding) - len(coding.lstrip("Nn"))
    written = "" if allele == "-" else allele

    def on_strand(bases: str) -> str:
        return bases.translate(COMPLEMENT)[::-1] if reverse else bases

    start, end = int(printed.group(1)), int(printed.group(2) or printed.group(1))
    insertion = span is not None and span[0] == span[1] + 1
    if insertion and end == start + 1 and (model is None or model.cds_pair(*span, reverse, padding) != (start, end)):
        start, end = end, start
    if not 1 <= start <= end + 1 <= len(coding) + 1:
        return None
    own = shifted = (start, end, on_strand(written))
    offset = perl_row.get("HGVS_OFFSET", "")
    if model is not None and span is not None and NON_ZERO_INTEGER.fullmatch(offset):
        moved = model.cds_pair(span[0] + int(offset), span[1] + int(offset), reverse, padding)
        if moved is not None:
            shifted = (*moved, on_strand(hgvs_shifted_allele(written, int(offset), reverse)))
    sequences = (coding, utr, codon_table, peptide_length or 0)
    reference_own = downstream_plugin_values(*sequences, *own)
    reference = downstream_plugin_values(*sequences, *shifted)
    if (reference is None or reference_own is None or vep != reference[index]
            or vep_rs != downstream_own_values(*sequences, *own)[index]):
        return None
    own_class = DOWNSTREAM_UTR_FRAME_PAST_CODING_END if reference_own[2] else DOWNSTREAM_START_PAST_OWN_CODON
    if reference[index] != reference_own[index]:
        both_causes = not index and vep_rs != reference_own[index]
        return PLUGIN_HGVS_SHIFTED_CDS_COORDINATES, (own_class if both_causes else None)
    if index or vep == vep_rs:
        return None
    return own_class, None


def cdna_coordinate_order(coordinate: str) -> tuple[int, int, int]:
    """A key ordering HGVS cDNA coordinates along the transcript: the 5' UTR (`-N`), the coding
    sequence (`N`), the 3' UTR (`*N`), each by its position, then an intron offset (`+M` after
    the exon base, `-M` before it). `*-M`, an intron distance before the stop codon's last base
    that `_get_cDNA_position` prints without that base's coordinate
    (TranscriptVariationAllele.pm:2854-2859), orders after every coding coordinate and before
    `*1`."""
    region, number, offset = 0, coordinate, 0
    if coordinate.startswith("*"):
        region, number = 1, coordinate[1:]
    elif coordinate.startswith("-"):
        region, number = -1, coordinate[1:]
    m = re.search(r"[+-]\d+$", number)
    if m:
        offset, number = int(m.group()), number[: m.start()]
    position = int(number or 0)
    return (region, -position if region < 0 else position, offset)


def names_transcript(hgvsc: str, feature: str) -> re.Match | None:
    """The indel body of an HGVSc that names the row's transcript (its stable id, with or without
    a version), or `None` for another transcript's notation or a notation of another kind."""
    name, _, body = hgvsc.partition(":")
    if name.split(".")[0] != feature.split(".")[0]:
        return None
    return HGVSC_INDEL_BODY.fullmatch(body)


def shifted_indel_written_span_gate(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                    perl_row: dict[str, str] | None, rust_row: dict[str, str] | None) -> bool:
    """Whether a differing `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset plugin pair is Ensembl VEP's
    empty field against vep-rs's notation of a shifted insertion or deletion:
    `_var2transcript_slice_coords` reads the written span (TranscriptVariationAllele.pm:2708-2713),
    so one wholly before the transcript that the shift carries inside gets no notation (:2729)
    and one the shift carries over the last base is refused once the offset is added (:1472).
    vep-rs describes the shifted span: `HGVSc` an indel notation of the row's transcript,
    `HGVS_OFFSET` its non-zero shift and the plugin's two intron offsets integers, where the
    reference has none of them; the row's other compared fields show the same shape, and the row
    is annotated inside the transcript (no upstream or downstream term)."""
    _location, _allele, feature, _feature_type, consequences = key
    if vep not in ("", "-") or name not in SHIFTED_INDEL_GATE_FIELDS:
        return False
    if OUTSIDE_TRANSCRIPT_TERMS & set(consequences.split(",")):
        return False
    perl_row, rust_row = perl_row or {}, rust_row or {}
    if any(perl_row.get(f, "") not in ("", "-") for f in SHIFTED_INDEL_GATE_FIELDS):
        return False
    notation = vep_rs if name == "HGVSc" else rust_row.get("HGVSc", "")
    offset = vep_rs if name == "HGVS_OFFSET" else rust_row.get("HGVS_OFFSET", "")
    if name in HGVS_INTRON_OFFSET_FIELDS and not re.fullmatch(r"-?\d+", vep_rs):
        return False
    return names_transcript(notation, feature) is not None and NON_ZERO_INTEGER.fullmatch(offset) is not None


def shifted_deletion_clamped_span(key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                  perl_row: dict[str, str] | None) -> bool:
    """Whether a differing `HGVSc` pair is a shifted deletion's span clamped before the offset
    was added: Ensembl VEP clamps the written span onto the transcript's first base and then
    moves it by the shift (TranscriptVariationAllele.pm:2732-2733, 1477-1478), so its deletion
    ends where vep-rs's does and starts later, describing fewer bases than the row deletes; the
    row is a deletion's (allele `-`), both notations name its transcript, and it carries the
    shift as a non-zero `HGVS_OFFSET` on both engines. The manifest records vep-rs's exact
    notation, which the golden harness compares, so no cDNA cross-check is made here."""
    _location, allele, feature, _feature_type, _consequences = key
    if allele != "-" or not NON_ZERO_INTEGER.fullmatch((perl_row or {}).get("HGVS_OFFSET", "")):
        return False
    vep_name = vep.partition(":")[0]
    vep_rs_name = vep_rs.partition(":")[0]
    reference, candidate = names_transcript(vep, feature), names_transcript(vep_rs, feature)
    if vep_name != vep_rs_name or reference is None or candidate is None:
        return False
    if reference.group(4) != "del" or candidate.group(4) != "del" or reference.group(1) != candidate.group(1):
        return False
    reference_start, reference_end = reference.group(2), reference.group(3) or reference.group(2)
    candidate_start, candidate_end = candidate.group(2), candidate.group(3) or candidate.group(2)
    return reference_end == candidate_end and cdna_coordinate_order(candidate_start) < cdna_coordinate_order(reference_start)


def clipped_reference_against_whole_alternate(key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                              perl_row: dict[str, str] | None) -> bool:
    """Whether a differing `HGVSc` pair is Ensembl VEP's whole alternate allele read against the
    reference bases the clamp left inside the transcript: `_var2transcript_slice_coords` clamps a
    span running over the transcript's edge onto the bases inside (TranscriptVariationAllele.pm
    :2732-2733) and `hgvs_transcript` passes that span with the allele whole (:1474-1483), so a
    substitution of several bases over an edge, or a multi-allelic record's allele still carrying
    an outside base, is written as the substitution, `delins`, `ins`, `dup` or shorter `del` the
    clip leaves of every alternate base against the bases inside, as an `ins` of the outside
    bases at an intron coordinate before exon 1 (:2802) where every base inside is unchanged, and
    not at all where the whole alternate equals them (Utils/Sequence.pm:529-531); vep-rs cuts the
    alternate as the clamp cut the reference and describes the change inside, none when there is
    none. The row is a sequence allele's whose reference cDNA position has an open side (`?`).
    Both notations lie within the allele's span over the transcript, Ensembl VEP's at most one
    base before it and an `ins` of either at most one base after it, so every coordinate of
    vep-rs's lies within the span's length plus one of every coordinate of Ensembl VEP's; the two
    name the same transcript and numbering; against an empty reference vep-rs's is a `del`, and
    against an empty vep-rs notation Ensembl VEP's is an `ins` or a `dup`, what the outside bases
    read as once the clip has stripped the unchanged bases inside, or, on the row of an allele
    shorter than its span, a `del`, a `delins` or a substitution: vep-rs describes no deletion
    that runs over the transcript's edge and can be written nowhere wholly inside it. The kinds
    and bases are otherwise unconstrained; the manifest records vep-rs's exact notation, which
    the golden harness compares."""
    location, allele, feature, _feature_type, _consequences = key
    if allele == "-":
        return False
    position = (perl_row or {}).get("cDNA_position", "")
    start_text, _, end_text = position.partition("/")[0].partition("-")
    if "?" not in (start_text, end_text):
        return False
    reference = clipped_notation(vep, feature)
    candidate = clipped_notation(vep_rs, feature)
    if vep in ("", "-"):
        return candidate is not None and candidate.group("plain") == "del"
    if reference is None:
        return False
    span_start, _, span_end = location.split(":", 1)[-1].partition("-")
    span_bases = int(span_end) - int(span_start) + 1 if span_end.isdigit() else 1
    if vep_rs in ("", "-"):
        if reference.group("kind") == "ins" or reference.group("plain") == "dup":
            return True
        deletes = reference.group("plain") == "del" or reference.group("kind") == "delins" or reference.group("alternate_base") is not None
        return deletes and 0 < len(allele) < span_bases
    if candidate is None or vep.partition(":")[0] != vep_rs.partition(":")[0]:
        return False
    if candidate.group("numbering") != reference.group("numbering"):
        return False
    reach = span_bases + 1
    return all(
        abs(c - r) <= reach
        for c in notation_coordinate_indexes(candidate)
        for r in notation_coordinate_indexes(reference)
    )


def insertion_without_reference_notation(key: tuple[str, str, str, str, str], candidate: re.Match,
                                         perl_row: dict[str, str]) -> bool:
    """Whether the row is one on which vep-rs prints an `ins` where Ensembl VEP prints no
    `HGVSc`: a sequence allele's row on a reverse-strand transcript (`STRAND` `-1`) whose
    reference cDNA position is open at its start (`?-N`), the transcript's first base lying at
    its genomic end, vep-rs's notation an `ins` of the bases a longer allele adds inside the
    span (`insertion_inside_the_span`). There the clip leaves Ensembl VEP an insertion at slice
    position 0, for which `_get_cDNA_position` returns no coordinate
    (TranscriptVariationAllele.pm:2834), so its row carries no notation and the HGVSIntronOffset
    plugin writes no field beside it. The test places vep-rs's two fields; the row's `HGVSc` is
    no pair of the class. A forward-strand transcript's row and a position open at its end are
    refused."""
    location, allele, _feature, _feature_type, _consequences = key
    start_text, _, end_text = perl_row.get("cDNA_position", "").partition("/")[0].partition("-")
    if allele == "-" or start_text != "?" or perl_row.get("STRAND", "") != "-1":
        return False
    return insertion_inside_the_span(location, allele, candidate, end_text)


def insertion_inside_the_span(location: str, allele: str, candidate: re.Match, closed: str) -> bool:
    """Whether vep-rs's `ins` notation (an `HGVSC_CLIPPED_BODY` match) writes the bases a longer
    allele adds inside the span of a row whose cDNA position is open at its start: the allele
    is a sequence longer than its Location, the notation inserts as many bases as it is longer,
    between two exon coordinates next to each other, and under `n.` numbering, whose
    coordinates are cDNA positions, at or beside a base inside: its first coordinate at most
    `closed`, the position's closed end. Under `c.` numbering the coordinates' place is not
    read."""
    first, last = candidate.group("start"), candidate.group("end")
    offset = re.compile(r"\d[+-]\d+$")
    if candidate.group("kind") != "ins" or last is None or any(offset.search(c) or c.startswith("*-") for c in (first, last)):
        return False
    if cdna_coordinate_index(last) - cdna_coordinate_index(first) != 1 or re.fullmatch(r"[ACGTN]+", allele) is None:
        return False
    span_start, _, span_end = location.split(":", 1)[-1].partition("-")
    if not span_start.isdigit() or (span_end and not span_end.isdigit()):
        return False
    span_bases = int(span_end) - int(span_start) + 1 if span_end else 1
    if len(candidate.group("bases")) != len(allele) - span_bases:
        return False
    if candidate.group("numbering") != "n":
        return True
    if not closed.isdigit() or not first.isdigit() or not last.isdigit():
        return False
    return int(first) <= int(closed)


def notation_intron_offsets(notation: re.Match) -> tuple[str, str] | None:
    """The two values the HGVSIntronOffset plugin writes beside a notation (an
    `HGVSC_CLIPPED_BODY` match), the start's and the end's: the intron parts of the two
    coordinates as `hgvs_transcript` reads them (TranscriptVariationAllele.pm:1530-1533, an
    absent part read as 0) before it orders the pair for printing (:1535-1538). The routine
    holds the `ins` that `_clip_alleles` leaves, its start one above its end (:2245-2248),
    from its 3' flank to its 5' flank and turns the two unless the end holds a star, so the
    start is the coordinate printed second, or the one printed first where the 5' flank lies
    past the stop codon; it holds a `dup` and every other kind in its printed order, one
    coordinate standing for both. `None` for a pair the routine would have printed the other way round, an `ins` of
    one coordinate, and a coordinate before the stop codon's last base (`*-5`)."""
    first, last = notation.group("start"), notation.group("end")
    insertion = notation.group("kind") == "ins"
    if (insertion and last is None) or any(c.startswith("*-") for c in (first, last or first)):
        return None
    first_parts, last_parts = (
        tuple(int(part or 0) for part in COORDINATE_PARTS_AS_READ.search(c).groups()) for c in (first, last or first)
    )
    if insertion and last_parts > first_parts and "*" not in first:
        return str(last_parts[1]), str(first_parts[1])
    if first_parts > last_parts and "*" not in (last or first):
        return None
    return str(first_parts[1]), str(last_parts[1])


def insertion_before_the_first_base_against_duplication(key: tuple[str, str, str, str, str], reference: re.Match | None,
                                                        candidate: re.Match, pairs: list[tuple[str, ...]],
                                                        perl_row: dict[str, str], rust_row: dict[str, str]) -> bool:
    """Whether the two rows are Ensembl VEP's `ins` from one base before a forward-strand
    transcript's first base against vep-rs's duplication of the bases the allele adds, the
    fourth shape of `clipped_notation_intron_offsets`, read from the two rows alone in four
    parts, all of which hold. The row: its allele holds more bases than its Location spans, `A`
    bases over `S` positions; its `STRAND` is `1`; its reference cDNA position is open at its
    start and closed at its end, `?-N`, with `N` above 0 and below `S`. Ensembl VEP's side: an
    `ins` of the row's transcript whose two coordinates are one coordinate `X` carrying the
    intron part `-1` and the same `X` carrying none, writing the allele's bases less the last
    `N`, with no `HGVS_OFFSET` and the plugin's pair 0 and -1: the clip has consumed the `N`
    bases the clamp at the first base left it (TranscriptVariationAllele.pm:2231-2242, :2271,
    :2732-2733), so `X` is the first base's coordinate. vep-rs's side: a `dup` of the same
    transcript and numbering whose coordinates carry no intron part and cover `A` - `S` bases
    (`exonic_run`), an `HGVS_OFFSET` that is a whole number `k`, an empty one read as 0, and
    the plugin's pair 0 and 0. The tie: from `X` to the duplication's last coordinate are
    `N` + `k` exonic bases of one region. The `ins` the clip leaves at a true intron, which
    writes another count of bases, a duplication across the start codon from `X`, a row
    without `STRAND`, a reverse-strand transcript's row and a pair either of whose rows
    carries `BAM_EDIT`, on which `hgvs_transcript` respells the notation's alleles after the
    clip (:1510-1514), are refused."""
    location, allele, _feature, _feature_type, _consequences = key
    blank = ("", "-")
    if reference is None or re.fullmatch(r"[ACGT]+", allele) is None:
        return False
    if perl_row.get("STRAND", "") != "1" or any(row.get("BAM_EDIT", "") not in blank for row in (perl_row, rust_row)):
        return False
    span_start, _, span_end = location.split(":", 1)[-1].partition("-")
    span_end = span_end or span_start
    if not (span_start.isdigit() and span_end.isdigit()):
        return False
    span_bases = int(span_end) - int(span_start) + 1
    added = len(allele) - span_bases
    start_text, _, end_text = perl_row.get("cDNA_position", "").partition("/")[0].partition("-")
    if added < 1 or start_text != "?" or not end_text.isdigit() or not 0 < int(end_text) < span_bases:
        return False
    inside = int(end_text)
    first, last = reference.group("start"), reference.group("end")
    if reference.group("kind") != "ins" or last is None:
        return False
    first_base = next((m.group(1) for c, o in ((first, last), (last, first))
                       if (m := ONE_BEFORE_COORDINATE.fullmatch(c)) and m.group(1) == o), None)
    if first_base is None or reference.group("bases") != allele[: len(allele) - inside]:
        return False
    if perl_row.get("HGVS_OFFSET", "") not in blank or pairs[0] != ("0", "-1"):
        return False
    if candidate.group("plain") != "dup" or candidate.group("numbering") != reference.group("numbering"):
        return False
    duplicated_last = candidate.group("end") or candidate.group("start")
    if exonic_run(candidate.group("start"), duplicated_last) != added or pairs[1] != ("0", "0"):
        return False
    carried = rust_row.get("HGVS_OFFSET", "")
    if carried not in blank and re.fullmatch(r"[0-9]+", carried) is None:
        return False
    return exonic_run(first_base, duplicated_last) == inside + (0 if carried in blank else int(carried))


def clipped_notation_intron_offsets(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                    perl_row: dict[str, str] | None, rust_row: dict[str, str] | None) -> bool:
    """Whether a differing HGVSIntronOffset plugin field is one
    `hgvs_clipped_reference_against_whole_alternate` takes, in four shapes: Ensembl VEP alone
    writes the plugin's two fields, on a row whose `HGVSc` is the `ins` or `dup`, or the `del`, `delins` or
    substitution of an allele shorter than its span, that the class takes against a vep-rs row
    that prints none (`clipped_reference_against_whole_alternate`);
    vep-rs alone writes them, on a row whose Ensembl VEP `HGVSc` is empty where vep-rs prints an
    `ins` inside the span at a reverse-strand transcript's first base
    (`insertion_without_reference_notation`), the `HGVSc` of such a row staying outside the
    class; or both write them, on a row whose two notations the class takes as a pair, each an
    `ins`, and the field named differs. In each of these three a
    side's two values are the intron parts of its own notation's coordinates as the plugin
    reads them (`notation_intron_offsets`), and beside vep-rs's fields neither row carries an
    `HGVS_OFFSET`. Or both write them and the end's field alone is named, on a forward-strand
    transcript's row whose Ensembl VEP `HGVSc` is an `ins` from one base before the
    transcript's first base, beside the pair 0 and -1, where vep-rs prints the duplication of
    the bases the allele adds, beside the pair 0 and 0 and the bases the HGVS 3' rule carried
    it in `HGVS_OFFSET` (`insertion_before_the_first_base_against_duplication`); the
    `HGVS_OFFSET` of such a row stays outside the class. The row's `HGVSc`, its `HGVS_OFFSET`
    and the plugin's other field are read from `perl_row` and `rust_row`, so a run that does
    not compare them names no row."""
    _location, _allele, feature, _feature_type, _consequences = key
    if name not in HGVS_INTRON_OFFSET_FIELDS:
        return False
    perl_row, rust_row = perl_row or {}, rust_row or {}
    blank = ("", "-")
    pairs = [tuple(row.get(f, "") for f in HGVS_INTRON_OFFSET_FIELDS) for row in (perl_row, rust_row)]
    filled = [all(value not in blank for value in pair) for pair in pairs]
    empty = [all(value in blank for value in pair) for pair in pairs]
    hgvsc, candidate_hgvsc = perl_row.get("HGVSc", ""), rust_row.get("HGVSc", "")
    reference, candidate = clipped_notation(hgvsc, feature), clipped_notation(candidate_hgvsc, feature)
    if filled[0] and empty[1]:
        if vep in blank or vep_rs not in blank or perl_row.get(name) != vep or hgvsc in blank or candidate_hgvsc not in blank:
            return False
        if not clipped_reference_against_whole_alternate(key, hgvsc, "", perl_row):
            return False
        offsets = notation_intron_offsets(reference)
        return offsets is not None and offsets == pairs[0]
    if not filled[1] or vep_rs in blank or rust_row.get(name) != vep_rs:
        return False
    if candidate is not None and candidate.group("plain") == "dup":
        if name != HGVS_INTRON_OFFSET_FIELDS[1] or perl_row.get(name) != vep:
            return False
        return insertion_before_the_first_base_against_duplication(key, reference, candidate, pairs, perl_row, rust_row)
    if candidate is None or candidate.group("kind") != "ins":
        return False
    if any(row.get("HGVS_OFFSET", "") not in blank for row in (perl_row, rust_row)):
        return False
    if notation_intron_offsets(candidate) != pairs[1]:
        return False
    if empty[0]:
        return vep in blank and hgvsc in blank and insertion_without_reference_notation(key, candidate, perl_row)
    if not clipped_reference_against_whole_alternate(key, hgvsc, candidate_hgvsc, perl_row):
        return False
    if not filled[0] or vep in blank or perl_row.get(name) != vep or vep == vep_rs or reference is None:
        return False
    if reference.group("kind") != "ins":
        return False
    return notation_intron_offsets(reference) == pairs[0]


def notation_coordinate_indexes(notation: re.Match) -> list[int]:
    """The `cdna_coordinate_index` of an `HGVSC_CLIPPED_BODY` match's one or two coordinates."""
    return [cdna_coordinate_index(c) for c in (notation.group("start"), notation.group("end")) if c]


def cdna_coordinate_index(coordinate: str) -> int:
    """An HGVS cDNA coordinate as a position on one integer scale along the transcript: the 5'
    UTR `-N` at `-N`, the coding sequence `N` at `N - 1` (`c.-1` and `c.1` are neighbours), the 3'
    UTR `*N` a billion further, since its distance from the coding sequence is the coding length
    the notation does not carry, each with its intron offset added (`-49-1` is one base before
    `-49`, `83+1` one after `83`)."""
    region, position, offset = cdna_coordinate_order(coordinate)
    if region < 0:
        return position + offset
    if region == 0:
        return position - 1 + offset
    return 10**9 + position + offset


def exonic_run(first: str, last: str) -> int | None:
    """The count of bases from the cDNA coordinate `first` to `last`, both exonic and in one
    region of the transcript (5' UTR, coding sequence or 3' UTR), `None` where either carries an
    intron offset or the two lie in different regions: the scale `cdna_coordinate_index` reads
    carries no distance across the stop codon or an intron."""
    (first_region, _, first_offset), (last_region, _, last_offset) = cdna_coordinate_order(first), cdna_coordinate_order(last)
    if first_region != last_region or first_offset or last_offset:
        return None
    return cdna_coordinate_index(last) - cdna_coordinate_index(first) + 1


def clipped_notation(hgvsc: str, feature: str) -> re.Match | None:
    """The `HGVSC_CLIPPED_BODY` match of an HGVSc that names the row's transcript (its stable id,
    with or without a version), or `None` for an empty value, another transcript's notation or a
    notation of another shape."""
    if hgvsc in ("", "-"):
        return None
    name, _, body = hgvsc.partition(":")
    if name.split(".")[0] != feature.split(".")[0]:
        return None
    return HGVSC_CLIPPED_BODY.fullmatch(body)


def added_bases_written_outside_the_transcript(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                               strand: str, perl_row: dict[str, str] | None,
                                               rust_row: dict[str, str] | None) -> bool:
    """Whether a differing `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset plugin pair is vep-rs's
    notation of inserted bases on a row both engines annotate outside the transcript, its one
    term `upstream_gene_variant` or `downstream_gene_variant`: Ensembl VEP fills the HGVSc
    column only on a row within the feature (OutputFactory.pm:1721), so it prints no `HGVSc` and
    no `HGVS_OFFSET` there and writes the plugin's two fields or neither
    (HGVSIntronOffset.pm:62), where vep-rs writes the `dup` or `ins` the HGVS 3' rule places on
    the row's transcript. The allele is an insertion's (a Location of two neighbouring
    positions) or holds more bases than its Location spans, and the notation is of as many bases
    as the allele holds or adds to its span where its coordinates count them; a downstream row's
    is a `dup` beside no `HGVS_OFFSET`, an upstream row's stands beside a non-zero `HGVS_OFFSET`,
    negative on a reverse strand (`strand`), larger in size than the reference row's `DISTANCE`
    on an insertion's row and no smaller than it on the row of an allele that holds more bases
    than its span, whose `DISTANCE` counts from that span,
    and vep-rs's row carries both plugin fields, the intron offsets of the notation's
    coordinates in their printed order, or neither where Ensembl VEP's row carries neither, as
    in a run without the plugin. The row's other fields are read from `perl_row` and `rust_row`,
    so a run that does not compare `HGVSc` names no row, one that does not compare
    `HGVS_OFFSET` no upstream row, and one that does not compare the plugin's fields reads as a
    run without the plugin."""
    location, allele, feature, _feature_type, consequences = key
    if consequences not in OUTSIDE_TRANSCRIPT_TERMS or re.fullmatch(r"[ACGT]+", allele) is None:
        return False
    span_start, _, span_end = location.split(":", 1)[-1].partition("-")
    span_end = span_end or span_start
    if not (span_start.isdigit() and span_end.isdigit()):
        return False
    positions = int(span_end) - int(span_start) + 1
    counts = ({len(allele)} if positions == 2 else set()) | ({len(allele) - positions} if len(allele) > positions else set())
    perl_row, rust_row = perl_row or {}, rust_row or {}
    if not counts or any(perl_row.get(f, "") not in ("", "-") for f in ("HGVSc", "HGVS_OFFSET")):
        return False
    reference_offsets = [perl_row.get(f, "") for f in HGVS_INTRON_OFFSET_FIELDS]
    if (reference_offsets[0] in ("", "-")) != (reference_offsets[1] in ("", "-")):
        return False
    notation = clipped_notation(vep_rs if name == "HGVSc" else rust_row.get("HGVSc", ""), feature)
    if notation is None or not (notation.group("plain") == "dup" or notation.group("kind") == "ins"):
        return False
    first, last = notation.group("start"), notation.group("end")
    if notation.group("kind") == "ins":
        if last is None:
            return False
        counted = len(notation.group("bases"))
    else:
        counted = exonic_run(first, last or first)
    if counted is not None and counted not in counts:
        return False
    if any(c.startswith("*-") for c in (first, last or first)):
        return False
    written = tuple(rust_row.get(f, "") for f in HGVS_INTRON_OFFSET_FIELDS)
    # A run without the plugin writes the two fields on neither row.
    without_plugin = all(v in ("", "-") for v in (*reference_offsets, *written))
    if not without_plugin and written != tuple(str(cdna_coordinate_order(c)[2]) for c in (first, last or first)):
        return False
    offset = vep_rs if name == "HGVS_OFFSET" else rust_row.get("HGVS_OFFSET", "")
    if consequences == "downstream_gene_variant":
        if offset not in ("", "-") or notation.group("plain") != "dup":
            return False
    else:
        distance = perl_row.get("DISTANCE", "")
        if NON_ZERO_INTEGER.fullmatch(offset) is None or not distance.isdigit():
            return False
        # An insertion's DISTANCE counts from the nearer base beside its point, a longer allele's from its span.
        added_to_a_span = len(allele) > positions and (positions != 2 or counted == len(allele) - positions)
        if int(offset.lstrip("-")) < int(distance) + (0 if added_to_a_span else 1):
            return False
        if strand not in ("", "-") and offset.startswith("-") != (strand == "-1"):
            return False
    if name in HGVS_INTRON_OFFSET_FIELDS:
        return vep_rs not in ("", "-") and rust_row.get(name) == vep_rs
    return name in ("HGVSc", "HGVS_OFFSET") and vep in ("", "-")


def stop_codon_protein_comparison(key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                  perl_row: dict[str, str] | None) -> bool:
    """`hgvsp_stop_codon_protein_comparison` on an `HGVSp` the engines write differently, as the
    format comparator's class of that name reads it: both values name one protein, the
    reference row's `FLAGS` (`perl_row`) name no `cds_end_NF`, and with an escaped `=` set aside
    Ensembl VEP's description is the stop replaced by a stop against vep-rs's stop unchanged at
    the same residue above 1; one residue and a stop inserted between two neighbouring
    residues, the residue the second neighbour's, on a row whose terms hold
    `stop_retained_variant`, against vep-rs's stop unchanged at the residue after that
    neighbour; a range from a residue to the stop replaced by residues, none of
    them a stop, on a row whose terms hold `stop_lost`, against vep-rs's frameshift named at
    that residue with the first of them, its count unknown, or above their number and not the
    one that leaves the new stop at the stop's own residue; or, against vep-rs's extension
    named at the stop: the stop alone replaced by the first residues that follow the protein
    (the extension's first residue the first named, its count unknown, above the residues
    named, or their number where no `Ter` closes them), the stop replaced by every added
    residue and then a stop (the count their number), residues inserted between the last
    residue and the stop (the count their number), the stop unchanged on a row whose terms
    hold `stop_lost`, the same extension with an unknown count on such a row, or the stop
    written as deleted before the extension's own count on such a row."""
    if vep in ("", "-") or vep_rs in ("", "-") or "cds_end_NF" in (perl_row or {}).get("FLAGS", ""):
        return False
    vep_name, separator, reference = vep.partition(":p.")
    vep_rs_name, _, candidate = vep_rs.partition(":p.")
    if not separator or vep_name != vep_rs_name:
        return False
    reference, candidate = reference.replace("%3D", "="), candidate.replace("%3D", "=")
    replaced, stays = STOP_REPLACED.fullmatch(reference), STOP_UNCHANGED.fullmatch(candidate)
    # The stop replaced by a stop.
    if replaced and stays and replaced.group(2) == "Ter" and replaced.group(1) == stays.group(1) != "1":
        return True

    # One residue and a stop inserted between two neighbouring residues, against the stop unchanged.
    placed = RESIDUE_AND_STOP_INSERTED.fullmatch(reference)
    if placed and stays:
        left, before, right, after, residue = placed.groups()
        return ("Ter" not in (left, right) and residue == right and int(before) + 1 == int(after)
                and int(stays.group(1)) == int(after) + 1 and "stop_retained_variant" in key[4].split(","))

    # A range from a residue to the stop replaced, against a frameshift named at that residue.
    ranged, frameshift = STOP_RANGE_REPLACED.fullmatch(reference), STOP_FRAMESHIFT.fullmatch(candidate)
    if ranged and frameshift:
        residue, at, stop, residues = ranged.groups()
        new = [residues[i:i + 3] for i in range(0, len(residues), 3)]
        shifted = frameshift.group(4)
        return (residue != "Ter" and "Ter" not in new and int(at) < int(stop)
                and frameshift.group(1) == residue and frameshift.group(2) == at and frameshift.group(3) == new[0]
                and (shifted == "?" or (int(shifted) > len(new) and int(shifted) != int(stop) - int(at) + 1))
                and "stop_lost" in key[4].split(","))

    extension = STOP_EXTENSION.fullmatch(candidate)
    if extension is None or extension.group(1) == "1" or extension.group(2) == "Ter":
        return False
    number, first, count = extension.groups()
    # An extension with an unknown count against the same extension with a count.
    unknown = STOP_EXTENSION.fullmatch(reference)
    if unknown and unknown.groups() == (number, first, "?") and count != "?" and "stop_lost" in key[4].split(","):
        return True
    # The stop written as deleted with an extension count, against an extension with that count.
    deleted = STOP_DELETED_EXTENSION.fullmatch(reference)
    if deleted and deleted.groups() == (number, count) and "stop_lost" in key[4].split(","):
        return True
    named = [replaced.group(2)[i:i + 3] for i in range(0, len(replaced.group(2)), 3)] if replaced else []
    closed = bool(named) and named[-1] == "Ter"
    named = named[:-1] if closed else named
    opens = bool(named) and "Ter" not in named and replaced.group(1) == number and named[0] == first
    # The stop alone replaced by the first residues that follow the protein.
    if opens and (count == "?" or int(count) > len(named) or (int(count) == len(named) and not closed)):
        return True

    # The stop replaced by every added residue and then the stop the inserted bases bring.
    if opens and closed and count == str(len(named)):
        return True

    # Residues inserted between the last residue and the stop.
    inserted = STOP_INSERTION.fullmatch(reference)
    if inserted:
        between = [inserted.group(3)[i:i + 3] for i in range(0, len(inserted.group(3)), 3)]
        if (int(inserted.group(1)) + 1 == int(inserted.group(2)) and inserted.group(2) == number
                and "Ter" not in between and between[0] == first and count == str(len(between))):
            return True

    # The stop written as unchanged where the row's terms say it is lost.
    unchanged = STOP_UNCHANGED.fullmatch(reference)
    return unchanged is not None and unchanged.group(1) == number and "stop_lost" in key[4].split(",")


def partial_terminal_codon_models(cache_dir: Path) -> frozenset[str]:
    """The cache's transcripts whose coding sequence is not a whole number of codons, by stable
    id: a coding sequence that ends inside a codon, which no stop codon closes. Empty when the
    cache is absent."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc

    partial: set[str] = set()
    for _chrom, tx in bmc.iter_cache_transcripts(cache_dir):
        coding = (tx.get("variation_effect_feature_cache") or {}).get("translateable_seq") or ""
        if len(coding) % 3:
            partial.add(str(tx.get("stable_id", "")))
    partial.discard("")
    return frozenset(partial)


def stop_named_on_incomplete_coding_sequence(key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                                             perl_row: dict[str, str] | None,
                                             partial_codon_models: frozenset[str] | None = None) -> bool:
    """`hgvsp_stop_named_on_incomplete_coding_sequence` on an `HGVSp` the engines write
    differently, as the format comparator's class of that name reads it: vep-rs writes none;
    Ensembl VEP's description, an escaped `=` and a pair of prediction parentheses set aside,
    names the stop alone as unchanged, deleted, replaced by one residue or extended; and the row
    shows the incomplete coding sequence: the reference row's `FLAGS` (`perl_row`) name
    `cds_end_NF`, its terms hold `incomplete_terminal_codon_variant`, or the description is
    `Ter1=` and its `FLAGS` name `cds_start_NF`. Beyond the row, the corpus cache shows it: the
    row's transcript is one of `partial_codon_models` (`partial_terminal_codon_models`)."""
    if vep in ("", "-") or vep_rs not in ("", "-"):
        return False
    _, separator, reference = vep.partition(":p.")
    reference = reference.replace("%3D", "=")
    if reference.startswith("(") and reference.endswith(")"):
        reference = reference[1:-1]
    named = STOP_NAMED_ALONE.fullmatch(reference)
    if not separator or named is None:
        return False
    flags = (perl_row or {}).get("FLAGS", "")
    return ("cds_end_NF" in flags or "incomplete_terminal_codon_variant" in key[4].split(",")
            or (named.group(1) == "1" and named.group(2) == "=" and "cds_start_NF" in flags)
            or key[2] in (partial_codon_models or ()))


def copy_of_the_intron_last_base(body: str) -> bool:
    """The copy form of `hgvsp_intron_end_insertion`: an HGVSc body duplicating the one position
    before an exon's first coding base `N`, `c.N-1dup`."""
    return INTRON_END_COPY_BODY.fullmatch(body) is not None


def duplication_of_the_intron_last_bases(body: str) -> bool:
    """The duplication form of `hgvsp_intron_end_insertion`: an HGVSc body duplicating two or more
    intron positions that end on the one before an exon's first coding base `N`, its first
    position counted back from `N` (`c.N-K_N-1dup`, `K` at least 2) or, in the intron's first
    half, on from the coding base before `N` (`c.M+J_N-1dup`, `M` one below `N`)."""
    m = INTRON_END_DUPLICATION_BODY.fullmatch(body)
    if m is None:
        return False
    first, sign, distance, last = m.groups()
    if sign == "-":
        return first == last and int(distance) >= 2
    return int(first) == int(last) - 1


def insertion_that_cannot_end_the_intron(body: str) -> bool:
    """The insertion form of `hgvsp_intron_end_insertion`: an HGVSc body inserting, between the
    position before an exon's first coding base `N` and that base (`c.N-1_Nins`), one base, or
    bases that do not end in one of `INTRON_END_PAIRS`."""
    m = INTRON_END_INSERTION_BODY.fullmatch(body)
    if m is None or m.group(1) != m.group(2):
        return False
    bases = m.group(3)
    return len(bases) == 1 or bases[-2:] not in INTRON_END_PAIRS


def insertion_that_ends_in_an_intron_end_pair(body: str) -> bool:
    """The pair form of `hgvsp_intron_end_insertion`: an HGVSc body inserting, between the
    position before an exon's first coding base `N` and that base (`c.N-1_Nins`), two or more
    bases that end in one of `INTRON_END_PAIRS`."""
    m = INTRON_END_INSERTION_BODY.fullmatch(body)
    if m is None or m.group(1) != m.group(2):
        return False
    bases = m.group(3)
    return len(bases) >= 2 and bases[-2:] in INTRON_END_PAIRS


def intron_end_insertion_form(body: str) -> str | None:
    """The form of `hgvsp_intron_end_insertion` an HGVSc body has, each form one helper and one
    block here, or `None` for any other notation."""
    if copy_of_the_intron_last_base(body):
        return "copy"

    if duplication_of_the_intron_last_bases(body):
        return "duplication"

    if insertion_that_cannot_end_the_intron(body):
        return "insertion"

    if insertion_that_ends_in_an_intron_end_pair(body):
        return "pair"

    return None


def intron_end_insertion(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                         perl_row: dict[str, str] | None, rust_row: dict[str, str] | None,
                         intron_ends: dict[tuple[str, str], frozenset[int]] | None) -> bool:
    """Whether a differing `HGVSp` pair is Ensembl VEP's notation against none from vep-rs on an
    insertion the input writes between an intron's last base and the first coding base of the
    next exon: Ensembl VEP reads the inserted bases as the exon's first bases and writes the
    notation of that reading, and vep-rs writes none, because the HGVS nomenclature gives a
    change at an intron's end no predicted protein consequence. The two rows carry one `HGVSc`,
    a notation of the row's transcript of one of the four forms (`intron_end_insertion_form`),
    and no `HGVS_OFFSET`, both read from `perl_row` and `rust_row`, so a run that does not
    compare the two fields beside `HGVSp` names nothing; the row's Location names two adjacent
    positions, the lower of them a gap of the row's transcript in `intron_ends`
    (`intron_end_gaps`), under the Location's chromosome as written or without a leading `chr`
    and the row's Feature as written or without a trailing version. A record that holds the
    insertion with reference bases it shares is placed by the bases an `ins` notation inserts
    (`location_in_gaps`, `held_insertion_gaps`), and a row that gives such a record's span with a
    minimised allele by those bases or by the count of positions a `dup` notation duplicates
    (`lone_insertion_gaps`)."""
    location, allele, feature, _feature_type, _consequences = key
    if name != "HGVSp" or ":p." not in vep or vep_rs not in ("", "-") or not intron_ends:
        return False
    perl_row, rust_row = perl_row or {}, rust_row or {}
    hgvsc = perl_row.get("HGVSc", "")
    if hgvsc in ("", "-") or rust_row.get("HGVSc") != hgvsc:
        return False
    if perl_row.get("HGVS_OFFSET") not in ("", "-") or rust_row.get("HGVS_OFFSET") not in ("", "-"):
        return False
    transcript, _, body = hgvsc.partition(":")
    stable_id = feature.split(".")[0]
    if transcript.split(".")[0] != stable_id or intron_end_insertion_form(body) is None:
        return False
    row = (allele.upper(), notation_inserted_bases(INTRON_END_INSERTION_BODY, body), notation_duplicated_positions(body),
           perl_row.get("STRAND", ""))
    return location_in_gaps(location, feature, intron_ends, row=row)


def notation_inserted_bases(pattern: re.Pattern, body: str) -> str:
    """The bases an HGVSc body inserts where it is an `ins` notation of `pattern`
    (`INTRON_END_INSERTION_BODY`, `INTRON_START_INSERTION_BODY`); empty for any other body."""
    m = pattern.fullmatch(body)
    return m.group(3) if m is not None else ""


def reverse_complement(bases: str) -> str:
    """`bases` as the other strand reads them."""
    return bases.translate(BASE_COMPLEMENT)[::-1]


def held_insertion_gaps(first: int, length: int, allele: str, inserted: str, strand: str) -> list[int]:
    """The gaps in which an insertion can lie that a row's record holds with the reference bases
    it shares, each the lower of the two positions either side of the gap. The record's reference
    spans `length` positions from `first` and `allele` is its alternate allele on the forward
    strand, those bases with the inserted ones; `inserted` are the bases the row's HGVSc inserts,
    on the transcript's strand, read on the forward strand by `strand` (`1`, `-1`; any other value
    admits either reading). The allele holds as many bases more than the span as the notation
    inserts. For each count of span bases clipped from the allele's front, the rest clipped from
    its back, the bases between must be the inserted ones, and the gap's lower position is `first`
    plus that count less one. A row states no reference bases, so every clipping that leaves the
    inserted bases names a gap."""
    if not inserted or length < 1 or len(allele) - length != len(inserted):
        return []
    if strand == "1":
        readings = {inserted}
    elif strand == "-1":
        readings = {reverse_complement(inserted)}
    else:
        readings = {inserted, reverse_complement(inserted)}
    gaps = []
    for front in range(length + 1):
        back = len(allele) - (length - front)
        if allele[front:back] in readings:
            gaps.append(first + front - 1)
    return gaps


def notation_duplicated_positions(body: str) -> int:
    """The count of positions an HGVSc body duplicates where the notation states it: one for a
    copy of the intron's last base (`c.N-1dup`), `K` for a duplication counted back from the
    exon's first coding base (`c.N-K_N-1dup`). 0 for any other body, a duplication that starts in
    the intron's first half among them, whose count depends on the intron's length."""
    if INTRON_END_COPY_BODY.fullmatch(body) is not None:
        return 1
    m = INTRON_END_DUPLICATION_BODY.fullmatch(body)
    if m is not None and m.group(2) == "-" and m.group(1) == m.group(4):
        return int(m.group(3))
    return 0


def lone_insertion_gaps(first: int, length: int, allele: str, inserted: str, duplicated: int = 0, moved: int = 0) -> list[int]:
    """The gaps a text row's Location names where its Allele is the insertion alone, each the
    lower of the two positions either side of the gap: every gap inside and beside the span of
    `length` positions from `first`, from `first` less one to the span's last position. A record
    that carries another allele beside the insertion keeps its reference span as the row's
    Location, and a minimised allele spells the inserted bases alone. `allele` is the insertion
    the row's HGVSc describes: the bases an `ins` notation inserts (`inserted`), on either strand,
    as they stand or turned back by the bases the 3' rule moved the record (`moved`), since a run
    prints a moved record's allele as the rule moved it or as the record writes it; or, for a
    `dup` notation, which names no bases, as many bases as the positions it duplicates
    (`duplicated`) on an unmoved record. None for any other allele."""
    if length < 1 or any(ord(base) not in BASE_COMPLEMENT for base in allele):
        return []
    if inserted:
        turns = dict.fromkeys((0, abs(moved) % len(inserted)))
        readings = (allele, reverse_complement(allele))
        if all(bases[turn:] + bases[:turn] != inserted for bases in readings for turn in turns):
            return []
    elif duplicated < 1 or moved != 0 or len(allele) != duplicated:
        return []
    return list(range(first - 1, first + length))


def location_in_gaps(location: str, feature: str, gaps: dict[tuple[str, str], frozenset[int]], moved: int = 0,
                     row: tuple[str, str, int, str] | None = None) -> bool:
    """Whether a Location names two adjacent positions the lower of which, moved `moved` positions
    along the chromosome or as named, is a gap of the transcript `feature` in `gaps`, under the
    Location's chromosome as written or without a leading `chr` and the Feature as written or
    without a trailing version: a run that prints the moved location names the gap the insertion
    is left in. With `row`, the row's Allele, the bases its HGVSc inserts, the count of
    positions it duplicates and its STRAND, the Location is read as well as the reference span of
    a record that carries another allele beside the insertion: an unmoved Allele that holds the
    span's bases with the inserted ones, as written or reverse-complemented, since the row spells
    it on the record's strand, names the gaps its clippings leave (`held_insertion_gaps`); and
    an Allele that is the insertion alone, a minimised allele, names every gap of the span, read
    as named or moved alike (`lone_insertion_gaps`), a span of two adjacent positions among
    them, which is an insertion written alone or a reference allele of two bases."""
    chrom, _, span = location.rpartition(":")
    start, _, end = span.partition("-")
    end = end or start
    if not (start.isdigit() and end.isdigit()):
        return False
    first, length = int(start), int(end) - int(start) + 1
    positions = [first] if length == 2 else []
    if row is not None:
        allele, inserted, duplicated, strand = row
        positions += lone_insertion_gaps(first, length, allele, inserted, duplicated, moved)
        if moved == 0:
            for bases in (allele, reverse_complement(allele)):
                positions += held_insertion_gaps(first, length, bases, inserted, strand)
    names = (feature, re.sub(r"\.\d+$", "", feature))
    return any(position + move in gaps.get((contig, name), ())
               for position in positions for move in dict.fromkeys((moved, 0))
               for contig in (chrom, chrom.removeprefix("chr")) for name in names)


def insertion_after_the_exon_last_base(body: str) -> bool:
    """The form of `hgvsp_intron_start_insertion`: an HGVSc body inserting bases between an exon's
    last coding base `N` and the first position of the intron after it, `c.N_N+1ins`."""
    m = INTRON_START_INSERTION_BODY.fullmatch(body)
    return m is not None and m.group(1) == m.group(2)


def intron_start_insertion(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                           perl_row: dict[str, str] | None, rust_row: dict[str, str] | None,
                           intron_starts: dict[tuple[str, str], frozenset[int]] | None) -> bool:
    """Whether a differing `HGVSp` pair is Ensembl VEP's notation against none from vep-rs on an
    insertion the 3' rule leaves between a coding exon's last base and the first base of the
    intron after it: Ensembl VEP reads the inserted bases as the exon's last bases and writes the
    notation of that reading, and vep-rs writes none, because the HGVS nomenclature gives a
    change at an intron's start no predicted protein consequence. The two rows carry one `HGVSc`,
    a notation of the row's transcript of the class's form (`insertion_after_the_exon_last_base`),
    and either no `HGVS_OFFSET` or one non-zero integer, the bases the 3' rule moved the record,
    all read from `perl_row` and `rust_row`, so a run that does not compare the two fields beside
    `HGVSp` names nothing; the lower of the two adjacent positions the row's Location names,
    moved by the offset as signed, is a gap of the row's transcript in `intron_starts`
    (`intron_start_gaps`, `location_in_gaps`). A row that carries the offset is the class as
    well where its Location names such a gap unmoved: a run under `--shift_genomic 1` prints the
    Location the 3' rule carried the record to. A record that holds the insertion with reference
    bases it shares carries no offset and is placed by the bases the notation inserts
    (`held_insertion_gaps`), and a row that gives such a record's span with a minimised allele by
    those bases, turned by the offset where the 3' rule moved the record (`lone_insertion_gaps`).
    A run under `--shift_3prime 1` prints a moved record's allele as the rule moved it, so on a
    row that carries the offset a minimised allele is the notation's bases as they stand or
    turned back by the offset, its span's gaps read as named or moved."""
    location, allele, feature, _feature_type, _consequences = key
    if name != "HGVSp" or ":p." not in vep or vep_rs not in ("", "-") or not intron_starts:
        return False
    perl_row, rust_row = perl_row or {}, rust_row or {}
    hgvsc = perl_row.get("HGVSc", "")
    if hgvsc in ("", "-") or rust_row.get("HGVSc") != hgvsc:
        return False
    offsets = (perl_row.get("HGVS_OFFSET"), rust_row.get("HGVS_OFFSET"))
    if None in offsets:
        return False
    moved = 0
    if not all(offset in ("", "-") for offset in offsets):
        if offsets[0] != offsets[1] or NON_ZERO_INTEGER.fullmatch(offsets[0]) is None:
            return False
        moved = int(offsets[0])
    transcript, _, body = hgvsc.partition(":")
    if transcript.split(".")[0] != feature.split(".")[0] or not insertion_after_the_exon_last_base(body):
        return False
    row = (allele.upper(), notation_inserted_bases(INTRON_START_INSERTION_BODY, body), 0, perl_row.get("STRAND", ""))
    return location_in_gaps(location, feature, intron_starts, moved, row)


def classify_field_divergence(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                              strand: str, flags: list[str],
                              edited_models: dict[str, frozenset[int]] | None = None,
                              codon_tables: dict[str, int] | None = None,
                              coding_sequences: dict[str, tuple[str, str, int]] | None = None,
                              perl_row: dict[str, str] | None = None,
                              rust_row: dict[str, str] | None = None,
                              minimal_split_record: bool = False,
                              coding_models: dict[tuple[str, str], CodingModel] | None = None,
                              genomic_span: tuple[int, int] | None = None,
                              partial_codon_models: frozenset[str] | None = None,
                              intron_ends: dict[tuple[str, str], frozenset[int]] | None = None,
                              intron_starts: dict[tuple[str, str], frozenset[int]] | None = None,
                              two_windows: dict[tuple[str, str], frozenset[int]] | None = None,
                              spliceai_table: object | None = None) -> tuple[str, str | None] | None:
    """The documented class of a field divergence on a key whose consequence sets agree, with
    the class of a second cause that changes the value (`None` but for a Downstream row two
    causes change, `downstream_field_class`), or
    `None`: `minimal_plugin_sees_the_rejoined_record` for a field of a covered plugin
    (`rejoined_record_plugin_field`) on a row of a multi-allelic record `--minimal` split
    (`minimal_split_record`, read from the run's flags, the row's `MINIMISED` and the record's
    ALT column), where
    Ensembl VEP's output pass hands the plugin the rejoined record and vep-rs the allele's own
    span and alleles, `motif_score_reverse_strand_mnv` for `MOTIF_SCORE_CHANGE` on a MotifFeature row of
    `STRAND` -1 whose allele is two or more bases over as many positions (Ensembl VEP writes
    the allele from the motif position of the variant's genomic start, the allele's last
    position along a reverse-strand motif, MotifFeatureVariationAllele.pm:245),
    `cell_type_comma_names` for a `CELL_TYPE` Ensembl VEP leaves empty under `--cell_type`
    where vep-rs prints an activity (its comma split never matches a name holding a comma,
    Config.pm:760-762 and OutputFactory.pm:1959-1969), `refseq_edited_model_hgvsc_offset`
    for an `HGVSc` on a transcript of `edited_models` (`indel_edited_models`) whose two
    notations differ in their cDNA coordinates alone and by one of the model's edit distances
    (some integer of the one notation lies at such a distance from some integer of the other):
    vep-rs numbers the position against the edited model, where Ensembl VEP walks the exons in
    their unedited numbering, its cache dumper having dropped the edit-aware mapper, and
    corrects an exonic `NM_` or `XM_` position alone by the misalignment offset
    (TranscriptVariationAllele.pm:1517, 1524-1525, 1608-1610), so an exon coordinate or an
    intronic distance of its notation sits a whole edit away from vep-rs's,
    `hgvs_added_bases_written_outside_the_transcript` for an `HGVSc`, `HGVS_OFFSET` or
    HGVSIntronOffset plugin field on the row of inserted bases both engines annotate upstream or
    downstream of the transcript, where Ensembl VEP prints no notation (OutputFactory.pm:1721)
    and vep-rs the `dup` or `ins` the HGVS 3' rule places on the transcript
    (`added_bases_written_outside_the_transcript`, which takes the row's other fields and the
    reference's `DISTANCE` from `perl_row` and `rust_row`),
    `refseq_edited_model_hgvs_offset` for an `HGVS_OFFSET` on a transcript of `edited_models`
    that Ensembl VEP leaves empty and vep-rs fills with a non-zero integer (on an edited `NM_`
    or `XM_` model `_return_3prime` copies the first shift hash attached to the variant,
    whatever transcript or strand made it, TranscriptVariationAllele.pm:164-187, so a row on
    the other strand inherits that transcript's shift, a zero that prints nothing; vep-rs
    shifts the row on its own transcript), and
    `organelle_peptide_paths_translation_table` for an `HGVSp` on a transcript `codon_tables`
    files under a code other than the standard one, both engines naming the same protein,
    where the row's consequence set and the two descriptions place the difference in one of
    the four paths of `organelle_peptide_path_divergence` (Ensembl VEP translates a sequence of
    each with the standard code, the alternate at TranscriptVariationAllele.pm:2333, 2492 and
    2555 and the reference at 2450, against a reference peptide of the transcript's code),
    `hgvs_shifted_indel_written_span_gate` for an `HGVSc`, `HGVS_OFFSET` or HGVSIntronOffset
    plugin field Ensembl VEP leaves empty on a shifted insertion or deletion vep-rs describes at
    its shifted span
    (`shifted_indel_written_span_gate`, reading the row's other field from `perl_row` and
    `rust_row`, the compared fields of the two rows), `hgvs_shifted_deletion_clamped_span`
    for an `HGVSc` whose two deletions end together and Ensembl VEP's starts later, on a row
    both engines give a non-zero `HGVS_OFFSET` (`shifted_deletion_clamped_span`),
    `hgvs_clipped_reference_against_whole_alternate` for an `HGVSc` on a sequence allele's row
    whose cDNA position has an open side, Ensembl VEP's a notation of the row's transcript read
    from the whole alternate allele over the bases inside the transcript, or none where that
    allele equals them, and vep-rs's the change inside, within the row's span of Ensembl VEP's,
    or none (`clipped_reference_against_whole_alternate`), and for the HGVSIntronOffset plugin's
    two fields that Ensembl VEP alone writes beside an `ins` or `dup` the class takes against a
    vep-rs row without an `HGVSc`, that vep-rs alone writes beside an `ins` inside the span at
    a reverse-strand transcript's first base where Ensembl VEP prints no `HGVSc`, or that both
    write beside two `ins` notations the class
    takes as a pair, each the intron part of its own side's notation's coordinate it is read
    from, and for the plugin's end offset that Ensembl VEP writes as -1 beside an `ins` from
    one base before a forward-strand transcript's first base where vep-rs writes 0 beside the
    duplication of the bases the allele adds (`clipped_notation_intron_offsets`, reading the
    row's other fields, each side's `BAM_EDIT` among them, from `perl_row` and
    `rust_row`), and
    `hgvsp_stop_codon_protein_comparison` for an `HGVSp` at the stop codon of one protein,
    Ensembl VEP's description read from the codon window and vep-rs's the one the two proteins
    give (`stop_codon_protein_comparison`, reading the reference row's `FLAGS` from
    `perl_row`), and `hgvsp_stop_named_on_incomplete_coding_sequence` for an `HGVSp` in which
    Ensembl VEP names the stop alone on a row that shows an incomplete coding sequence, or on
    a transcript of `partial_codon_models`, and vep-rs writes none
    (`stop_named_on_incomplete_coding_sequence`), and one of three classes for a field
    of the Downstream plugin on a frameshift row of a transcript of `coding_sequences`
    (`downstream_sequences`) whose sequences give exactly the two values under the transcript's
    codon table (`downstream_field_class`, reading the span from `perl_row`).
    `plugin_hgvs_shifted_cds_coordinates` is a `DownstreamProtein` or `ProteinLengthChange`
    Ensembl VEP computes, under `--hgvs`, on the CDS span and the allele its HGVS shift leaves
    (TranscriptVariationAllele.pm:844-858, :1330-1345, :2406-2410; Downstream.pm:108-110) where
    that changes its value and vep-rs computes on the record's own; the shifted span is the
    record's (`genomic_span`, the start and end of `row_record_span`) moved by the row's
    `HGVS_OFFSET` on the transcript's model in `coding_models` (`downstream_coding_models`).
    `downstream_start_past_own_codon` and `downstream_utr_frame_past_coding_end` are a
    `DownstreamProtein` on a row whose own span opens on a codon's third base: Ensembl VEP's
    plugin translates from the codon after the one the change opens in (Downstream.pm:112-115),
    and translates the 3' UTR from its first base when the change inserts no base and runs to
    the last coding base (Downstream.pm:115-119), where vep-rs translates from the changed
    codon. The class is
    `hgvsp_intron_end_insertion` for
    an `HGVSp` Ensembl VEP writes and vep-rs leaves empty on an insertion written between an
    intron's last base and the first coding base of the next exon of the row's transcript, a gap
    of `intron_ends` (`intron_end_gaps` of the corpus cache), the two rows carrying one `HGVSc`
    of one of the class's four forms, a duplication ending on the intron's last base or any
    insertion after it, and no `HGVS_OFFSET` (`intron_end_insertion`, reading both fields from
    `perl_row` and `rust_row`); a row with an `HGVS_OFFSET` and the end of an intron of 13 bases
    or fewer name nothing, and `hgvsp_intron_start_insertion` for an `HGVSp` Ensembl VEP writes
    and vep-rs leaves empty on an insertion the 3' rule leaves between a coding exon's last base
    and the first base of the intron after it, a gap of `intron_starts` (`intron_start_gaps` of
    the corpus cache), the two rows carrying one `HGVSc` that inserts bases between that coding
    base and the intron's first position, and no `HGVS_OFFSET` or one non-zero integer, the
    row's Location naming the gap or the offset carrying it there (`intron_start_insertion`,
    reading both fields from `perl_row` and `rust_row`); either intron class names as well a row
    whose record holds the insertion with a reference base it shares, the Allele that base and
    the inserted bases at a Location of the reference span, by the gap the allele clipped to the
    bases the `HGVSc` inserts lies in (`held_insertion_gaps`), the Allele read on either strand
    of the record, or, where the Allele is the insertion alone at such a Location (a minimised
    allele), by any gap of the span (`lone_insertion_gaps`). `maxentscan_overlapping_windows` for a
    `MaxEntScan_` field of a row whose one-base record lies in two canonical windows of its
    transcript in `two_windows` (`maxentscan_overlapping_windows`). `spliceai_multi_entry_cell` is a
    SpliceAI field on a row whose score record, among the files of `spliceai_table` (the format
    comparator's `SpliceAiTable`), holds several entries (`spliceai_multi_entry_cell`, reading
    the rows' SpliceAI fields from `perl_row` and `rust_row`); nothing takes it without the
    files."""
    location, allele, feature, feature_type, consequences = key
    if minimal_split_record and rejoined_record_plugin_field(name):
        return MINIMAL_PLUGIN_SEES_THE_REJOINED_RECORD, None
    if name == "MOTIF_SCORE_CHANGE" and feature_type == "MotifFeature" and strand == "-1":
        start, _, end = location.split(":", 1)[-1].partition("-")
        bases = int(end) - int(start) + 1 if start.isdigit() and end.isdigit() else 1
        if len(allele) >= 2 and bases == len(allele) and re.fullmatch(r"[ACGT]+", allele):
            return MOTIF_SCORE_REVERSE_STRAND_MNV, None
    if (name == "CELL_TYPE" and vep == "" and vep_rs not in ("", "-")
            and any(f == "--cell_type" or f.startswith("--cell_type=") for f in flags)):
        return CELL_TYPE_COMMA_NAMES, None
    edit_distances = (edited_models or {}).get(feature)
    if name == "HGVSc" and edit_distances and ":" in vep and ":" in vep_rs:
        vep_name, vep_body = vep.split(":", 1)
        vep_rs_name, vep_rs_body = vep_rs.split(":", 1)
        same_shape = HGVSC_COORDINATE.sub("#", vep_body) == HGVSC_COORDINATE.sub("#", vep_rs_body)
        vep_integers = [int(x) for x in HGVSC_INTEGER.findall(vep_body)]
        vep_rs_integers = [int(x) for x in HGVSC_INTEGER.findall(vep_rs_body)]
        an_edit_apart = any(abs(a - b) in edit_distances for a in vep_integers for b in vep_rs_integers)
        if vep_name == vep_rs_name and same_shape and an_edit_apart:
            return REFSEQ_EDITED_MODEL_HGVSC_OFFSET, None
    if added_bases_written_outside_the_transcript(name, key, vep, vep_rs, strand, perl_row, rust_row):
        return HGVS_ADDED_BASES_WRITTEN_OUTSIDE_THE_TRANSCRIPT, None
    if name == "HGVS_OFFSET" and edit_distances and vep in ("", "-") and NON_ZERO_INTEGER.fullmatch(vep_rs):
        return REFSEQ_EDITED_MODEL_HGVS_OFFSET, None
    if name == "HGVSp" and feature in (codon_tables or {}):
        vep_name, _, vep_description = vep.rpartition(":p.")
        vep_rs_name, _, vep_rs_description = vep_rs.rpartition(":p.")
        if (vep_name and vep_name == vep_rs_name
                and organelle_peptide_path_divergence(consequences, vep_description, vep_rs_description)):
            return ORGANELLE_PEPTIDE_PATHS_TRANSLATION_TABLE, None
    if shifted_indel_written_span_gate(name, key, vep, vep_rs, perl_row, rust_row):
        return HGVS_SHIFTED_INDEL_WRITTEN_SPAN_GATE, None
    if name == "HGVSc" and shifted_deletion_clamped_span(key, vep, vep_rs, perl_row):
        return HGVS_SHIFTED_DELETION_CLAMPED_SPAN, None
    if name == "HGVSc" and clipped_reference_against_whole_alternate(key, vep, vep_rs, perl_row):
        return HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE, None
    if clipped_notation_intron_offsets(name, key, vep, vep_rs, perl_row, rust_row):
        return HGVS_CLIPPED_REFERENCE_AGAINST_WHOLE_ALTERNATE, None
    if name == "HGVSp" and stop_codon_protein_comparison(key, vep, vep_rs, perl_row):
        return HGVSP_STOP_CODON_PROTEIN_COMPARISON, None
    if name == "HGVSp" and stop_named_on_incomplete_coding_sequence(key, vep, vep_rs, perl_row,
                                                                       partial_codon_models):
        return HGVSP_STOP_NAMED_ON_INCOMPLETE_CODING_SEQUENCE, None
    if name in DOWNSTREAM_FIELDS:
        return downstream_field_class(name, key, vep, vep_rs, strand, perl_row, rust_row, coding_sequences,
                                      coding_models, genomic_span)
    if intron_end_insertion(name, key, vep, vep_rs, perl_row, rust_row, intron_ends):
        return HGVSP_INTRON_END_INSERTION, None
    if spliceai_multi_entry_cell(name, key, vep, vep_rs, perl_row, rust_row, spliceai_table):
        return SPLICEAI_MULTI_ENTRY_CELL, None
    if intron_start_insertion(name, key, vep, vep_rs, perl_row, rust_row, intron_starts):
        return HGVSP_INTRON_START_INSERTION, None
    if maxentscan_overlapping_windows(name, key, vep, vep_rs, two_windows):
        return MAXENTSCAN_OVERLAPPING_WINDOWS, None
    return None


def lists_a_zero_frequency(freqs: str) -> bool:
    """Whether a FREQS value (`POPULATION:allele:frequency` items joined by `,`) gives any
    allele a frequency of exactly 0; `NA`, an absent column (`-`) and a text that is no
    number are not zero."""
    for item in freqs.split(","):
        value = item.rsplit(":", 1)[-1]
        try:
            if float(value) == 0.0:
                return True
        except ValueError:
            continue
    return False


def zero_frequency_keys(path: Path) -> set[tuple[str, str, str, str]]:
    """The (Location, Allele, Feature, Feature_type) keys of a vep-rs default-format output
    whose rows list a frequency of exactly 0 in FREQS (`lists_a_zero_frequency`)."""
    return {key[:4] for key, fields in load_fields(path, ["FREQS"]).items() if lists_a_zero_frequency(fields["FREQS"])}


# A VCF REF column: bases or `.`, never the `/`-joined allele string or the class word of an
# Ensembl-format line, so a line of five or more columns with one is a VCF record.
VCF_REF = re.compile(r"[ACGTNacgtn.]+")


def first_five_columns(line: str) -> str:
    """A VCF line's CHROM, POS, ID, REF and ALT, the columns a twin of the line (another
    encoding, sample columns added or dropped) shares with it."""
    return "\t".join(line.split("\t")[:5])


def carried_records(records: list[dict], json_path: Path) -> set[int]:
    """Ordinals of the records the reference carried with no row: the JSON output at
    `json_path`, when the file exists, has an object for the record whose
    `most_severe_consequence` is `?`, the object VEP writes for a feature with no overlap
    allele (every ALT a `*` or `<DEL:*>`, equal to the reference or a repeat of an earlier
    one; VariationFeatureOverlap.pm:459-491) and for a record a row filter (`--coding_only`,
    `--no_intergenic`) left no consequence of (JSON.pm:387-388). A record is the object whose
    `input` opens with its first five columns, so a twin input carrying sample columns names
    it too. A file that is not gzip, or holds a line that is not a JSON object, is refused
    with its path."""
    if not json_path.is_file():
        return set()
    by_head: dict[str, list[int]] = defaultdict(list)
    for i, rec in enumerate(records):
        by_head[first_five_columns(rec["vcf"])].append(i)
    carried: set[int] = set()
    try:
        with gzip.open(json_path, "rt", encoding="utf-8", errors="replace") as fh:
            for number, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                obj = json.loads(line)
                if not isinstance(obj, dict):
                    raise failure(f"{json_path} is not a readable JSON output (line {number} is not a JSON object)")
                if obj.get("most_severe_consequence") != "?":
                    continue
                head = first_five_columns(str(obj.get("input", "")))
                carried.update(by_head.get(head, []))
    except (OSError, EOFError, zlib.error, json.JSONDecodeError) as exc:
        raise failure(f"{json_path} is not a readable JSON output ({exc})") from None
    return carried


def run_records(corpus: Path, manifest: dict, target: dict) -> list[dict]:
    """The records of a run, numbered as the golden harness numbers them (every line of the
    run's input that is not a header, in order): the manifest's `records` for a run over
    `variants.vcf`; for a flag set over another `input`, that file's lines. A line sharing
    its first five columns with a manifest record is that record (a twin of `variants.vcf`
    in another encoding keeps its `reference_warning` and takes `reference_rows`). A file
    holding a VCF record (`vcf_record`: five tab-separated columns or more, a POS of digits,
    a REF of bases) is a VCF whose other lines are records of their own, which carry no
    `reference_warning`, so a record the reference skipped in such an input is refused as
    unexplained. A file holding none is `variants.vcf` re-encoded in another format or
    delimiter (`--format hgvs`, an Ensembl-format line, a space-delimited VCF line), record
    for record, and each line is the manifest record at its ordinal; such a file whose
    record count is not the manifest's is refused with its path, as is a line of a VCF input
    that is neither a twin nor a VCF record, and a file that cannot be read."""
    input_name = target.get("input") or DEFAULT_INPUT
    if input_name == DEFAULT_INPUT:
        return manifest["records"]
    by_head = {first_five_columns(rec["vcf"]): rec for rec in manifest["records"]}
    path = corpus / input_name
    try:
        with opener(path)(path, "rt", encoding="utf-8", errors="replace") as handle:
            lines = [line.rstrip("\n") for line in handle if not line.startswith("#")]
    except (OSError, EOFError, zlib.error) as exc:
        raise failure(f"{path} is not a readable text file ({exc})") from None

    def vcf_record(line: str) -> bool:
        columns = line.split("\t")
        return len(columns) >= 5 and columns[1].isdecimal() and VCF_REF.fullmatch(columns[3]) is not None

    if any(vcf_record(line) for line in lines):
        records: list[dict] = []
        for i, line in enumerate(lines):
            manifest_record = by_head.get(first_five_columns(line))
            if manifest_record is not None:
                records.append(manifest_record)
            elif vcf_record(line):
                records.append({"vcf": line})
            else:
                raise failure(f"{path} record {i} ({line[:40]!r}) is no VCF record in a VCF input")
        return records
    if len(lines) != len(manifest["records"]):
        raise failure(f"{path} holds no VCF record and its {len(lines)} records are not the manifest's "
                      f"{len(manifest['records'])} re-encoded")
    return manifest["records"]


def reference_skipped_records(
    records: list[dict],
    perl_names: dict[tuple[str, str, str, str], set[str]],
    by_name: dict[str, list[int]],
    carried: set[int] = frozenset(),
) -> tuple[set[int], list[str]]:
    """Ordinals of the records no VEP row names (VEP dropped them before annotation), and the
    records that contradict their `reference_warning`: a zero-row record without one (its name
    unresolved, or the reason unrecorded) and a warned record with rows. A record that carries
    no row by its nature is kept, not skipped, when it has no warning: a non-variant record (ALT
    `.`, written under `--allow_non_variant` as its VCF line and JSON object and no default row,
    Parser/VCF.pm:259-267) and a record of `carried`, one the reference's JSON names as a
    feature with no allele to annotate (`carried_records`), one a row filter or the frequency
    filter of the run left without a row on either side, or one the frequency filter dropped
    and vep-rs kept for a zero frequency (`cmd_classify`)."""
    with_rows = {i for names in perl_names.values() for n in names for i in by_name.get(n, [])}
    skipped = set(range(len(records))) - with_rows
    problems: list[str] = []
    for i, rec in enumerate(records):
        warned = isinstance(rec.get("reference_warning"), str)
        columns = rec["vcf"].split("\t")
        label = " ".join(columns[:5])
        non_variant = len(columns) > 4 and columns[4] == "."
        if i in skipped and not warned:
            if non_variant or i in carried:
                skipped.discard(i)
            else:
                problems.append(f"record {i} ({label}) has no VEP row and no reference_warning")
        elif warned and i not in skipped:
            problems.append(f"record {i} ({label}) has VEP rows and reference_warning {rec['reference_warning']!r}")
    return skipped, problems


def format_comparator():
    """The format comparator, `scripts/concordance/compare_vep_formats.py`, as a module."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
    import compare_vep_formats

    return compare_vep_formats


def spliceai_multi_entry_cell(name: str, key: tuple[str, str, str, str, str], vep: str, vep_rs: str,
                              perl_row: dict[str, str] | None, rust_row: dict[str, str] | None,
                              table: object | None) -> bool:
    """Whether the two values of the SpliceAI field `name` on a key are the difference a score
    record of several entries makes, as the format comparator's class of the same name decides
    it: from the two rows' SpliceAI fields and `SYMBOL` (`SPLICEAI_ROW_FIELDS`), the corpus
    cache's symbol of the row's transcript where neither row writes a `SYMBOL`, and the records
    of `table` at the row's Location and Allele. Nothing without the table."""
    if table is None:
        return False
    formats = format_comparator()

    def entry(row: dict[str, str] | None):
        fields = {field: value for field, value in (row or {}).items() if field in SPLICEAI_ROW_FIELDS and value != "-"}
        return formats.Entry(key="", record=None, location=key[0], allele=key[1], feature=key[2], consequence=key[4], fields=fields)

    return formats.field_class_accepts(SPLICEAI_MULTI_ENTRY_CELL, entry(perl_row), name, vep, vep_rs, entry(rust_row), spliceai_table=table)


def spliceai_members(row: dict[str, str]) -> list[str] | None:
    """The nine members of a default-format row's SpliceAI value: the first nine `|`-joined
    parts of its `SpliceAI_pred`, else its nine `SpliceAI_pred_` fields; `None` for a row that
    writes neither whole."""
    whole = row.get("SpliceAI_pred", "-")
    members = whole.split("|")[:9] if whole != "-" else [row.get(field, "-") for field in SPLICEAI_MEMBER_FIELDS]
    return members if len(members) == 9 and "-" not in members else None


def spliceai_json_object(members: list[str] | None, flag: str) -> str:
    """A row's JSON `spliceai` object as the golden harness prints an object: no spaces, its
    members in the byte order of their names, the symbol a string, any other member a number
    where its text is one and a string otherwise, and the cutoff flag among them where the row
    carries one; `-` for a row without the value. Both engines write a number by its value,
    so a score written `0.00` is the integer `0` there and `0.10` is `0.1`."""
    if members is None:
        return "-"

    def number(text: str) -> int | float:
        value = float(text)
        return int(value) if value.is_integer() else value

    names = [field.removeprefix("SpliceAI_pred_") for field in SPLICEAI_MEMBER_FIELDS]
    written: dict[str, object] = {
        name: number(member) if name != "SYMBOL" and SPLICEAI_NUMBER.fullmatch(member) else member
        for name, member in zip(names, members)
    }
    if flag != "-":
        written["SpliceAI_cutoff"] = flag
    return json.dumps(written, sort_keys=True, separators=(",", ":"))


def spliceai_other_format_entries(written: dict, perl_row: dict[str, str], rust_row: dict[str, str],
                                  compared: Collection[str]) -> list[dict]:
    """The entries the VCF and JSON formats need for a row `classify` writes under
    `spliceai_multi_entry_cell` (`written`, one of its entries), each value derived from the
    two rows' default-format values: one per member field of `SPLICEAI_MEMBER_FIELDS` whose two
    members differ, left out where `compared` names the field and the run's default format
    writes it itself, and one for the JSON object (`spliceai_json_object`). The golden harness
    reads a documented field by the name its format gives it."""
    theirs, ours = spliceai_members(perl_row), spliceai_members(rust_row)
    absent = ["-"] * len(SPLICEAI_MEMBER_FIELDS)
    entries = [
        {**written, "field": field, "vep": vep, "vep_rs": vep_rs}
        for field, vep, vep_rs in zip(SPLICEAI_MEMBER_FIELDS, theirs or absent, ours or absent)
        if vep != vep_rs and field not in compared
    ]
    objects = (spliceai_json_object(members, row.get("SpliceAI_cutoff", "-")) for members, row in ((theirs, perl_row), (ours, rust_row)))
    entries.append({**written, "field": SPLICEAI_JSON_FIELD, **dict(zip(("vep", "vep_rs"), objects))})
    return entries


def classify_pair(perl_csq: str, rust_csq: str, allele: str) -> str | None:
    """First matching rule of the key comparator's registry, by taxonomy class."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
    import compare_vep_outputs as cvo  # noqa: E402

    for rule in cvo.EXCLUSION_REGISTRY:
        if rule.matches(perl_csq, rust_csq, allele):
            return rule.taxonomy_class
    return None


def classify_mate_side_pair(
    location: str, allele: str, feature: str, perl_csq: str, rust_csq: str, models: dict
) -> str | None:
    """`breakend_mate_local_read` when VEP's set is the local-coordinate read of a point
    breakend's mate-side Transcript row and vep-rs's set is the mate-coordinate read, the
    two differing; the derivation is the SV comparator's `breakend_mate_context`, over the
    corpus's own cache."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc  # noqa: E402

    mate = bmc.bracket_mate(allele)
    if mate is None:
        return None
    model = models.get((mate[0], feature))
    if model is None:
        return None
    derived = bmc.mate_side_derivation(location, allele, model)
    if derived is None:
        return None
    mate_set, local_set = derived
    if mate_set == local_set:
        return None
    if set(perl_csq.split(",")) == local_set and set(rust_csq.split(",")) == mate_set:
        return "breakend_mate_local_read"
    return None


# The five classes of the stop codon tests behind `stop_retained` (Utils/VariationEffect.pm
# :1290-1324) that `classify_stop_codon_pair` names. `stop_retained_without_stop_codon_test`:
# `_ins_del_stop_altered_cil` (:1488-1539) returns 0 before it edits anything for an allele
# holding a letter that names no base (:1498) and for a span without a cDNA coordinate at one
# end or a CDS start (:1508), and :1318 reads that 0 as a stop codon kept; vep-rs reads a span
# the test cannot be made on as keeping nothing. `stop_codon_overlap_genomic_span`:
# `_overlaps_stop_codon_cil` (:1418-1422) takes the last three genomic positions of the coding
# region as the stop codon, which on a codon an intron splits miss the stop bases that end the
# exon before; vep-rs reads the codon's own bases and makes the test.
# `stop_retained_three_base_stop_codon_test`: where the test is made it reads the three bases
# at the old stop position alone (:1528-1535), so a stop codon there is a stop codon kept
# whatever the edit did to the codons before it; vep-rs keeps the term only where each of those
# codons reads the reference's residue as well, and `frameshift` (:1541-1567) then reads away
# an in-frame or protein-altering term. `stop_codon_removed_without_stop_codon_test`: a span
# that holds the whole stop codon and has no coordinates goes through the same returns, and
# through their twin behind `stop_lost` (:1376, :1453); vep-rs writes `stop_lost` unless the
# sequence the record gives shows the stop codon kept. `stop_codon_overlap_insertion_range`:
# `_overlaps_stop_codon_cil` moves an insertion's start by the inserted length (:1408-1415),
# so on the forward strand an insertion before or inside the stop codon that keeps the codon
# and every base before it misses the overlap; vep-rs reads the insertion on the coding
# sequence and writes `stop_retained_variant`.
# `classify_stop_codon_pair` names one more class of the key comparator's row-aware rules.
# `splice_term_past_the_record_span`: On a record whose alternate allele is longer than its
# reference allele, Ensembl VEP writes a splice donor, splice acceptor, donor fifth base, donor
# region, polypyrimidine tract or intron term from positions past the end of the record's
# span; where no base of the span lies in a window of that term, vep-rs does not write it.
# `_get_differing_regions` (VariationFeatureOverlapAllele.pm:405-460) compares the two alleles
# to the longer one's length, and `_intron_effects` places each differing region from the
# record's start with no bound at the span's end (BaseTranscriptVariationAllele.pm:117-118),
# where it sets the flags of the six terms (:143-204); the intron is the intron term's window.
# A row is read with its record's alleles where the builder holds them, in the reference row's
# first column or among the run's records by position (`record_reference`): the span may then
# lie inside an exon, across its end or in an intron, and may be two positions. A row no
# record of the run writes is read without them, on three or more positions that lie in the
# transcript's sequence.
STOP_RETAINED_WITHOUT_STOP_CODON_TEST = "stop_retained_without_stop_codon_test"
STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST = "stop_codon_removed_without_stop_codon_test"
STOP_LOST_TERM = "stop_lost"
STOP_CODON_OVERLAP_GENOMIC_SPAN = "stop_codon_overlap_genomic_span"
STOP_RETAINED_THREE_BASE_STOP_CODON_TEST = "stop_retained_three_base_stop_codon_test"
STOP_CODON_OVERLAP_INSERTION_RANGE = "stop_codon_overlap_insertion_range"
SPLICE_TERM_PAST_THE_RECORD_SPAN = "splice_term_past_the_record_span"
STOP_RETAINED_TERM = "stop_retained_variant"
CODING_SEQUENCE_TERM = "coding_sequence_variant"
FRAMESHIFT_TERM = "frameshift_variant"
UNTESTED_STOP_CODON_REPLACEMENTS = (frozenset(), frozenset({CODING_SEQUENCE_TERM}), frozenset({FRAMESHIFT_TERM}))
FRAMESHIFT_REPLACEMENTS = (frozenset({FRAMESHIFT_TERM}), frozenset({"stop_gained", FRAMESHIFT_TERM}))
RETAINED_INSERTION_TERMS = (frozenset({STOP_RETAINED_TERM}), frozenset({"inframe_insertion", STOP_RETAINED_TERM}))


def load_positions(path: Path, masked: frozenset[str] = frozenset()) -> dict[tuple[str, str, str, str, str], list[tuple[str, str, str]]]:
    """(Location, Allele, Feature, Feature_type, consequence set) -> every distinct row that
    carries the tuple, in the file's order, as its `cDNA_position` and `CDS_position` columns
    and its `Uploaded_variation`, leaving out the rows whose Feature_type is in `masked`."""
    rows: dict[tuple[str, str, str, str, str], list[tuple[str, str, str]]] = defaultdict(list)
    for r in iter_rows(path, "x"):
        if r.feature_type in masked:
            continue
        row = (r.positions[0], r.positions[1], r.uploaded_variation)
        carried = rows[(r.location, r.allele, r.feature, r.feature_type, r.consequence_set)]
        if row not in carried:
            carried.append(row)
    return dict(rows)


def record_reference(records: Iterable[dict], location: str, allele: str) -> str | None:
    """The reference bases a row's `allele` stands in for at its `location`, read by position
    from the alleles of a run's `records` (`run_records`) as the key comparator reads them from
    an input VCF (`read_record_references` of `compare_vep_outputs.py`): a record on the
    Location's chromosome, a `chr` prefix set aside, that writes the row gives them
    (`replaced_reference`), the empty string for a pure insertion's row. `None` where no record
    writes the row and where two records read it differently."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
    import compare_vep_outputs as cvo

    def unprefixed(chrom: str) -> str:
        return chrom[3:] if chrom.lower().startswith("chr") else chrom

    chrom, sep, _ = location.partition(":")
    if not sep:
        return None
    found: set[str] = set()
    for rec in records:
        fields = rec["vcf"].split("\t", 5)
        if len(fields) < 5 or unprefixed(fields[0]) != unprefixed(chrom):
            continue
        position = cvo.vep_rs_pos(fields[1])
        if position is None:
            continue
        bases = cvo.replaced_reference(position, fields[3], cvo.perl_split(fields[4], ","), location, allele)
        if bases is not None:
            found.add(bases)
    return found.pop() if len(found) == 1 else None


def comparator_stop_codon_class(
    location: str, allele: str, perl_csq: str, rust_csq: str, positions: tuple[str, str], candidate_positions: tuple[str, str],
    uploaded_variation: str, records: Collection[dict] = (),
) -> str | None:
    """The class whose row-aware predicate in the key comparator's 116.2 registry takes the
    pair of Transcript rows, asked in the registry's order: the reference's row with its first
    column and its `cDNA_position` and `CDS_position`, vep-rs's with its two positions. The
    reference's row carries the reference bases its allele stands in for where its first column
    spells the alleles of a record that writes the row; else both rows carry them as the run's
    `records` give them by position (`record_reference`), and neither where no record writes
    the row. The five stop codon classes and `splice_term_past_the_record_span` are defined
    once, in `compare_vep_outputs.py`."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "concordance"))
    import compare_vep_outputs as cvo

    reference = cvo.pair_row([uploaded_variation, location, allele, "", "", "Transcript", perl_csq, *positions])
    candidate = cvo.pair_row(["", location, allele, "", "", "Transcript", rust_csq, *candidate_positions])
    if reference.reference is None:
        reference.reference = candidate.reference = record_reference(records, location, allele)
    for rule in cvo.REFERENCE_116_2_ROW_PAIR_RULES:
        if rule.predicate(reference, candidate):
            return rule.taxonomy_class
    return None


def classify_stop_codon_pair(
    location: str, allele: str, feature_type: str, perl_csq: str, rust_csq: str, positions: tuple[str, str] | None,
    untested_records: Collection[tuple[str, str]] = (),
    retained_insertion_locations: Collection[str] = (),
    uploaded_variation: str = "",
    candidate_positions: tuple[str, str] | None = None,
    records: Collection[dict] = (),
) -> str | None:
    """The class of one reference row and one vep-rs row of a tuple whose two consequence sets
    differ by a stop codon term or by a splice or intron term the reference alone holds,
    `positions` being the reference row's `cDNA_position` and `CDS_position`,
    `uploaded_variation` its first column, `candidate_positions` the two positions of vep-rs's
    row, for which the reference's stand in where a caller has none, and `records` the run's
    records (`run_records`), whose alleles a row under an identifier is read with.
    A Transcript row takes the class whose predicate in the key comparator takes the pair of
    rows (`comparator_stop_codon_class`), on an insertion's or deletion's row as that
    comparator reads one. `stop_retained_without_stop_codon_test`: the reference's
    `stop_retained_variant` against the same set without it, with nothing,
    `coding_sequence_variant` or `frameshift_variant` in its place, on rows that show the test
    is not made. `stop_codon_overlap_genomic_span`: the reference's `coding_sequence_variant`
    against vep-rs's `stop_retained_variant`, both beside `splice_region_variant`, on rows
    whose test is made. `stop_codon_removed_without_stop_codon_test`: the reference's
    `stop_retained_variant` against vep-rs's `stop_lost` on rows without a cDNA position at
    one end or a CDS start. `stop_retained_three_base_stop_codon_test`: on rows whose test is
    made, the reference's `stop_retained_variant` against `coding_sequence_variant`, or it,
    alone or beside one in-frame or protein-altering term, against `frameshift_variant` with
    or without `stop_gained`. `stop_codon_overlap_insertion_range`: on the row of an insertion
    between two neighbouring positions with both cDNA positions and a CDS start on both
    sides, the reference's `coding_sequence_variant` away from `splice_region_variant`, or
    its `frameshift_variant` with or without `stop_gained`, against vep-rs's
    `stop_retained_variant`, the second with or without `inframe_insertion` beside it.
    `splice_term_past_the_record_span`: on rows with the same cDNA and CDS positions, the
    reference's `splice_donor_variant`, `splice_acceptor_variant`,
    `splice_donor_5th_base_variant`, `splice_donor_region_variant`,
    `splice_polypyrimidine_tract_variant` and `intron_variant`, one or more, against a set of
    vep-rs's without them. With the record's alleles, from the first column where it spells
    them and else from `records` by position (`record_reference`), the record is no pure
    insertion, its reference bases are two or more and as many as the Location's positions,
    the allele is longer, and one intron stands beside the span, at a place the cDNA positions
    allow, at which the record's differing regions give the reference's terms and the span
    gives vep-rs's, `splice_donor_region_variant` and `splice_region_variant` coming where the
    terms that withhold them went: the span may lie inside an exon, across its end or in an
    intron. Without them the Location holds three or more positions under a longer allele,
    both rows carry two cDNA positions as far apart as the Location's ends, vep-rs's set
    holds none of the six, and `splice_region_variant` alone may come, where a donor,
    acceptor, donor-region or fifth-base term went. With or without them the two sets may
    differ by the terms of a stop codon class as well, where that class's rule takes the rest
    of the difference. A row that names no feature, the one
    row a record has under `--summary` or `--most_severe`, carries the terms of all its features
    and no position: it takes `stop_retained_without_stop_codon_test` or
    `stop_codon_removed_without_stop_codon_test` when its sets differ in that
    class's way and its Location and Allele are in `untested_records`, the records whose
    Transcript rows the manifest's own run documents under one of the two; and
    `stop_codon_overlap_insertion_range` when its sets differ in that class's way and its
    Location is in `retained_insertion_locations`, the records documented there under it. It
    takes `splice_term_past_the_record_span` nowhere and stays unexplained: the class reads
    one transcript's cDNA and CDS positions, which such a row does not carry."""
    if feature_type == "-":
        perl_only = set(perl_csq.split(",")) - set(rust_csq.split(","))
        rust_only = set(rust_csq.split(",")) - set(perl_csq.split(","))
        if (location in retained_insertion_locations and rust_only in RETAINED_INSERTION_TERMS
                and (perl_only in FRAMESHIFT_REPLACEMENTS or (perl_only, rust_only) == ({CODING_SEQUENCE_TERM}, {STOP_RETAINED_TERM}))):
            return STOP_CODON_OVERLAP_INSERTION_RANGE
        if (location, allele) not in untested_records or perl_only != {STOP_RETAINED_TERM}:
            return None
        if rust_only in UNTESTED_STOP_CODON_REPLACEMENTS:
            return STOP_RETAINED_WITHOUT_STOP_CODON_TEST
        if rust_only == {STOP_LOST_TERM}:
            return STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST
        return None
    if feature_type != "Transcript" or positions is None:
        return None
    return comparator_stop_codon_class(
        location, allele, perl_csq, rust_csq, positions, candidate_positions or positions, uploaded_variation, records)


def classify_stop_codon_rows(
    location: str, allele: str, feature_type: str, perl_csq: str, rust_csq: str,
    reference_rows: Collection[tuple[str, str, str]], candidate_rows: Collection[tuple[str, str, str]],
    untested_records: Collection[tuple[str, str]] = (),
    retained_insertion_locations: Collection[str] = (),
    records: Collection[dict] = (),
) -> str | None:
    """The class `classify_stop_codon_pair` names for every reference row that carries a tuple
    against every vep-rs row that carries the other (`load_positions`), as the key comparator's
    rule takes a pair only where every such pair of rows meets it (`ExclusionRule.matches`):
    `None` where a side has no row or two pairs of rows are named differently. A row that
    names no feature carries no position and is read by its terms alone. `records` are the
    run's records, whose alleles a row is read with where its first column spells none
    (`record_reference`)."""
    if feature_type == "-":
        return classify_stop_codon_pair(location, allele, feature_type, perl_csq, rust_csq, None, untested_records,
                                        retained_insertion_locations)
    if not reference_rows or not candidate_rows:
        return None
    classes = {
        classify_stop_codon_pair(location, allele, feature_type, perl_csq, rust_csq, reference[:2], untested_records,
                                 retained_insertion_locations, reference[2], candidate[:2], records)
        for reference in reference_rows for candidate in candidate_rows
    }
    return classes.pop() if len(classes) == 1 else None


def corpus_cache_dir(corpus: Path, manifest: dict) -> Path:
    """The corpus's cache directory: `json_cache/`, or under the Ensembl layout
    `cache_root/<species>/<release>_<assembly>` (`homo_sapiens` when the manifest names no
    species; the species suffixed `_refseq` and then `_merged` for each of `--refseq` and
    `--merged` the manifest's flags carry, as Ensembl VEP names the directory, CacheDir.pm:312
    at release/116.2), the directory `crates/vep-cli/tests/common/mod.rs` resolves."""
    if manifest.get("cache_layout") == "ensembl":
        species = manifest.get("species", DEFAULT_SPECIES)
        flags = manifest.get("flags") or []
        if "--refseq" in flags:
            species += "_refseq"
        if "--merged" in flags:
            species += "_merged"
        return corpus / "cache_root" / species / f"{manifest['release']}_{manifest['assembly']}"
    return corpus / "json_cache"


def load_corpus_models(corpus: Path) -> dict:
    """The corpus's pruned cache, located by the manifest's `cache_layout`, as
    `breakend_mate_context` transcript models."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
    import breakend_mate_context as bmc  # noqa: E402

    return bmc.load_transcript_models(corpus_cache_dir(corpus, read_manifest(corpus)))


def cmd_classify(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    target = flag_set_named(manifest, args.flag_set) if args.flag_set is not None else manifest
    if args.flag_set is not None and run_options(args):
        print("ERROR: [build_golden_corpus] --flag and --fasta-name describe the default run; "
              "a flag set's flags are its own entry's", file=sys.stderr)
        return 1
    masks = run_masks(manifest, target)
    masked = masked_feature_types(masks)
    perl, perl_names = load_tuples(Path(args.vep_default), masked)
    perl_positions = load_positions(Path(args.vep_default), masked)
    rust_positions = load_positions(Path(args.vep_rs_output))
    rust, rust_names = load_tuples(Path(args.vep_rs_output))
    records = run_records(corpus, manifest, target)
    by_name = record_name_index(records)
    # The default run's JSON stands in for a flag set whose expected directory has none:
    # `carried_records` names a record by its first five columns, which a twin input shares.
    reference_json = corpus / "expected" / "json.jsonl.gz"
    if args.flag_set is not None:
        expected_dir = target.get("expected_dir", f"expected_{args.flag_set}")
        own_json = corpus / expected_dir / "json.jsonl.gz"
        if own_json.is_file():
            reference_json = own_json
    rowless = carried_records(records, reference_json)
    if args.flag_set is not None:
        flags = target.get("flags") or []
    else:
        flags = run_options(args).get("flags") or manifest.get("flags") or []
    rust_keys_of: dict[int, set[tuple[str, str, str, str]]] = defaultdict(set)
    for key, names in rust_names.items():
        for name in names:
            for i in by_name.get(name, []):
                rust_keys_of[i].add(key)
    frequency_filtered = bool(FREQUENCY_FILTER_FLAGS.intersection(flags))
    if ROW_FILTERS.intersection(flags) or frequency_filtered:
        # A row filter leaves a record with no coding (or no intergenic) row rowless on both
        # sides, and the frequency filter drops a record whole on both; the vep-rs side is
        # read here, the reference side by the skip rule.
        rowless |= set(range(len(records))) - set(rust_keys_of)
    # The records Ensembl VEP dropped by frequency and vep-rs kept for a zero frequency: no
    # reference row, every vep-rs key listing the zero. They are rowless on the reference
    # side by the run's nature, and their vep-rs-only keys take the class.
    zero_records: set[int] = set()
    if frequency_filtered:
        zero_keys = zero_frequency_keys(Path(args.vep_rs_output))
        with_perl_rows = {i for names in perl_names.values() for n in names for i in by_name.get(n, [])}
        zero_records = {i for i, keys in rust_keys_of.items() if i not in with_perl_rows and keys <= zero_keys}
        rowless |= zero_records
    skipped, problems = reference_skipped_records(records, perl_names, by_name, rowless)
    if problems:
        for p in problems:
            print(f"ERROR: [build_golden_corpus] {p}", file=sys.stderr)
        return 1
    for i, rec in enumerate(records):
        rec.pop("reference_rows", None)
        if i in skipped:
            # The two markers are the record's last two keys, in this order.
            rec["reference_rows"] = 0
            rec["reference_warning"] = rec.pop("reference_warning")
    models = load_corpus_models(corpus)

    def records_for(key: tuple[str, str, str, str]) -> list[int]:
        found: set[int] = set()
        for n in perl_names.get(key, set()) | rust_names.get(key, set()):
            found.update(by_name.get(n, []))
        return sorted(found)

    # The records the manifest's own run documents under an untested stop codon class, which
    # a flag set's record-level row (`--summary`, `--most_severe`) is read against.
    untested_records = set() if target is manifest else {
        (d["location"], d["allele"]) for d in manifest.get("divergences") or []
        if d.get("expected_divergence") in (STOP_RETAINED_WITHOUT_STOP_CODON_TEST, STOP_CODON_REMOVED_WITHOUT_STOP_CODON_TEST)
    }
    retained_insertion_locations = set() if target is manifest else {
        d["location"] for d in manifest.get("divergences") or []
        if d.get("expected_divergence") == STOP_CODON_OVERLAP_INSERTION_RANGE
    }
    divergences: list[dict] = []
    for key, perl_sets in sorted(perl.items()):
        rust_sets = rust.get(key, set())
        for cs in sorted(perl_sets - rust_sets):
            klass = None
            for rs in sorted(rust_sets):
                klass = (
                    classify_pair(cs, rs, key[1])
                    or classify_mate_side_pair(key[0], key[1], key[2], cs, rs, models)
                    or classify_stop_codon_rows(key[0], key[1], key[3], cs, rs, perl_positions.get((*key, cs), ()),
                                                rust_positions.get((*key, rs), ()), untested_records, retained_insertion_locations,
                                                records)
                )
                if klass:
                    break
            divergences.append({
                "location": key[0], "allele": key[1], "feature": key[2], "feature_type": key[3],
                "record_indices": records_for(key),
                "vep_consequence_set": cs, "vep_rs_consequence_sets": sorted(rust_sets),
                "expected_divergence": klass or UNEXPLAINED,
            })
    extra_rust: list[dict] = []
    for k, v in sorted(rust.items()):
        if k in perl:
            continue
        indices = records_for(k)
        if indices and skipped.issuperset(indices):
            klass = REFERENCE_SKIPPED
        elif indices and zero_records.issuperset(indices):
            klass = FREQUENCY_FILTER_ZERO_VALUE
        else:
            klass = UNEXPLAINED
        extra_rust.append({
            "location": k[0], "allele": k[1], "feature": k[2], "feature_type": k[3],
            "record_indices": indices, "vep_rs_consequence_sets": sorted(v),
            "expected_divergence": klass,
        })
    # A with-reference key the target run holds is carried onto the entry written for the same
    # Location, Allele, Feature, records and Ensembl VEP set; one with no such entry, or of
    # another shape beside the entry's regenerated sets, stops the run before anything is written.
    def carried_key(d: dict) -> tuple:
        return (d.get("location"), d.get("allele"), d.get("feature"), tuple(d.get("record_indices") or ()), d.get("vep_consequence_set"))

    stated = {carried_key(d): d for d in target.get("divergences") or [] if WITH_REFERENCE_SETS in d}
    for d in divergences:
        held = stated.pop(carried_key(d), None)
        if held is not None:
            d[WITH_REFERENCE_SETS] = held[WITH_REFERENCE_SETS]
    problems = [
        f"the entry on {entry_key_text(d)} carries {WITH_REFERENCE_SETS}, and this run writes no entry of its "
        "Location, Allele, Feature, records and Ensembl VEP set to carry it on: remove the key, or classify the "
        "outputs of a run without a reference sequence"
        for d in stated.values()
    ]
    problems += [
        f"the entry on {entry_key_text(d)}: {WITH_REFERENCE_SETS} {reason}"
        for d in divergences for reason in [with_reference_sets_error(d)] if reason
    ]
    if problems:
        for problem in problems:
            print(f"ERROR: [build_golden_corpus] {problem}", file=sys.stderr)
        return 1
    # Classify's own keys are rewritten from scratch, in one order, so a key a previous run
    # wrote with other inputs (another `--field` list) does not survive this one.
    for key in CLASSIFY_KEYS:
        target.pop(key, None)
    target["divergences"] = divergences
    target["vep_rs_only_tuples"] = extra_rust
    target["divergence_summary"] = dict(Counter(d["expected_divergence"] for d in divergences))
    target["vep_rs_only_summary"] = dict(Counter(t["expected_divergence"] for t in extra_rust))
    field_note = ""
    if args.field:
        # Field values are compared only on keys whose consequence sets agree; a key whose
        # terms differ is already a documented (or unexplained) divergence above and the
        # golden test skips its fields, which is where the start co-emission keys land (VEP
        # prints `Met1?` wherever its start_lost predicate holds, vep-rs the peptide change
        # wherever it keeps start_retained_variant). A field difference on an agreeing key is
        # documented when `classify_field_divergence` knows its class, else it is an
        # `unexplained_residual`: the golden test applies an entry whatever its class, and
        # `test_golden_scripts.py` refuses a committed manifest that carries the class.
        # The reference's STRAND, FLAGS, DISTANCE and BAM_EDIT and vep-rs's BAM_EDIT ride along
        # with the compared fields: `classify_field_divergence` reads them.
        # The score files the run's SpliceAI plugin read, and the rows' SpliceAI fields with them.
        spliceai_table = None
        if args.spliceai_table:
            cache = corpus_cache_dir(corpus, manifest)
            spliceai_table = format_comparator().SpliceAiTable(
                [Path(path) for path in args.spliceai_table], cache if (cache / "transcripts").is_dir() else None,
            )
        spliceai_fields = SPLICEAI_ROW_FIELDS if spliceai_table is not None else ()
        perl_fields = load_fields(Path(args.vep_default), [*args.field, "STRAND", "MINIMISED", "FLAGS", "DISTANCE", "BAM_EDIT", *spliceai_fields], masked)
        rust_fields = load_fields(Path(args.vep_rs_output), [*args.field, "BAM_EDIT", *spliceai_fields])
        edited_models = indel_edited_models(corpus_cache_dir(corpus, manifest))
        codon_tables = organelle_codon_tables(corpus)
        downstream = any(name in DOWNSTREAM_FIELDS for name in args.field)
        coding_sequences = downstream_sequences(corpus) if downstream else None
        coding_models = downstream_coding_models(corpus) if downstream else None
        partial_codon_models = partial_terminal_codon_models(corpus_cache_dir(corpus, manifest))
        intron_ends = intron_end_gaps(corpus_cache_dir(corpus, manifest))
        intron_starts = intron_start_gaps(corpus_cache_dir(corpus, manifest))
        two_windows = None
        if any(name.startswith(MAXENTSCAN_FIELD_PREFIX) for name in args.field):
            two_windows = two_window_positions(corpus_cache_dir(corpus, manifest))
        field_divergences: list[dict] = []
        minimal = "--minimal" in flags
        for key, pvals in sorted(perl_fields.items()):
            rvals = rust_fields.get(key)
            if rvals is None:
                continue
            indices = records_for(key[:4])
            # A multi-allelic record `--minimal` split, its rows marked MINIMISED (a record none
            # of whose alleles trims is annotated unsplit by both engines): the split hands
            # every plugin the rejoined record in Ensembl VEP, the allele's own span and
            # alleles in vep-rs.
            minimal_split_record = (minimal and pvals["MINIMISED"] == "1"
                                    and any("," in records[i]["vcf"].split("\t")[4] for i in indices))
            several_entries = None
            for name in args.field:
                if not fields_agree(name, pvals[name], rvals[name], ORDER_ONLY in masks):
                    genomic_span = None
                    if name in DOWNSTREAM_FIELDS:
                        genomic_span = row_record_span(key[0], (records[i]["vcf"] for i in indices))
                    klass, second = classify_field_divergence(
                        name, key, pvals[name], rvals[name], pvals["STRAND"], flags, edited_models, codon_tables,
                        coding_sequences, pvals, rvals, minimal_split_record, coding_models, genomic_span,
                        partial_codon_models, intron_ends, intron_starts, two_windows, spliceai_table,
                    ) or (UNEXPLAINED, None)
                    entry = {
                        "location": key[0], "allele": key[1], "feature": key[2], "feature_type": key[3],
                        "record_indices": indices,
                        "consequence_set": key[4], "field": name, "vep": pvals[name], "vep_rs": rvals[name],
                        "expected_divergence": klass,
                    }
                    if second is not None:
                        entry[SECOND_EXPECTED_DIVERGENCE] = second
                    field_divergences.append(entry)
                    if klass == SPLICEAI_MULTI_ENTRY_CELL:
                        several_entries = entry
            if several_entries is not None:
                # The row's value under the names the VCF and JSON formats give it.
                field_divergences += spliceai_other_format_entries(several_entries, pvals, rvals, args.field)
        target["fields_compared"] = list(args.field)
        target["field_divergences"] = field_divergences
        target["field_divergence_summary"] = dict(Counter(d["expected_divergence"] for d in field_divergences))
        field_note = f"; field divergences {len(field_divergences)} ({target['field_divergence_summary']})"
    manifest.update(run_options(args))
    write_manifest(corpus, manifest)
    print(f"divergent VEP rows {len(divergences)} ({target['divergence_summary']}); "
          f"vep-rs-only tuples {len(extra_rust)} ({target['vep_rs_only_summary']}); "
          f"reference-skipped records {sorted(skipped)}{field_note}")
    return 0


# ------------------------------------------------------------------------------ check
def sample_column_errors(corpus: Path, manifest: dict) -> list[str]:
    """One line per way `variants.vcf` fails to carry the manifest's kept samples: a `#CHROM`
    line whose columns are not the eight fixed ones, `FORMAT` and the manifest's `samples`, or
    a record whose column count differs from that header's; a file that cannot be read is
    refused by `vcf_lines`."""
    keep = manifest["keep_samples"]
    expected = [*FIXED_COLUMNS, "FORMAT", *manifest["samples"]]
    errors: list[str] = []
    has_header = False
    for is_header, f in vcf_lines(corpus / "variants.vcf"):
        if is_header:
            has_header = True
            if f != expected:
                errors.append(f"variants.vcf header names {' '.join(f[8:]) or 'no sample columns'}, the manifest keeps FORMAT {' '.join(manifest['samples'])}")
        elif len(f) != 9 + keep:
            errors.append(f"variants.vcf record {' '.join(f[:5])} has {len(f)} columns, the header {9 + keep}")
    if not has_header:
        errors.append("variants.vcf has no #CHROM header line")
    return errors


def genotype_coverage(
    suites: list[tuple[str, Path, Path]], manifest: dict, required: list[str]
) -> tuple[set[str], set[str], list[str]]:
    """The genotype classes the suite inputs show among the manifest's kept samples, the ones
    the corpus records show, and those of the first set or of `required` the corpus lacks."""
    keep = manifest["keep_samples"]
    observed: set[str] = set()
    for _sid, _ref, inp in suites:
        samples: list[str] = []
        for is_header, f in vcf_lines(inp, 9 + keep):
            if is_header:
                samples = kept_samples(f, keep, inp, manifest["samples"], "the manifest")
            else:
                observed |= record_genotype_classes(f, keep)
        if not samples:
            raise failure(f"{inp} has no #CHROM header line naming its samples")
    covered = {c for rec in manifest["records"] for c in record_genotype_classes(rec["vcf"].split("\t"), keep)}
    return observed, covered, sorted((observed | set(required)) - covered)


def cmd_check(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    errors: list[str] = []
    if args.matrix is not None:
        matrix = Path(args.matrix)
        items = matrix_items(matrix)
        for owner, names in covers_lists(manifest):
            errors += names_outside_matrix(names, items, matrix, owner)
    suites = [parse_suite_arg(s) for s in args.suite]
    family_specs = parse_family_columns(args.family_columns)
    required = genotype_class_names(args.require_genotype_class)
    if required and not manifest.get("keep_samples"):
        raise failure("--require-genotype-class wants a corpus built with --keep-samples")
    by_name = record_name_index(manifest["records"])
    records = manifest["records"]

    def is_corpus_row(r: Row) -> bool:
        chrom = r.location.split(":", 1)[0]
        return any(records[i]["suite"] == r.suite and records[i]["vcf"].split("\t", 1)[0] == chrom
                   for i in by_name.get(r.uploaded_variation, []))

    observed: set[str] = set()
    strata_observed: set[str] = set()
    strata_covered: set[str] = set()
    for sid, ref, inp in suites:
        try:
            # A family stratum ends in the row's variant class, which is read with the input.
            substitutions = substitution_alts(inp) if family_specs else None
        except SystemExit:
            # The refusal that ends the run comes after the errors found before it.
            for e in errors:
                print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
            raise
        for r in iter_rows(ref, sid, family_specs, substitutions):
            observed.add(r.consequence_set)
            if family_specs:
                strata_observed.add(r.family_stratum)
                if is_corpus_row(r):
                    strata_covered.add(r.family_stratum)
    covered = {c for rec in records for c in rec["combinations"]}
    missing = sorted(observed - covered)
    print(f"observed {len(observed)}, covered by corpus {len(covered & observed)}, missing {len(missing)}")
    for m in missing[:50]:
        print("  MISSING", m)
    missing_strata = sorted(strata_observed - strata_covered)
    if family_specs:
        print(f"family strata observed {len(strata_observed)}, covered by corpus {len(strata_covered)}, "
              f"missing {len(missing_strata)}")
        for m in missing_strata[:50]:
            print("  MISSING family stratum", m)
    missing_classes: list[str] = []
    if manifest.get("keep_samples"):
        try:
            errors += sample_column_errors(corpus, manifest)
            observed_classes, covered_classes, missing_classes = genotype_coverage(
                suites, manifest, [*manifest.get("required_genotype_classes", []), *required]
            )
        except SystemExit:
            # The refusal that ends the run comes after the errors found before it.
            for e in errors:
                print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
            raise
        print(f"genotype classes observed {len(observed_classes)}, covered by corpus {len(covered_classes)}, "
              f"missing {len(missing_classes)}")
        for m in missing_classes:
            print("  MISSING genotype class", m)
    missing_record_classes: list[str] = []
    if manifest.get("record_class_key"):
        classes_of, _samples, class_seen = record_class_pass(suites, manifest.get("seed", 0), manifest["record_class_key"])
        covered_record_classes: set[str] = set()
        for (sid, name, location), classes in classes_of.items():
            chrom = location.split(":", 1)[0]
            if any(records[i]["suite"] == sid and records[i]["vcf"].split("\t", 1)[0] == chrom
                   for i in by_name.get(name, [])):
                covered_record_classes |= classes
        missing_record_classes = sorted(set(class_seen) - covered_record_classes)
        print(f"record classes observed {len(class_seen)}, covered by corpus {len(covered_record_classes & set(class_seen))}, "
              f"missing {len(missing_record_classes)}")
        for m in missing_record_classes:
            print("  MISSING record class", m)
    for e in errors:
        print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
    return 1 if missing or missing_strata or missing_classes or missing_record_classes or errors else 0


def cmd_covers(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    target = flag_set_named(manifest, args.flag_set) if args.flag_set is not None else manifest
    covers = list(target.get("covers") or [])
    if args.add:
        added = covers_names(args.add)
        matrix = Path(args.matrix)
        errors = names_outside_matrix(added, matrix_items(matrix), matrix, covers_owner(args.flag_set))
        for e in errors:
            print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
        if errors:
            return 1
        covers = list(dict.fromkeys([*covers, *added]))
        target["covers"] = covers
        write_manifest(corpus, manifest)
    if args.list or not args.add:
        for name in covers:
            print(name)
    return 0


def cmd_flag_set(args: argparse.Namespace) -> int:
    corpus = Path(args.corpus)
    manifest = read_manifest(corpus)
    sets = manifest.setdefault("flag_sets", [])
    if args.list:
        for fs in sets:
            columns = [
                fs["name"],
                " ".join(fs.get("flags") or []),
                fs.get("input", DEFAULT_INPUT),
                str(bool(fs.get("stdin"))).lower(),
                str(bool(fs.get("unordered"))).lower(),
                ",".join(fs.get("formats") or OUTPUT_FORMATS),
                fs.get("expected_dir", f"expected_{fs['name']}"),
                ",".join(fs.get("covers") or []),
            ]
            print("\t".join(columns))
        return 0
    if not args.name:
        print("ERROR: [build_golden_corpus] flag-set wants --name (or --list)", file=sys.stderr)
        return 1
    if not is_flag_set_name(args.name):
        print(f"ERROR: [build_golden_corpus] flag set name {args.name!r}: letters, digits, _ and -, "
              "and not the manifest's own run, default", file=sys.stderr)
        return 1
    formats = [f.strip() for f in (args.format or ",".join(OUTPUT_FORMATS)).split(",") if f.strip()]
    bad = [f for f in formats if f not in OUTPUT_FORMATS]
    if bad or not formats:
        print(f"ERROR: [build_golden_corpus] --format wants a non-empty subset of {','.join(OUTPUT_FORMATS)}, "
              f"got {args.format!r}", file=sys.stderr)
        return 1
    flags = list(args.flag or [])
    if args.stdin and not names_format(flags):
        # Ensembl VEP detects no format on STDIN and dies without --format, so a stdin run's
        # reference exists only for a flag list that names one.
        print("ERROR: [build_golden_corpus] a stdin run names its input format: add --flag=--format --flag=<format>",
              file=sys.stderr)
        return 1
    path_inside_corpus(args.input, "--input")
    if not (corpus / args.input).is_file():
        print(f"ERROR: [build_golden_corpus] input {args.input} is not a file in {corpus}", file=sys.stderr)
        return 1
    expected_dir = path_inside_corpus(args.expected_dir or f"expected_{args.name}", "--expected-dir")
    covers = covers_names(args.covers)
    if covers:
        matrix = Path(args.matrix)
        errors = names_outside_matrix(covers, matrix_items(matrix), matrix, covers_owner(args.name))
        for e in errors:
            print(f"ERROR: [build_golden_corpus] {e}", file=sys.stderr)
        if errors:
            return 1
    written = {
        "name": args.name,
        "flags": flags,
        "input": args.input,
        "stdin": bool(args.stdin),
        "unordered": bool(args.unordered),
        "formats": formats,
        "expected_dir": expected_dir,
        "covers": covers,
    }
    if args.perl_mode != "offline":
        written["perl_mode"] = args.perl_mode
    for i, fs in enumerate(sets):
        if fs["name"] == args.name:
            # The options set their keys where they sit and the entry's other keys (`masks`,
            # `reference_warnings`, any hand-written key) stay. The divergence keys go, since
            # they describe the run before the rewrite (`classify --flag-set` rewrites them), and
            # so does `perl_mode`, which `written` carries only for a mode other than offline.
            entry = {k: v for k, v in fs.items() if k not in CLASSIFY_KEYS and k != "perl_mode"}
            entry.update(written)
            sets[i] = entry
            break
    else:
        entry = written
        sets.append(entry)
    write_manifest(corpus, manifest)
    print(f"flag set {args.name}: {len(entry['flags'])} flags, input {entry['input']}"
          f"{' on stdin' if entry['stdin'] else ''}{', unordered' if entry['unordered'] else ''}"
          f"{', reference ' + entry['perl_mode'] if 'perl_mode' in entry else ''}, "
          f"formats {','.join(formats)}, expected {entry['expected_dir']}, "
          f"covers {len(covers)}")
    return 0


def add_run_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--flag", action="append",
                        help="a command-line flag both engines ran with, written --flag=--hgvs (repeatable), recorded as `flags`")
    parser.add_argument("--fasta-name", default=None,
                        help="file name of the gzipped reference beside the manifest, recorded as `fasta`")


MATRIX_HELP = "the parity matrix whose `item` column every `covers` name must be in"
DEFAULT_MATRIX = "tests/parity/matrix.tsv"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--release", required=True)
    s.add_argument("--assembly", required=True)
    s.add_argument("--species", default=None,
                   help=f"the Ensembl species of the corpus, written to the manifest (default: an existing manifest's, else {DEFAULT_SPECIES})")
    s.add_argument("--suite", action="append", required=True, help="sid=reference.txt=input.vcf[.gz]")
    s.add_argument("--out", required=True)
    s.add_argument("--k", type=int, default=3)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--max-records", type=int, default=1500)
    s.add_argument("--multi-allelic-per-suite", type=int, default=4,
                   help="multi-allelic sequence records sampled per suite and stratum (SNV-only, indel)")
    s.add_argument("--max-span", type=int, default=1_000_000,
                   help="prefer exemplars no wider than this many bases; 0 disables the preference")
    s.add_argument("--focus-term", action="append",
                   help="a consequence term; combinations containing it take --focus-k exemplars")
    s.add_argument("--focus-k", type=int, default=None, help="exemplars per focus combination")
    s.add_argument("--include", action="append", metavar="SID=NAME",
                   help="a record of suite SID's input, by its ID or VEP's name for an ID-less record, carried whatever "
                        "the sampling chooses (repeatable), recorded as `included_records`")
    s.add_argument("--keep-records", action="store_true",
                   help="keep every record of the manifest being regenerated, resolved and rowed again from its suite; "
                        "the count is recorded as `kept_records`")
    s.add_argument("--cache", default=None,
                   help="the JSON cache the fixture cache is pruned from; exemplars whose transcripts add fewer "
                        "bytes to it come first, ahead of the flag and position preferences")
    s.add_argument("--flank", type=int, default=5000,
                   help="the transcript flank the pruner keeps around each record, in bases (its --flank)")
    s.add_argument("--covers", action="append",
                   help="a parity-matrix row the corpus exercises (repeatable), recorded as `covers`; "
                        "required unless the manifest being regenerated already carries a list")
    s.add_argument("--family-columns", default=None,
                   help="Extra columns of the family key, comma-separated: NAME stratifies by presence, NAME= by value, "
                        "custom_match=EXACT/OVERLAP by how a --custom file's exact and overlap sources answered, "
                        "A/B by whether two columns agree; every stratum seen keeps an exemplar")
    s.add_argument("--keep-samples", type=int, default=0,
                   help="keep FORMAT and the first N sample columns of every record and of the #CHROM header "
                        "(default: none), and sample records by the genotype classes those samples show")
    s.add_argument("--genotype-exemplars", type=int, default=8,
                   help="records showing each genotype class the inputs show among the kept samples, when they have that many")
    s.add_argument("--require-genotype-class", action="append", metavar="CLASS",
                   help=f"a genotype class the corpus must carry, shown by the inputs or not (repeatable; one of {', '.join(GENOTYPE_CLASSES)})")
    s.add_argument("--record-classes", default=None, metavar="KEY",
                   help=f"class every record by its reference rows under this key (one of {', '.join(RECORD_CLASS_KEYS)}) "
                        "and add records until every class the references show has --record-class-exemplars records")
    s.add_argument("--record-class-exemplars", type=int, default=2,
                   help="records showing each record class the references show, when they have that many")
    s.add_argument("--matrix", default=None, help=MATRIX_HELP)
    add_run_options(s)
    s.set_defaults(fn=cmd_select)
    c = sub.add_parser("classify")
    c.add_argument("--corpus", required=True)
    c.add_argument("--vep-default", required=True, help="VEP default-format output for variants.vcf")
    c.add_argument("--vep-rs-output", required=True, help="vep-rs default-format output for variants.vcf")
    c.add_argument("--field", action="append",
                   help="an Extra field (HGVSc, HGVSp) compared on every key whose consequence sets agree")
    c.add_argument("--flag-set", default=None,
                   help="write the divergence keys into this flag set's entry instead of the manifest's top level")
    c.add_argument("--spliceai-table", action="append", metavar="FILE",
                   help="a score file the run's SpliceAI plugin read, its tabix index beside it (repeatable: the snv= file "
                        "and the indel= file), whose records the class spliceai_multi_entry_cell reads; the gene of a row "
                        "that writes no SYMBOL is read from the corpus's cache")
    add_run_options(c)
    c.set_defaults(fn=cmd_classify)
    k = sub.add_parser("check")
    k.add_argument("--corpus", required=True)
    k.add_argument("--suite", action="append", required=True)
    k.add_argument("--matrix", default=None, help=MATRIX_HELP)
    k.add_argument("--family-columns", default=None,
                   help="the family key (see select); every stratum the suites show must have an exemplar")
    k.add_argument("--require-genotype-class", action="append", metavar="CLASS",
                   help="a genotype class the corpus must carry, beside the manifest's required_genotype_classes (repeatable)")
    k.set_defaults(fn=cmd_check)
    v = sub.add_parser("covers")
    v.add_argument("corpus")
    v.add_argument("--add", action="append", help="a parity-matrix row to add (repeatable)")
    v.add_argument("--list", action="store_true", help="print the list, one name per line")
    v.add_argument("--flag-set", default=None, help="the flag set whose list is read or extended")
    v.add_argument("--matrix", default=DEFAULT_MATRIX, help=f"{MATRIX_HELP} (default: {DEFAULT_MATRIX})")
    v.set_defaults(fn=cmd_covers)
    f = sub.add_parser(
        "flag-set",
        description="Write one flag_sets entry of the corpus manifest: a run of the corpus under its own flags. "
                    "A set that exists has the keys the options below set rewritten and keeps every other key "
                    "of its entry (its hand-written `masks` and `reference_warnings` among them); its divergence "
                    "keys go, so `classify --flag-set` runs after it.",
    )
    f.add_argument("corpus")
    f.add_argument("--name", default=None, help="the set's name; `expected_<name>` holds its outputs unless --expected-dir")
    f.add_argument("--flag", action="append", help="a command-line flag of the run, written --flag=--pick (repeatable)")
    f.add_argument("--input", default=DEFAULT_INPUT, help=f"the input file inside the corpus (default: {DEFAULT_INPUT})")
    f.add_argument("--stdin", action="store_true", help="pipe the input on standard input and pass -i STDIN")
    f.add_argument("--unordered", action="store_true",
                   help="the reference's row order within a record is Ensembl VEP's hash order (--per_gene gene "
                        "buckets, --pick_allele allele buckets), so the runner compares the entry set and not the order")
    f.add_argument("--perl-mode", default="offline", choices=PERL_MODES,
                   help="how the set's reference was generated: offline (Ensembl VEP on its cache; the default, "
                        "which writes no key, so the set takes the manifest's own mode, offline when the manifest "
                        "names none), database (connected to the Ensembl database, for a column Ensembl VEP writes only with "
                        "one; the default and tab headers then compare by their column descriptions alone), or "
                        "gtf or gff (reading its transcripts from the file given to --gtf or --gff, whose SOURCE "
                        "column and header lines the golden harness sets aside)")
    f.add_argument("--format", default=None, help=f"the output formats compared, comma-separated (default: all of {','.join(OUTPUT_FORMATS)})")
    f.add_argument("--expected-dir", default=None, help="the directory of the set's expected outputs")
    f.add_argument("--covers", action="append", help="a parity-matrix row the set exercises (repeatable)")
    f.add_argument("--matrix", default=DEFAULT_MATRIX, help=f"{MATRIX_HELP} (default: {DEFAULT_MATRIX})")
    f.add_argument("--list", action="store_true", help="print the sets, one tab-separated line each")
    f.set_defaults(fn=cmd_flag_set)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
