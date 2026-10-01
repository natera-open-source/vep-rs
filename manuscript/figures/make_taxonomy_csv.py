"""Build manuscript/data/discordance_taxonomy.csv.

Emits one row per (suite, discordance class) pair documenting how many
tuples that class accounts for in the suite, classified as a Perl VEP defect
that vep-rs corrects (excluded from adjusted F1), a representation difference
in which neither engine is wrong, or an open vep-rs gap.

Every count is measured at the vep-rs binary MEASURED_BINARY names, in the sweep
VEP_RS_SWEEP_ID in generate_figures.py names, and the ``provenance`` column
records that binary. The SNP/indel classes come from ``adjusted.excluded_categories``
in each suite's ``report/summary.json``, whose per-bucket ``{count, excluded}``
pair supplies ``tuples_per_suite`` and ``masked_pairs_per_suite``; the SV classes
come from the ``adjusted.excluded_*`` counters in the two SV suites'
``report/concordance_report.json`` (the transcript-selection class from
``excluded_transcript_selection_rust``, the layer's own size, which the report
carries beside the union ``excluded_rust_extra``). The ``tuples_per_suite`` entries
below carry the pinned sweep's measured counts.
"""

from __future__ import annotations

import csv
from pathlib import Path

# Every count below is read from this binary's clone reports unless the entry names a
# later `measured_binary`; the CSV's `provenance` column carries whichever applies.
MEASURED_BINARY = "23a294cf"

# One counting convention for every row. A pair class (the SNP/indel classes and the
# <CNV:TR> class) counts PAIRS: one per (Location, Allele, Feature) the two engines both
# wrote with different consequence sets. A one-sided SV class counts TUPLES removed from
# the adjusted denominators by its predicates on BOTH engines' sides, on the report's
# union basis, each tuple attributed to one class: a Perl tuple the transcript-selection
# Perl arm and the cross-chromosome mask both claim is counted in the cross-chromosome
# class (the report writes the overlap as
# `transcript_selection_cross_chromosome_overlap_perl`).

# The two SV classes keyed `sv_grch37` / `sv_grch38` are per-assembly aggregates over
# the assembly's whole 16-file SV set (13 synthetic + 3 real-world files), not over the
# real-world files alone.

# `category` is a coarse class grouping, not a published vocabulary: the supplementary
# class table prints `mask_status` and no category column, and a class carrying
# `perl_defect_*` here need not meet the main text's stricter Perl-defect test, that
# Ensembl VEP contradicts itself. Its prefix says whether the class is a Perl VEP defect
# (`perl_defect_`) or a difference in which neither engine contradicts itself or its input
# (`representation_difference_`); its suffix says how far the adjusted-F1 mask reaches the
# class, as measured, and is the class-level mask state. Permitted values:
#   - "perl_defect_fixed": a Perl VEP defect vep-rs corrects; every pair is excluded.
#   - "perl_defect_partially_masked": a Perl VEP defect masked only where an Ensembl VEP
#     self-contradiction accounts for the pair; the rest stay charged to vep-rs.
#   - "representation_difference_fully_masked": a symmetric filter pairs every divergence
#     of the class and excludes it from both sides.
#   - "representation_difference_partially_masked": the same, reaching part of the class.
#   - "vep_rs_gap": vep-rs differs from Perl and nothing is masked.
# `mask_status` carries the machine-checkable per-suite mask state for the
# counted-vs-filtered reconciliation:
#   - "fully_masked":     every swap pair in the class is excluded from adjusted F1.
#   - "partially_masked": a named filter excludes some pairs; the rest remain.
#   - "not_masked":       nothing excluded.
MASK_STATUS_BY_CATEGORY = {
    "perl_defect_fixed": "fully_masked",
    "perl_defect_partially_masked": "partially_masked",
    "representation_difference_fully_masked": "fully_masked",
    "representation_difference_partially_masked": "partially_masked",
    "vep_rs_gap": "not_masked",
}
TAXONOMY = [
    {
        # The layer's own size (`adjusted.excluded_transcript_selection_rust`), which
        # compare_sv_concordance.py writes beside the union of the vep-rs-side layers
        # (`excluded_rust_extra`, transcript selection with <CNV:TR>) that the adjusted
        # denominator subtracts; the two differ by the <CNV:TR> layer's pairs, and the
        # per-file split (`excluded_transcript_selection_by_file`) puts every pair of this
        # class in the gnomAD-SV file of each assembly. A report without the field yields
        # the same figure as the union less the <CNV:TR> count, the two layers touching
        # disjoint files.
        #
        # Two predicates make the class: `filter_intended_divergences` (the local allele of
        # a record above --max_sv_size, or of an unsupported symbolic type Perl flags
        # `vep_skip` the same way, `CPX` and `CTX`; both engines' sides) and
        # `filter_giant_breakend_mate_divergences` (the mate allele of a breakend above
        # --max_sv_size that Perl never wrote, vep-rs side). The count is
        # `excluded_transcript_selection_rust` + `excluded_transcript_selection_mate_rust` +
        # `excluded_transcript_selection_perl` less
        # `transcript_selection_cross_chromosome_overlap_perl` (10,934 + 56 + 0 - 0 on
        # GRCh37; 183,964 + 310 + 1,823 - 1,823 on GRCh38, every Perl-side tuple there being
        # a cross-chromosome one). The per-file split puts every tuple in the gnomAD-SV file
        # of each assembly.
        "class": "sv_nondeterministic_transcript_selection",
        "category": "perl_defect_fixed",
        "measured_binary": "0db781b4",
        "perl_source": "ensembl-vep VEP/AnnotationSource.pm:get_all_features_by_InputBuffer",
        # `annotate_batch` annotates every variant of a batch against the complete
        # transcript set of its chromosome, a breakend of any span included, which makes
        # its output batch-independent where Perl's InputBuffer-scoped lookup above
        # --max_sv_size is not. transcript_index.rs is the index that lookup reads, not
        # the dispatch.
        "rust_source": "vep-cli/src/runner.rs:annotate_batch",
        "tuples_per_suite": {
            "sv_grch37": 10990,
            "sv_grch38": 184274,
        },
    },
    # The taxonomy lists exactly the classes that move an adjusted denominator; a divergence
    # excluded from nothing has no entry. A class whose measured count reaches zero is
    # deleted rather than emitted as a zero row.
    {
        # Ensembl VEP annotates a structural variant against transcripts on a different
        # chromosome: a `<DEL>` on chromosome 21 is annotated against transcripts on
        # chromosome 22, output that raw F1 charges to vep-rs. It fires on both assemblies
        # and the mechanism is not established: the GRCh38 gnomAD SV VCF carries
        # `CHR2=chr21` against an unprefixed CHROM while the GRCh37 file carries `CHR2=21`,
        # yet the class fires on GRCh37 too, so contig naming is not the whole story, which
        # is why `perl_source` names the entry point and no mechanism.
        #
        # `filter_cross_chromosome_divergences` excludes it on an objectively checkable
        # criterion, whether the transcript is on the variant's chromosome, judged against
        # the vep-rs JSON cache both engines annotated against (`--vep-rs-cache`). Without
        # the cache the mask is inactive and the report records
        # `cross_chromosome_check_ran: false`. The filter raises rather than excluding if
        # vep-rs produces the shape, since that would be a vep-rs bug.
        #
        # The mechanism: a neighbouring breakend's mate loads the other chromosome's cache
        # regions (AnnotationSource.pm:249-257) and the input buffer's interval tree matches
        # transcripts to records by numeric coordinate with no chromosome test
        # (InputBuffer.pm:284-330, 345-366); a non-breakend record then gets an overlap
        # allele without the breakend distance gate (StructuralVariationOverlap.pm:64-72).
        # Alone, every such record writes `intergenic_variant`, which is what vep-rs wrote,
        # so the class has a vep-rs side too: `filter_cross_chromosome_orphan_intergenic`
        # removes that row where every Perl tuple of the record was removed
        # (`excluded_cross_chromosome_rust`). The count is `excluded_cross_chromosome_perl`
        # + `excluded_cross_chromosome_rust` (1,838 + 154; 22,481 + 469); the per-file
        # counts of the vep-rs side are 155 and 518, the difference being records the
        # 1000 Genomes SV file repeats at the same coordinates, where both engines write
        # `intergenic_variant`.
        "class": "sv_cross_chromosome_transcript_annotation",
        "category": "perl_defect_fixed",
        "measured_binary": "0db781b4",
        "perl_source": "ensembl-vep VEP/AnnotationSource.pm:get_all_regions_by_InputBuffer (mate regions, 249-257); VEP/InputBuffer.pm:get_overlapping_vfs (numeric interval tree, 284-330)",
        "rust_source": "scripts/validation/compare_sv_concordance.py:filter_cross_chromosome_divergences, filter_cross_chromosome_orphan_intergenic",
        "tuples_per_suite": {
            "sv_grch37": 1992,
            "sv_grch38": 22950,
        },
    },
    {
        # A `<NON_REF>` record (a gVCF reference block) has no Sequence Ontology term, so Perl
        # keeps it with `vep_skip`, loads no cache region for it and annotates it against
        # whatever regions its batch loaded; alone in a batch it writes one row,
        # `intergenic_variant`. vep-rs writes that one row for every such record.
        # `filter_non_ref_batch_divergences` excludes Perl's Transcript tuples on a record
        # whose vep-rs set is exactly the intergenic row, and that row where Perl wrote
        # Transcript tuples and no intergenic row. The count is
        # `excluded_non_ref_batch_perl` + `excluded_non_ref_batch_rust`: 1,838 + 180 on
        # GRCh37 (Perl also wrote the intergenic row on 20 of the 200 records, where the
        # two engines agree) and 5,403 + 200 on GRCh38, all in `10_special_alleles`.
        "class": "sv_non_ref_batch_annotation",
        "category": "perl_defect_fixed",
        "measured_binary": "0db781b4",
        "perl_source": "ensembl-vep VEP/Parser/VCF.pm:create_StructuralVariationFeatures (477-481, vep_skip for a type without a Sequence Ontology term); VEP/AnnotationSource.pm:get_all_regions_by_InputBuffer (238)",
        "rust_source": "vep-cli/src/vcf_parser.rs:classify_symbolic_alt (a reference block yields one intergenic_variant row); scripts/validation/compare_sv_concordance.py:filter_non_ref_batch_divergences",
        "tuples_per_suite": {
            "sv_grch37": 2018,
            "sv_grch38": 5603,
        },
    },
    {
        # Perl builds a Transcript row for every transcript within 5 kb of a point
        # breakend's mate on the mate chromosome and then evaluates every positional
        # predicate on it with the LOCAL variation feature; only `feature_truncation` reads
        # the mate. The row carries the region term the local coordinate lands on when
        # applied numerically to the mate transcript, or `intergenic_variant` on a
        # Transcript row when nothing fires. vep-rs evaluates the mate row at the mate
        # coordinate. `filter_breakend_mate_local_read_pairs` derives both sets from the
        # cache (`breakend_mate_context.py`) and excludes a pair only when Perl's set is the
        # local read, vep-rs's the mate read, and the two differ; the count is pairs
        # (`excluded_breakend_mate_local_read_perl`, equal on the vep-rs side): 1,181
        # cross-chromosome + 226 same-chromosome on GRCh37, 6,840 + 418 on GRCh38, in
        # `09_breakends`, `13_vcf45_features` and the gnomAD-SV file of GRCh38.
        "class": "sv_breakend_mate_context_local_read",
        "category": "perl_defect_fixed",
        "measured_binary": "0db781b4",
        "perl_source": "ensembl-variation BaseVariationFeatureOverlapAllele.pm:get_all_OverlapConsequences (257, 273: the local feature reaches every predicate); Utils/VariationEffect.pm:feature_truncation (358, the one breakend read); StructuralVariationOverlap.pm:new (73-87, the mate row)",
        "rust_source": "vep-effects/src/sv/breakend.rs:calculate_paired_mate (the mate evaluated as a point breakend on its own chromosome); scripts/validation/breakend_mate_context.py:mate_side_derivation",
        "tuples_per_suite": {
            "sv_grch37": 1407,
            "sv_grch38": 7258,
        },
    },
    {
        # Perl expands a <CNV:TR> record into its literal repeat sequence, trims the
        # reference over the run, and reads a gain as a whole-unit insertion at the run's
        # end and a loss as a deletion ending at the run's last base, then runs its
        # point-variant predicates over that allele; vep-rs keeps the allele symbolic,
        # classifies it as VariantClass::TandemRepeat and reads the same record as a
        # structural gain (insertion path, feature_elongation) or loss (deletion path,
        # feature_truncation) by comparing INFO/RB with the reference run. Neither engine
        # contradicts itself or its input, so this is a representation difference, not a
        # Perl defect. Synthetic 07_cnv_repeat only, whose records sit on tandem repeats
        # the chromosome 21 reference carries.
        #
        # `filter_cnv_tr_expansion_divergences` pairs the two readings of one record and
        # transcript and excludes the pair from both sides when they agree on every other
        # term; raw F1 is identical with and without the mask, which a swap-mode filter
        # requires by construction. `tuples_per_suite` is the class size (the report's
        # `cnvtr_swap_pairs_total`); the filter reaches every pair on both assemblies (the
        # `excluded_cnvtr` counters equal the totals), so the class is fully masked and the
        # row carries no separate masked count.
        "class": "sv_cnv_tr_literal_expansion",
        "category": "representation_difference_fully_masked",
        "perl_source": "ensembl-vep VEP/Parser/VCF.pm:_expand_tandem_repeat_allele_string (382-431; called from create_StructuralVariationFeatures:520)",
        "rust_source": "vep-effects/src/sv/mod.rs:calculate_sv_consequences (TandemRepeat -> insertion::calculate for a gain, deletion::calculate for a loss, by INFO/RB against the reference run; allele kept symbolic; class from vep-cli/src/vcf_parser.rs:classify_symbolic_alt, alternate length from vcf_parser.rs:tandem_repeat_alt_bases)",
        "tuples_per_suite": {
            "sv_synth_grch37": 1368,
            "sv_synth_grch38": 2644,
        },
    },
    {
        # One class folded from the comparator's four splice buckets
        # (intronic_overcall_swap, ppt_overcall_swap, frameshift_splice_cascade_swap,
        # splice_lastwrite_swap), because the covered/uncovered split is measured across
        # the family rather than per bucket; sibling rows would present one measurement
        # four times and invite summing them. `tuples_per_suite` is the family total and
        # `masked_pairs_per_suite` is what the adjusted F1 removes: of the fourteen writes
        # in Perl's `_intron_effects` only `splice_region` assigns a return value a later
        # boundary intron can clobber, so the cited mechanism can drop exactly one term,
        # and the unmasked remainder add a polypyrimidine, acceptor, donor, fifth-base or
        # donor-region term whose flag is set once. The masked counts are the reports'
        # per-bucket {count, excluded} pairs summed over the family. A suite whose count
        # is 0 is omitted rather than emitted as a zero row.
        "class": "splice_family_swap",
        "category": "perl_defect_partially_masked",
        "perl_source": "BaseTranscriptVariationAllele.pm:_intron_effects:215 (splice_region overwrite; the other 13 writes are set-once)",
        "rust_source": "vep-effects/src/consequences.rs",
        # The SV suites reach the class through `filter_registry_swap_pairs`, which
        # applies the same excluding rule to the SV comparator's pairs
        # (`excluded_registry_swap_by_bucket.covered_splice_region_swap`).
        "tuples_per_suite": {
            "clinvar_grch37": 440,
            "clinvar_grch38": 1444,
            "sv_grch38": 10,
        },
        "masked_pairs_per_suite": {
            "clinvar_grch37": 440,
            "clinvar_grch38": 1444,
            "sv_grch38": 10,
        },
        "measured_binary_per_suite": {"sv_grch38": "0db781b4"},
    },
    {
        # The reports' excluded_categories.start_cooccurrence_swap per suite.
        "class": "start_cooccurrence_swap",
        "category": "perl_defect_fixed",
        "perl_source": "VariationEffect.pm:start_lost+start_retained_variant+_ins_del_start_altered",
        # The co-emission gate keeps the two mutually exclusive terms from being
        # emitted for one allele: on a sequence variant `perl_coding_terms` drops
        # `start_lost` (Perl's retention check has found the start codon intact), on a
        # structural deletion the SV path keeps `start_lost` alone (that check is
        # skipped for structural alleles); cited by symbol, not by line.
        "rust_source": "vep-effects/src/coding.rs:perl_coding_terms (start co-emission gate, sequence variants); vep-effects/src/sv/deletion.rs (structural alleles)",
        # The SV suites reach the class through `filter_registry_swap_pairs`
        # (`excluded_registry_swap_by_bucket.start_cooccurrence_swap`), on the report's union
        # basis: 23 on GRCh38 is 24 per file, one deletion
        # recurring in the 1000 Genomes and gnomAD SV files at the same coordinates.
        "tuples_per_suite": {
            "clinvar_grch37": 261,
            "clinvar_grch38": 928,
            "gnomad_grch37": 3,
            "gnomad_grch38": 77,
            "1kg_grch38": 1,
            "sv_grch37": 3,
            "sv_grch38": 23,
        },
        "measured_binary_per_suite": {"sv_grch37": "0db781b4", "sv_grch38": "0db781b4"},
    },
]


def main() -> None:
    out_dir = Path(__file__).resolve().parent.parent / "data"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "discordance_taxonomy.csv"

    rows: list[dict[str, str]] = []
    for entry in TAXONOMY:
        # Only the partially masked classes carry masked_pairs_per_suite; a fully masked
        # class's masked count is n_tuples by definition, and an unmasked suite's is 0.
        masked = entry.get("masked_pairs_per_suite", {})
        class_mask_status = MASK_STATUS_BY_CATEGORY[entry["category"]]
        # An entry, or one of its suites, re-measured at a later binary names it; the
        # others carry the pinned one.
        per_suite_binary = entry.get("measured_binary_per_suite", {})
        for suite, n_tuples in entry["tuples_per_suite"].items():
            provenance = "measured_" + per_suite_binary.get(
                suite, entry.get("measured_binary", MEASURED_BINARY)
            )
            if class_mask_status == "fully_masked":
                masked_pairs = str(n_tuples)
            elif suite in masked:
                masked_pairs = str(masked[suite])
            else:
                masked_pairs = "0"
            # Derived per suite from the masked count, because a partially masked class can
            # have suites its filter never matches, and that row would otherwise carry
            # `masked_pairs=0` with `partially_masked`. A consumer that needs one status per
            # class aggregates these rows, taking the strongest state any row carries.
            mask_status = class_mask_status
            if mask_status == "partially_masked" and masked_pairs == "0":
                mask_status = "not_masked"
            rows.append(
                {
                    "class": entry["class"],
                    "category": entry["category"],
                    "perl_source": entry["perl_source"],
                    "rust_source": entry["rust_source"],
                    "suite": suite,
                    "n_tuples": str(n_tuples),
                    "masked_pairs": masked_pairs,
                    "mask_status": mask_status,
                    "provenance": provenance,
                    "count_semantics": entry.get("count_semantics", "exact"),
                }
            )

    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            lineterminator="\n",
            fieldnames=[
                "class",
                "category",
                "perl_source",
                "rust_source",
                "suite",
                "n_tuples",
                "masked_pairs",
                "mask_status",
                "provenance",
                "count_semantics",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} rows -> {out_path}")


if __name__ == "__main__":
    main()
