# vep-effects

The consequence-calling engine of vep-rs. Given one variant and one transcript, it
decides which Sequence Ontology (SO) consequence terms the variant earns on that
transcript, following the Perl implementation of the Ensembl Variant Effect Predictor
(VEP) predicate by predicate. The `vep` binary ([vep-cli](../vep-cli/README.md)) is its
one caller in the workspace, and the crate hands back a `TranscriptConsequence`
(defined in [vep-core](../vep-core/README.md)) holding the sorted terms, the impact,
the display columns, the exon and intron ordinals and, when asked, the Human Genome
Variation Society (HGVS) notations.

## How it works

Follow a single-base substitution inside a coding exon. `vep-cli` has already selected
the transcripts whose span, widened by the up- and downstream distances, reaches the
variant, and it calls `calculate_consequences` once per transcript with an
`EffectsConfig`. The calculator first checks the variant class: a symbolic allele
(`<DEL>`, `<DUP>`, a breakend) is handed to `sv::calculate_sv_consequences`, which
dispatches on the class; this substitution is not symbolic, so it stays on the
sequence-variant path. `mapper::map_genomic_to_transcript` reads the transcript's exon
mapper pairs and returns a `TranscriptPosition::Coding` carrying the complementary DNA
(cDNA), coding sequence (CDS), protein and within-codon positions. The position pass
takes that classification to `coding::get_codon_change`, which cuts the reference codon
out of the transcript's translateable sequence, substitutes the alternate base,
translates both codons with the transcript's codon table and compares the residues: the
same residue yields `synonymous_variant`, a different one `missense_variant`, a new stop
`stop_gained`, a change in the start codon `start_lost`. Whole-span passes then add what
a single endpoint cannot see (a splice-region window in a neighbouring intron, an
untranslated region (UTR) the span runs into) and the suppression rules remove the pairs
Perl VEP never writes together. The list is sorted by severity rank, the impact is the
most severe term's, and `hgvs::generate_hgvs` adds HGVSc and HGVSp when `compute_hgvs`
is set. Back in `calculate_consequences`, `display::apply_display_fields` renders the
five display columns (`cDNA_position`, `CDS_position`, `Protein_position`,
`Amino_acids`, `Codons`) and, under `compute_exon_intron_numbers`, the `EXON` and
`INTRON` ordinals. The finished `TranscriptConsequence` returns to `vep-cli`, which
attaches it to the variant for the plugin phase and the output formatters.

## Modules

| Module                     | Role                                                                                                                                                                                                                                                                  |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/lib.rs`               | Declares the modules and re-exports `calculate_consequences`, `EffectsConfig` and `find_overlapping_transcripts`; holds the test-only transcript builders                                                                                                             |
| `src/consequences.rs`      | `EffectsConfig` and `calculate_consequences`: classifies each end of a sequence variant, applies the position terms, the whole-span splice and UTR passes and the suppression rules, sorts the terms by severity and sets the impact                                  |
| `src/mapper.rs`            | Maps a genomic position to a `TranscriptPosition` (upstream, downstream, UTR, coding, intron, non-coding exon) and a genomic span to CDS and peptide bounds or to a segment projection                                                                                |
| `src/coding.rs`            | Codon extraction and translation, the frameshift test on CDS-mapped lengths, inframe indel classification by peptide comparison, full-CDS peptide alleles and the ported Perl coding predicates (`perl_coding_terms`, `perl_hgvs_protein`)                            |
| `src/display.rs`           | The five display columns and the `EXON` and `INTRON` ordinals, rendered through the transcript mapper as Ensembl VEP renders them, apart from the term analysis                                                                                                       |
| `src/hgvs.rs`              | HGVSc and HGVSp strings and the `HGVS_OFFSET` of an insertion or deletion shifted 3' along the transcript                                                                                                                                                             |
| `src/overlap.rs`           | `find_overlapping_transcripts`: filters a transcript slice by chromosome and by the strand-aware extended region                                                                                                                                                      |
| `src/sv/mod.rs`            | `calculate_sv_consequences`: dispatches a structural variant on its class, reads the `start = end + 1` insertion pair and the `<CNV:TR>` contraction or gain, and holds the shared overlap predicates (exon, CDS, UTR, intron, polypyrimidine tract, mature microRNA) |
| `src/sv/deletion.rs`       | `<DEL>` and mobile-element deletions: `transcript_ablation` when the deletion contains the transcript, otherwise `feature_truncation` with the region terms; `stop_lost` over a stop-codon base, `start_lost` over a start-codon base when both ends are exonic       |
| `src/sv/duplication.rs`    | `<DUP>` and `<DUP:TANDEM>`: `transcript_amplification` when the duplication contains the transcript, `feature_elongation` when it lies inside the transcript and overlaps exonic sequence, the region terms otherwise                                                 |
| `src/sv/insertion.rs`      | `<INS>`, `<INS:ME:*>` and the `<CNV:TR>` gain: ranged, point and pair forms; `feature_elongation` inside the transcript over exonic sequence, `start_lost` over a start-codon base                                                                                    |
| `src/sv/inversion.rs`      | `<INV>`: region terms only, never `transcript_ablation`; `coding_transcript_variant` for a contained protein-coding transcript                                                                                                                                        |
| `src/sv/cnv.rs`            | `<CNV>`, `<CN0>` to `<CN9>` and `<CPX>`: `<CN0>` takes the deletion rules, `<CN2>` the duplication rules, every other copy number the overlap terms alone                                                                                                             |
| `src/sv/breakend.rs`       | Breakends in bracket and dot notation: `calculate` for the record's own breakpoint, `calculate_paired_mate` for the row on the mate chromosome; `feature_truncation` only when the allele's own breakend lies inside the transcript body and names a mate             |
| `src/concordance_tests.rs` | Test-only module of targeted tests, each asserting the Perl-correct terms for one divergence shape; the one `#[ignore]` test pins a Perl behaviour the crate does not reproduce                                                                                       |

## Public API

`calculate_consequences(variant, transcript, config) -> Option<TranscriptConsequence>`
is the entry point. It returns `None` when no term applies, as for a variant that
reaches neither the transcript nor its up- or downstream window, so a caller never sees
an empty row. A structural variant takes the same entry point and is routed internally.

`EffectsConfig` carries the run-wide settings. `vep-cli` fills it from the flags listed
in [docs/cli-reference.md](../../docs/cli-reference.md).

| Field                                                 | Default       | Effect                                                                                                                                    |
| ----------------------------------------------------- | ------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `upstream_distance`, `downstream_distance`            | 5000, 5000    | How far beyond the transcript span a variant still earns `upstream_gene_variant` or `downstream_gene_variant` (`--distance`)              |
| `reference_fasta`                                     | `None`        | An indexed FASTA from [vep-fasta](../vep-fasta/README.md) for the sequence-aware steps, HGVS and indel 3' shifting among them (`--fasta`) |
| `enable_indel_3prime_shift`, `max_indel_3prime_shift` | `false`, 2000 | Shift an insertion or deletion 3' along the transcript before classifying it, up to the given number of bases (`--shift_3prime`)          |
| `compute_hgvs`                                        | `false`       | Build HGVSc, HGVSp and `HGVS_OFFSET` (`--hgvs`; `vep-cli` also sets it for the Variant Call Format (VCF), JSON and Parquet outputs)       |
| `compute_exon_intron_numbers`                         | `false`       | Fill the `EXON` and `INTRON` ordinals (`--numbers`, and whenever a plugin runs)                                                           |
| `protein_version`                                     | `false`       | Append the translation version to the protein identifier (`--protein_version`)                                                            |
| `populate_loftee_context`                             | `false`       | Read by the caller, which attaches the transcript model to each row for the LoFTEE plugin; the crate itself leaves it empty               |

The other callables a caller reaches for, each one phrase:

- `find_overlapping_transcripts(variant, transcripts, up, down)`: the transcripts of a slice whose extended region overlaps the variant, strand-aware.
- `sv::calculate_sv_consequences(variant, transcript, up, down)`: the structural-variant dispatcher, also reachable through `calculate_consequences`.
- `sv::breakend::calculate_paired_mate(variant, transcript, config)`: the mate-side row of a paired breakend, built at the mate coordinate.
- `mapper::map_genomic_to_transcript(pos, transcript, up, down) -> Option<TranscriptPosition>`: where one genomic base falls on a transcript.
- `mapper::map_genomic_span_to_cds_bounds(start, end, transcript) -> Option<CdsSpanBounds>`: the CDS and peptide bounds of a span, intronic ends snapped to the nearest exon boundary.
- `coding::get_codon_change(variant, transcript, cds_pos) -> Option<CodonChange>`: the reference and alternate codons and residues of one variant, the changed base in upper case.
- `display::display_fields(variant, transcript) -> DisplayFields`: the five display columns on their own.
- `hgvs::generate_hgvs(variant, transcript, fasta) -> HgvsNotation`: HGVSc, HGVSp and the shift offset in one call; `generate_hgvsc` and `generate_hgvsp` give each alone.
- `consequences::shift_indel_3prime_coords(variant, transcript, fasta, max_shift) -> Option<(u64, u64)>`: the shifted span of an insertion or deletion.

`cargo doc` renders the full item list:

```bash
cargo doc --no-deps -p vep-effects
```

## Run it

The crate is a library with no binary. Build and test it alone:

```bash
cargo build -p vep-effects
cargo test -p vep-effects
```

The test target holds the module unit tests and the targeted concordance tests. The
one test marked `#[ignore]` fails by design when run; run it with:

```bash
cargo test -p vep-effects -- --ignored
```

The criterion benchmark times `calculate_consequences` on a substitution, a deletion,
an insertion and an upstream variant against a three-exon transcript, and
`find_overlapping_transcripts` over one transcript:

```bash
cargo bench -p vep-effects
```

The repository's standard commands, the full workspace test run and the targeted run
for consequence changes (`cargo test -p vep-effects consequences && cargo test -p vep-cli`)
are in [CONTRIBUTING.md](../../CONTRIBUTING.md).

## Related

- [vep-core](../vep-core/README.md): the `InputVariant`, `Transcript`, `Consequence` and `TranscriptConsequence` types this crate reads and writes, and the codon tables.
- [vep-fasta](../vep-fasta/README.md): the indexed FASTA reader behind `EffectsConfig::reference_fasta`.
- [vep-cli](../vep-cli/README.md): the binary that selects the transcripts, calls this crate once per transcript and formats the rows.
- [docs/output-formats.md](../../docs/output-formats.md): how the terms, the display columns and the `IMPACT`, `DISTANCE`, `EXON` and `INTRON` fields appear in each output format.
- [docs/intended-divergences.md](../../docs/intended-divergences.md): every case where this crate's terms differ from Ensembl VEP's, each traced to the mechanism in Ensembl VEP's source.
- [docs/plugins.md](../../docs/plugins.md): the LoFTEE consequence filter that reads the exon and intron model this crate's rows carry.
- [Root README](../../README.md): what vep-rs is, the measured concordance and the quick start.
