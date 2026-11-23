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
translates both codons with the transcript's codon table (the genetic code its cache record
names, any table BioPerl's `Bio::Tools::CodonTable` ships or NCBI's code 26: table 1 for a nuclear transcript,
2 for a vertebrate's mitochondrial one, 5 for an invertebrate's, 3 for a yeast's, 11 for a
plastid's) and compares the residues, reading the table's own start and stop codons: the
same residue yields `synonymous_variant`, a different one `missense_variant`, a new stop
`stop_gained`, a change in the start codon `start_lost`. Whole-span passes then add what
a single endpoint cannot see (a splice-region window in a neighbouring intron, an
untranslated region (UTR) the span runs into) and the suppression rules remove the pairs
Ensembl VEP never writes together. The list is sorted by severity rank, the impact is the
most severe term's, and `hgvs::generate_hgvs` adds HGVSc and HGVSp when `compute_hgvs`
is set. Back in `calculate_consequences`, `display::display_fields` renders the
five display columns (`cDNA_position`, `CDS_position`, `Protein_position`,
`Amino_acids`, `Codons`) and, under `compute_exon_intron_numbers`,
`display::exon_intron_numbers` adds the `EXON` and `INTRON` ordinals. The finished
`TranscriptConsequence` returns to `vep-cli`, which attaches it to the variant for the
plugin phase and the output formatters.

## Modules

| Module                     | Role                                                                                                                                                                                                                                                                  |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/lib.rs`               | Declares the modules and re-exports `calculate_consequences`, `shift_indel_3prime`, `shift_motif`, `EffectsConfig`, `LofteeContextRows` and `find_overlapping_transcripts`; holds the test-only transcript builders                                                  |
| `src/consequences.rs`      | `EffectsConfig` and `calculate_consequences`: classifies each end of a sequence variant, applies the position terms, the whole-span splice and UTR passes and the suppression rules, sorts the terms by severity and sets the impact                                  |
| `src/mapper.rs`            | Maps a genomic position to a `TranscriptPosition` (upstream, downstream, UTR, coding, intron, non-coding exon, or the unmapped position of a base inside the transcript span that no exon pair and no intron holds, such as one a RefSeq model's alignment deleted) and a genomic span to CDS and peptide bounds or to a segment projection |
| `src/coding.rs`            | Codon extraction and translation, the frameshift test on CDS-mapped lengths, inframe indel classification by peptide comparison, full-CDS peptide alleles, the ported Ensembl VEP coding predicates (`perl_coding_terms`, `perl_hgvs_protein`) and `rotate_shifted_allele`, which rotates a 3'-shifted allele as read on the transcript's strand, by the shift modulo its length, for the terms, codons, `Allele`, HGVSp and `USED_REF` of a shifted row |
| `src/display.rs`           | The five display columns and the `EXON` and `INTRON` ordinals, rendered through the transcript mapper as Ensembl VEP renders them, apart from the term analysis                                                                                                       |
| `src/hgvs.rs`              | HGVSc and HGVSp strings and the `HGVS_OFFSET` of an insertion or deletion shifted 3' along the transcript                                                                                                                                                             |
| `src/hgvsg.rs`             | `HGVSg`, the genomic notation of a sequence variant on the reference FASTA's chromosome slice (`--hgvsg`), and `nc_accession` and `refseq_name`, the RefSeq names among a sequence's synonyms that name it                                                                        |
| `src/spdi.rs`              | `SPDI` and the `GA4GH_VRS` allele built from it (`--spdi`, `--ga4gh_vrs`) for each allele of a record, run after the transcript rows are known: an allele with no transcript row, one with regulatory or motif rows alone included, names its chromosome as the reference FASTA does, and as the cache does without one |
| `src/regulatory.rs`        | `regulatory_terms`, the terms of a variant against a regulatory or motif feature: the overlap term, plus the ablation term for a deletion covering the whole feature or the amplification term for a copy-number gain covering it; `motif_fields`, the `MotifFeature` row's motif position, informative-position flag and binding-score change |
| `src/overlap.rs`           | `find_overlapping_transcripts`: filters a transcript slice by chromosome and by the strand-aware extended region                                                                                                                                                      |
| `src/sv/mod.rs`            | `calculate_sv_consequences`: dispatches a structural variant on its class, reads the `start = end + 1` insertion pair and the symbolic `<CNV:TR>` gain, which names a truncation for an allele that states a loss and no length term for one marked exactly as long as its run, and holds the shared overlap predicates (exon, CDS, UTR, intron, polypyrimidine tract, mature microRNA) and the shared coding guards `perl_coding_pre_predicate` and `perl_partial_codon` |
| `src/sv/deletion.rs`       | `<DEL>` and mobile-element deletions: `transcript_ablation` when the deletion contains the transcript, otherwise `feature_truncation` with the region terms; `stop_lost` over a stop-codon base unless the span's first base in transcript order lies in an incomplete terminal codon, `start_lost` over a start-codon base when both ends are exonic |
| `src/sv/duplication.rs`    | `<DUP>` and `<DUP:TANDEM>`: `transcript_amplification` when the duplication contains the transcript, `feature_elongation` when it lies inside the transcript and overlaps exonic sequence, the region terms otherwise                                                 |
| `src/sv/insertion.rs`      | `<INS>`, `<INS:ME:*>` and the `<CNV:TR>` gain: ranged, point and pair forms; `feature_elongation` inside the transcript over exonic sequence, `start_lost` over a start-codon base                                                                                    |
| `src/sv/inversion.rs`      | `<INV>`: region terms only, never `transcript_ablation`; `coding_transcript_variant` for a contained protein-coding transcript                                                                                                                                        |
| `src/sv/cnv.rs`            | `<CNV>`, `<CN0>` to `<CN9>` and `<CPX>`: `<CN0>` takes the deletion rules, `<CN2>` the duplication rules, every other copy number the overlap terms alone                                                                                                             |
| `src/sv/breakend.rs`       | Breakends in bracket and dot notation: `calculate` for the record's own breakpoint, `calculate_paired_mate` for the row on the mate chromosome; `feature_truncation` only when the allele's own breakend lies inside the transcript body and names a mate             |
| `src/concordance_tests.rs` | Test-only module of targeted tests, each asserting the Ensembl VEP-correct terms for one divergence shape; the one `#[ignore]` test pins an Ensembl VEP behaviour the helper `mapper::map_genomic_span_to_cds_bounds` does not reproduce                              |
| `src/stop_codon_tests.rs`  | Test-only module of HGVSp tests at the stop codon, each laying a transcript model at its own coordinates over a reference of synthetic bases and asserting the protein notation |
| `src/intron_end_insertion_tests.rs` | Test-only module for an insertion written at either end of an intron, each test asserting the row's terms, HGVSc, positions and HGVSp; the files under `src/intron_end_insertion_tests/` each hold one form of it |

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
| `enable_indel_3prime_shift`                           | `false`       | Shift an insertion or deletion 3' along the transcript before classifying it, over Ensembl VEP's 1,000-base search area (`--shift_3prime 1`) |
| `hgvs`                                                | `HgvsOptions::default()` | The HGVSc and HGVSp flags: `shift_hgvs` (`true`; `--shift_hgvs 0` describes an insertion or deletion where the input placed it), `hgvsp_use_prediction` (`--hgvsp_use_prediction`, `p.(Arg97Ter)`), `remove_hgvsp_version` (`--remove_hgvsp_version`, no translation version on the HGVSp name) and `shift_3prime` (`false`; `--shift_3prime 1`, set beside `enable_indel_3prime_shift`: HGVSc is written when the span the shift annotated the row at overlaps the transcript, so an insertion or deletion the shift carries inside is described at its shifted span, and HGVSp reads the stop, start and frameshift predicates on the shifted allele) |
| `compute_hgvs`                                        | `false`       | Build HGVSc, HGVSp and `HGVS_OFFSET` (`--hgvs`, with a reference FASTA present, whatever the output format)                              |
| `compute_exon_intron_numbers`                         | `false`       | Fill the `EXON` and `INTRON` ordinals (`--numbers`, and whenever a plugin runs)                                                           |
| `protein_version`                                     | `false`       | Append the translation version to the protein identifier (`--protein_version`)                                                            |
| `loftee_context`                                      | `LofteeContextRows::None` | Read by the caller, which attaches the transcript model to the rows it names (every row for LoFTEE, MaxEntScan, NMD and NearestExonJB, the 5' UTR rows for UTRAnnotator); the crate itself leaves it empty |
| `context_sequences`                                   | `LofteeContextRows::None` | Read by the caller, which attaches the transcript's coding, 3' UTR and protein sequences to the rows it names (the `frameshift_variant` rows for Downstream); the crate itself leaves it empty |
| `hgvsg`                                               | `None`        | Build `HGVSg` for every sequence variant, the chromosome named by its RefSeq accession under `--hgvsg_use_accession` (`--hgvsg`); needs `reference_fasta`               |
| `spdi`                                                | `None`        | Build `SPDI`, and the `GA4GH_VRS` allele when `vrs_accessions` is set, for every sequence variant once its transcript rows are known (`--spdi`, `--ga4gh_vrs`)           |
| `domains`                                             | `false`       | List the protein features the variant's residues overlap on each coding row (`--domains`)                                                                                 |
| `use_transcript_ref`                                  | `false`       | Analyse each row of a sequence variant against the reference allele read from the transcript and carry it as `used_ref` (`--use_transcript_ref`)                        |
| `process_ref_homs`                                    | `false`       | The record's reference row is written and reads its codon before any alternate allele's, so the `used_ref` codon path rotates nothing under a 3' shift (`--process_ref_homs`) |

The other callables a caller reaches for, each one phrase:

- `find_overlapping_transcripts(variant, transcripts, up, down)`: the transcripts of a slice whose extended region overlaps the variant, strand-aware.
- `sv::calculate_sv_consequences(variant, transcript, up, down)`: the structural-variant dispatcher, also reachable through `calculate_consequences`.
- `sv::breakend::calculate_paired_mate(variant, transcript, config)`: the mate-side row of a paired breakend, built at the mate coordinate.
- `mapper::map_genomic_to_transcript(pos, transcript, up, down) -> Option<TranscriptPosition>`: where one genomic base falls on a transcript.
- `mapper::map_genomic_span_to_cds_bounds(start, end, transcript) -> Option<CdsSpanBounds>`: the CDS and peptide bounds of a span, intronic ends snapped to the nearest exon boundary.
- `coding::get_codon_change(variant, transcript, cds_pos) -> Option<CodonChange>`: the reference and alternate codons and residues of one variant, the changed base in upper case.
- `display::display_fields(variant, transcript) -> DisplayFields`: the five display columns on their own.
- `hgvs::generate_hgvs(variant, transcript, fasta, options) -> HgvsNotation`: HGVSc, HGVSp and the shift offset in one call under an `HgvsOptions`. The notations stay where the input placed the allele when its `shift_hgvs` is false; HGVSp takes the predicted form or drops its version under `hgvsp_use_prediction` and `remove_hgvsp_version`; `shift_3prime` reads the shifted span as the `hgvs` row above says. `generate_hgvsc` and `generate_hgvsp` give each alone under the defaults.
- `consequences::shift_indel_3prime_coords(variant, transcript, fasta) -> Option<(u64, u64)>`: the span an insertion or deletion reaches after the 3' shift along the transcript's strand, `None` when nothing moves.

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

The test target holds the module unit tests, the targeted concordance tests and the
HGVSp tests at the stop codon and at an intron's ends. The one test marked `#[ignore]`
fails by design when run; run it with:

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
- [docs/plugins.md](../../docs/plugins.md): the LoFTEE consequence filter, which reads the exon and intron model `vep-cli` attaches to this crate's rows.
- [Root README](../../README.md): what vep-rs is, the measured concordance and the quick start.
