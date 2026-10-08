# Changelog

All notable changes to this project are recorded in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `--protein_version` appends the translation version to the `ENSP` identifier (`ENSP00000286808.3`) when `--protein` prints it, as Ensembl VEP 116 does. (#N)
- `--mask_header_cache_path` prints the cache directory as `[PATH]/<leaf>` in the default and `--tab` headers, the VCF `cache=` token and the Parquet footer. (#N)
- `--regulatory_gff`, `--extended_promoters` and `--custom_suppress_filter` are refused with an error naming what each would do; regulatory and custom annotation are features of a later release. (#N)
- `--reference-release {115.2,116.2}` on both comparators selects the Ensembl VEP release scored against (116.2 default, sets aside the records it skips); the runners derive it from the reference's provenance. (#N)

### Changed

- Ensembl VEP 116.2 is the reference release: the output headers and the banner read `v116.2`, and `--cache_version` defaults to 116 (it names only the header's cache path). (#N)
- The golden corpora under `tests/golden/116/` carry Ensembl VEP 116.2's output and replace the release 115 corpora under `tests/golden/115/`. (#N)
- Stop-codon consequences follow Ensembl VEP release 116: an edit beginning inside the stop codon is `stop_lost` or `stop_retained_variant` without `frameshift_variant`; indels decide `stop_retained_variant` on the genomic span. (#N)
- A breakend record reaching no transcript at either breakend writes one `intergenic_variant` row per allele, and a bracket record's own breakend prints `N.` without a MATEID (Ensembl VEP 116). (#N)
- A multi-allelic structural variant takes one class from its joined ALTs (`<INS>,<INS>` an insertion, any DEL beside DUP a `copy_number_variation`) and writes one row per ALT (Ensembl VEP 116). (#N)
- `--max_sv_size -1` lifts the size limit, so a structural variant of any span is annotated in every output format. (#N)
- Records Ensembl VEP 116 drops before annotation (oversize, an unsupported type or ALT list, a `<DEL>` without `END` or `SVLEN`) are annotated, a documented divergence the golden manifests class `reference_skipped_record`. (#N)
- JSON `cdna_end`, `cds_end` and `protein_end` are omitted for a position whose end is undefined (`445-?` writes `cdna_start` alone), as Ensembl VEP 116 writes them. (#N)
- `Uploaded_variation` of a record whose ID is `.` is built from the input line, `CHROM_POS_REF/ALT1/ALT2` as written, as Ensembl VEP 116 names it (`21_43512967_AT/ATT/A`, `21_33867341_C/<CN2>`). (#N)
- `--hgvs` describes a deletion that runs over either end of a transcript over the bases inside it (`c.-49_77+1248del`), as Ensembl VEP 116 does, instead of omitting HGVSc. (#N)
- `--hgvs` writes a frameshift whose first changed residue is the stop codon as an extension, `p.Ter124IleextTer14` (count excluding the replaced stop), as Ensembl VEP 116 does, not `p.Ter124IlefsTer15`. (#N)
- `docs/intended-divergences.md` is rewritten for readability and organised by
  the five Ensembl VEP defects behind the classes the adjusted concordance sets
  aside. Each defect's section states the defect in one sentence, why it is a
  defect, what vep-rs does instead, where it is in the source and what
  correcting it would change; each class opens with a summary and walks through
  its record in short sections, with citations and counts as lists. Every
  record, row, number, citation and command is kept. The `<CNV:TR>`
  representation difference has a section of its own, the former class 9
  (classes 1 and 2 on the structural-variant sets) is split into one class per
  defect, and the classes are renumbered in page order: the former 1, 9 (splice
  half), 2, 9 (start half), 3, 7, 8, 10, 4, 6, 11 and 5 are now 1 to 12. The
  0.3.1 entry's "section 3" is section 5 of the reorganised page.

### Fixed

- `--protein` fills ENSP on every structural-variant transcript row (`transcript_ablation`, `feature_truncation`, breakends) with a translation, where it printed `-`. (#N)
- A breakend whose mate lies inside an NMD transcript while its own position lies outside carries `feature_truncation` without `NMD_transcript_variant`. (#N)
- `run_clone_measurement.sh --engine perl` structural-variant cells run to completion: the cache version reaches the per-VCF child as an argument, not as a variable unbound under `set -u`. (#N)

## [0.3.2] - 2026-10-06

### Changed

- `--fork 0` (the default) uses every logical CPU up to 32 instead of every
  logical CPU: one reader feeds the worker pool, so the wall time stops falling
  past about 32 workers and rises on wider hosts (gnomAD v4.1 chromosome 2 took
  1.8 times longer at 192 workers than at 32 on a 192-vCPU host). An explicit
  `--fork N` is never capped.
- A batch is annotated, filtered, rendered and freed in one pass over the worker
  pool, each finished variant freed by the thread that annotated it and the
  task size following the pool's thread count, instead of five pool barriers per
  batch with the coordinator freeing every variant serially between batches. The
  output is byte-identical. Together with the libdeflate inflate below, measured
  against 0.3.1 under the same whole-genome protocol (the median of 20 machines
  per chromosome at 16 threads, summed over the chromosomes), the whole of gnomAD
  v4.1 falls from 974 to 724 s on ARM Graviton4 and from 1,219 to 1,169 s on x86
  and the whole of gnomAD v2.1.1 from 586 to 282 s and from 835 to 476 s, with
  CPU time down 19 to 44 percent on those inputs and peak memory at most 2.4 GB
  on any chromosome against 4.9 GB.
- bgzip blocks inflate through libdeflate (`noodles-bgzf` with its `libdeflate`
  feature), and the cache builder's MySQL client takes its Rust deflate backend
  so that no crate of the workspace links C zlib: the single-inflater floor on
  gnomAD v4.1 chromosome 21 falls from 73.8 s to 34.5 s.

## [0.3.1] - 2026-10-02

### Added

- The SNP/indel comparator applies the transcript-selection class of
  `docs/intended-divergences.md` (section 3) to the records of a whole-genome
  input that Ensembl VEP marks skipped, under the structural-variant
  comparator's two guards, and reports it as `sv_transcript_selection`; the
  scope test is defined once and shared by the two comparators.

### Fixed

Each entry names the Ensembl VEP mechanism vep-rs follows.

- A symbolic insertion whose `END` is its position or whose `SVLEN` is 0 is the
  insertion between the two bases, `start = POS + 1, end = POS`, as ensembl-io's
  `get_end` reads it: its Location names both flanks, its predicates read the
  pair as `_intron_effects` reads a sequence insertion (both flanks inside a
  region, the two insertion special cases at the intron edges, the UTR special
  cases at the coding region's edges, `within_cdna` through `map_insert`), and
  the row carries no `OverlapBP` or `OverlapPC` (an overlap length of 0). A
  symbolic record with neither field stays the single base at `POS + 1`, which
  is also ensembl-io's reading.
- UTR term on an insertion at the edge of the coding region: `_before_coding`
  and `_after_coding` special-case an insertion whose start is the coding
  region start or whose end is the coding region end, so an insertion between
  the last base of an intron and the first base of a coding exon (or between
  the coding region's last base and the intron after it) is `5_prime_UTR_variant` or
  `3_prime_UTR_variant` by strand beside its splice term, where it carried the
  splice term alone.
- A frameshift intron on the UTR side of the coding region: `within_cds`
  requires the span to overlap the coding region, and `_bvfo_preds` sets `utr`
  rather than `coding` for a span that does not, so a variant inside an intron
  of twelve bases or fewer that lies before the coding region is
  `5_prime_UTR_variant` (or `3_prime_UTR_variant` after it) rather than
  `coding_sequence_variant`.
- Structural `stop_lost` behind the `coding` pre-predicate: the term carries
  `include {coding => 1}`, and `_bvfo_preds` sets `coding` only for a span that
  overlaps the coding region and an exon, so a deletion inside the intron that
  splits a stop codon is `non_coding` and carries no `stop_lost`, although the
  codon window drawn on genomic coordinates reaches into that intron.
- No UTR term on a structural span that starts on the first base of a
  transcript whose coding region starts there (`cds_start_NF`, and likewise at
  the last base with `cds_end_NF`): `_before_coding` is
  `overlap(start, end, transcript start, coding_region_start - 1)`, an inverted
  window that admits only a span starting before the transcript.

## [0.3.0] - 2026-10-01

### Added

- A container image, `ghcr.io/natera-open-source/vep-rs`, for each published
  release from 0.2.0: the release's attested x86_64 Linux binaries and the
  `duckdb` CLI on Debian 13, smoke-tested on the GRCh37 corpus before it is
  pushed, with a build-provenance attestation on the pushed digest
  (`.github/workflows/docker.yml`).
- `docs/intended-divergences.md`: every class the adjusted concordance sets
  aside, each with the record, both engines' rows, the mechanism in Ensembl
  VEP's source, and a command that reproduces it.
- Structural-variant comparator masks for the second side of two classes and
  for three more: the `intergenic_variant` row vep-rs writes for a record
  whose VEP tuples all name transcripts on another chromosome
  (`filter_cross_chromosome_orphan_intergenic`); the mate allele of a
  breakend above `--max_sv_size` that VEP never wrote, gated on the mate
  lying within 5 kb of the transcript span
  (`filter_giant_breakend_mate_divergences`); the SNP/indel registry's
  excluding rules applied to structural pairs (`filter_registry_swap_pairs`);
  VEP's batch-dependent Transcript tuples on a `<NON_REF>` record
  (`filter_non_ref_batch_divergences`); and a mate-side row VEP evaluated at
  the local coordinate (`filter_breakend_mate_local_read_pairs`), with
  `scripts/validation/breakend_mate_context.py` deriving the mate-coordinate
  and local-coordinate readings from the JSON cache. The transcript-selection
  mask's unsupported-type scope names `<CPX>` and `<CTX>`.
- The golden classifier names `breakend_mate_local_read`, computed by the
  same derivation over the corpus's own pruned cache.
- The comparators write the adjusted residual: `discordant_open.tsv` beside
  `discordant.tsv` for the SNP and indel sets and `discordant_adjusted.tsv`
  for the structural-variant sets, each the one-sided rows every mask leaves
  in place. Per-class mode caches the per-term VEP totals of a reference
  output (`--vep-totals-cache`).

### Changed

- JSON cache format: each transcript's `variation_effect_feature_cache`
  carries `seq_edits`, the translation's `Bio::EnsEMBL::SeqEdit` list
  (`start`, `end` and `alt_seq` in protein coordinates, with `code`, `name`
  and `length_diff`), which `scripts/data/storable_to_json.pl` writes from
  the Storable cache. vep-rs applies the edits overlapping a variant's
  translation span to the reference allele's peptide, as
  `TranscriptVariationAllele::peptide` does. A cache converted without the
  key still loads: no edit is applied, and a variant at codon 1 of a
  transcript with an alternative initiation codon is called from the cached
  peptide. Convert such a cache again to get VEP's calls there.
- The deletion-shaped evaluation of a same-chromosome breakend span borrows
  the variant for each overlapping transcript instead of copying it with the
  consequences accumulated so far, which grew with the square of the
  transcript count. Output is unchanged.

### Fixed

Each entry names the Ensembl VEP mechanism vep-rs follows.

- Start codon on an alternative initiation codon: with `seq_edits` in the
  cache, an `initial_met` or `amino_acid_sub` edit over a non-ATG start
  codon gives a reference `M` against the alternate allele's literal
  residue, so the variant is `start_lost` rather than `synonymous_variant`
  and `Amino_acids` prints the edited residue.
- Splice terms across a frameshift intron: `_intron_effects` skips an intron
  of twelve bases or fewer for every differing region that overlaps it, in
  both of its loops, so a deletion running from one exon across such an
  intron into the next carries no splice term; the coding predicates alone
  decide the row. The exonic splice windows apply the same skip, so a donor
  window that crosses such an intron yields no `splice_donor_5th_base_variant`.
- `transcript_ablation` requires Perl's `deletion` pre-predicate (the
  reference span longer than the ALT sequence) beside `complete_overlap`. A
  same-length substitution or a longer ALT covering a whole transcript takes
  the tier-3 predicates instead: the biotype's context term with the UTR,
  intron and splice terms of the span.
- `splice_polypyrimidine_tract_variant` on structural alleles: a structural
  allele whose span reaches the 15 bp window inside an intron's acceptor end
  carries the term unless the span touches an exon of the transcript, as
  Ensembl's predicate does for every BaseVariationFeature. A `<CNV:TR>`
  allele keeps its endpoint test.
- A paired breakend's bracket allele gets a row for a transcript only when
  the mate coordinate lies within `MAX_DISTANCE_FROM_TRANSCRIPT` (5 kb) of
  it (`StructuralVariationOverlap::_close_to_feature`); the local POS is
  carried by the single-breakend allele.
- A single-breakend allele above `--max_sv_size` that partially covers a
  coding transcript takes the region predicates like every other
  same-chromosome span; VEP requires `complete_overlap_feature` for
  `coding_transcript_variant` at any span.
- The context term of a transcript a structural allele engulfs follows the
  biotype, as VEP's include hashes do: `coding_transcript_variant` only for
  `protein_coding`, `NMD_transcript_variant` only for
  `nonsense_mediated_decay`, `non_coding_transcript_variant` only without a
  translation, and `intergenic_variant` for a translated transcript of any
  other biotype. One rule serves the inversion, copy-number, insertion and
  breakend arms.
- Structural insertions, mobile-element insertions and duplications whose
  span covers a start-codon base with both ends in exons emit `start_lost`,
  the genomic-overlap predicate the deletion arm applies.
  `start_retained_variant`, which VEP co-emits there without reading
  sequence, is not emitted.
- The mate-side Transcript row of a paired breakend describes the mate
  breakend at the mate coordinate: `feature_truncation` with the region term
  when the mate lies inside the transcript, `upstream_gene_variant` or
  `downstream_gene_variant` with the distance when it lies within the window
  outside, and no row beyond it. VEP evaluates every predicate but
  `feature_truncation` at the local coordinate there.
- A gVCF `<NON_REF>` record is a reference-confidence block, not an
  alternate allele. It is written as one `intergenic_variant` row spanning
  POS+1 to END with the allele as written and `IMPACT=MODIFIER`, with a
  warning per record; that is the row VEP writes when no other record in the
  batch loads the region. A `<NON_REF>` carrying a supported `INFO/SVTYPE`
  keeps that type.

## [0.2.0] - 2026-09-29

### Fixed

- HGVSp is Perl VEP's `hgvs_protein`, ported whole. In-frame insertions
  between codons, duplications, stop-loss extensions (`extTer{n}` and
  `extTer?`), insertions carrying a stop codon and deletions shifted at the
  peptide level print VEP's string where they printed nothing or a different
  one. The coding gate is read from the allele as annotated and the
  frameshift test from its 3'-shifted span, as VEP reads them, and the
  transcript peptide as `Transcript::translate` writes it (`Met1` for any
  start codon of the table). The one intended difference: on a variant that
  keeps the start codon intact and is reported as `start_retained_variant`,
  VEP prints `p.Met1?` from its own `start_lost` call while vep-rs prints
  the peptide change.
- HGVSc is Perl VEP's `hgvs_transcript`, ported whole: the 3' shift within
  1,000 bases of flank, the variant typing and duplication lookup of
  `hgvs_variant_notation`, allele clipping and `_get_cDNA_position`.
  Deletions that cross an exon boundary or the start or stop codon,
  reverse-strand `delins` ranges, insertions without a readable reference
  and variants outside the transcript span now match VEP.
- `HGVS_OFFSET` is emitted beside `HGVSc` and `HGVSp` when an insertion or
  deletion was shifted, signed by the transcript strand as VEP signs it.

### Added

- A GRCh37 chromosome 21 golden corpus generated with `--hgvs`
  (`tests/golden/115/GRCh37-hgvs/`), carrying the complete chromosome 21
  reference so the golden tests compare `HGVSc`, `HGVSp` and `HGVS_OFFSET`
  against Ensembl VEP's output in every format with no external file.
- Release measurements beside the paper's: the README's "Current release"
  section states the released version's ten-suite concordance and its wall
  times (20 independent machines per cell on both architectures) beneath the
  paper's tables, which stay as published, from the release's provenance
  record under `docs/concordance-provenance/`.

### Changed

- Releases: the workflow creates a draft that a maintainer publishes, with a
  body rendered from the repository's files (`scripts/release/render_release_notes.py`);
  the archives are `vep-X.Y.Z-<target>.tar.gz`, each unpacking to its own
  directory; a `git archive` source tarball ships with its `.sha256`; one
  `SHA256SUMS` covers every archive and every file in it carries a build
  provenance attestation (`gh attestation verify`); `MD5SUMS` is not published.
  `workflow_dispatch` runs the same build as a dry run without a release.
- Dependencies: noodles 0.109 (the last release declaring rust-version 1.88;
  the bgzf readers and writer and the tabix index loader moved to their
  `io` and `fs` modules), mysql 28 (vep-cache-builder), rand 0.10
  (vep-cache-converter), criterion 0.8 (benches), rustls 0.23.45, and the
  minor and patch releases of bigtools, bstr, clap, crossbeam-channel,
  flate2, indexmap, rayon, rustc-hash, serde, serde_json, smallvec, tempfile,
  thiserror and tracing-subscriber. The declared rust-version stays 1.88.
- Workflows: every action is pinned to a commit sha with its tag beside it,
  each rust-toolchain step names its Rust version through the `toolchain:`
  input, and the CI workflow's token is read-only. Dependabot ignores
  noodles 0.110 and later (rust-version 1.89) and the pinned toolchain refs.
- Golden corpus provenance records the digest of the released 0.1.0 binary
  that classified the manifests.
- The measurement harness (`scripts/concordance/run_clone_measurement.sh`)
  flushes and deletes the discarded warmup's output before every timed run,
  so the timed run starts with the memory the warmup had; the warmup's
  pages otherwise had to be reclaimed while the timed run wrote its own,
  and that reclaim landed in the timed wall time while user and system time
  stayed unchanged.

### Security

- Cleared from the dependency tree: RUSTSEC-2026-0235 (rkyv),
  RUSTSEC-2026-0173 (proc-macro-error2), RUSTSEC-2026-0002 (lru 0.12) and
  RUSTSEC-2026-0097 (rand 0.8), all through the mysql and rand updates, and
  RUSTSEC-2026-0285 (rustls) through the rustls update. None was reachable
  from the `vep` binary; all sat under vep-cache-builder or
  vep-cache-converter.

## [0.1.0] - 2026-09-22

First public release.

### Added

- `vep`, a from-scratch Rust implementation of the Ensembl Variant Effect
  Predictor algorithm. It reads a JSON transcript cache converted from the
  Ensembl VEP cache, accepts VCF, Ensembl default, HGVS and region input, and
  annotates SNVs, indels, MNPs and structural variants (deletions,
  duplications, insertions, inversions, copy-number variants and breakends)
  with Sequence Ontology consequence terms, HGVS notation, protein predictions
  and co-located variation.
- Output in the four Perl VEP formats (VEP default, VCF, JSON, Tab), each
  field for field Ensembl VEP 115.2's for the same input and flags, for every
  row outside the keys the golden manifests list as divergent (header lines carrying a time,
  a path or Perl API component versions excepted), plus Parquet, a vep-rs
  extension: five key columns (`chrom`, `pos`, `end`, `ref`, `alt`) then
  exactly the tab columns, `-` stored as NULL, Parquet v2, ZSTD 9, sorted row
  groups, a dictionary and a Bloom filter on every column, and the run's
  identity in the footer metadata.
- Golden corpora under `tests/golden/<release>/<assembly>/`: per Ensembl
  release and assembly, exemplar records for every consequence-set combination
  observed on the benchmark datasets, Ensembl VEP's output for them in the
  default, tab, VCF and JSON formats, a pruned JSON cache, and a manifest of
  the consequence keys documented to differ. `crates/vep-cli/tests/golden.rs`
  and `format_parity.rs` compare every column of every format against them on
  each `cargo test`, and the Parquet output round-trips through
  `scripts/adapters/parquet_to_vep_tab.py` to the tab output. Generator and
  pruner under `scripts/golden/`.
- `scripts/concordance/compare_vep_fields.py`: per-column and per-Extra-key
  agreement between two default-format outputs over their shared consequence
  keys, written as `fields.json` beside every F1 report the measurement
  harness produces.
- `--max_sv_size` (default 10,000,000): as in VEP, a wider structural variant
  keeps its VCF line without consequences and is absent from the JSON output;
  the default and tab formats list every transcript it overlaps.
- Parallel annotation with `--fork N`; output is deterministic and identical at
  every `--fork` value.
- Built-in plugins for CADD, REVEL, gnomADc, AlphaMissense, dbscSNV, LoFTEE,
  LoFtool and pLI on both GRCh37 and GRCh38, GWAS and SpliceAI on GRCh38, and
  dbNSFP, built in but unsupported;
  a binary annotation-store format for plugin data; and a dynamic-library
  plugin ABI (`vep-plugin`) for plugins built outside this repository.
- `vep-cache-builder`, which builds the JSON transcript cache directly from
  Ensembl source data (GFF3, protein FASTA, variation VCF, and with `--gtf`
  the release's GTF for the `cds_start_NF` / `cds_end_NF` transcript
  attributes and the assembly and genebuild versions) without a Perl
  installation, and `vep-cache-converter`, which validates an existing JSON
  cache for runtime use and converts tabix-indexed plugin data to the binary
  annotation-store format. Every JSON cache carries its `source_versions` in
  `info.json`.
- A concordance harness that compares vep-rs output against Perl VEP
  (release 115) tuple by tuple and classifies every divergence, and
  structural-variant validation VCFs for both assemblies under
  `tests/sv_validation/`.
- Documentation under `docs/` covering the CLI, output formats, plugins and
  cache setup, the runbook for reproducing the published concordance under
  `scripts/`, and the data and figure sources for the accompanying manuscript
  under `manuscript/`.
- `scripts/check_advisory_reachability.sh`, which fails when a known
  third-party advisory is reachable from the `vep` binary rather than merely
  present somewhere in the workspace, with an expiring allowlist for accepted
  risks.
- Continuous integration on GitHub Actions: tests, Clippy, formatting, a
  build on the declared minimum supported Rust version, a module-documentation
  check, a Developer Certificate of Origin check on pull
  requests, and a release workflow that publishes Linux (x86_64, aarch64) and
  macOS (aarch64) tarballs with sha256 and md5 checksum files.

[Unreleased]: https://github.com/natera-open-source/vep-rs/compare/v0.3.2...HEAD
[0.3.2]: https://github.com/natera-open-source/vep-rs/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/natera-open-source/vep-rs/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/natera-open-source/vep-rs/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/natera-open-source/vep-rs/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/natera-open-source/vep-rs/releases/tag/v0.1.0
