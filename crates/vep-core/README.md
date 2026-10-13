# vep-core

The data model of vep-rs, a Rust reimplementation of the Ensembl Variant Effect Predictor (VEP): the variant, transcript, consequence and coordinate types every other crate reads and writes, plus the shared structural-variant, codon and prediction-matrix rules. All eight sibling crates depend on it: `vep-io` builds an `InputVariant` from each input line, `vep-effects` fills it with `TranscriptConsequence` rows, and `vep-cli`, `vep-builtins` and `vep-plugin` read and extend those rows. The crate writes no files; its product is the in-memory record the writers print.

Throughout this README "Perl VEP" is the reference implementation, as in the [root README](../../README.md).

## How it works

Take the VCF (Variant Call Format) line `chr21 25000100 . ACG A . . .`, a two-base deletion. The reader in `vep-io` hands the `POS` text to `coordinate::parse_pos`, which accepts any spelling that names an exact integer, the way Perl numifies the column (`3e7` is a position, `30425949.5` is not), and the `INFO` text to `vcf_info::info_value`, which returns the last `SVTYPE` as Perl's hash would hold it. `sv_routing::written_alts` splits the `ALT` column as Perl's `split` does, and `sv_routing::routes_to_structural` answers whether the line is a structural record; here `REF` and `ALT` are plain `ACGT` runs, so the reader trims the shared `A`, computes the one-based inclusive span and calls `InputVariant::new`, which normalises `chr21` to `21` through `coordinate::normalize_chromosome`, keeps `chr21` in `original_chr` for output, builds the allele string `CG/-` and classifies the pair as `VariantClass::Deletion` through `classify_variant`. The variant then reaches `vep-effects`, which reads each overlapping `Transcript` (exons, introns, the coordinate mapper and translateable sequence of its `TranscriptVEFC`, loaded by `vep-cli` from the cache), takes the per-transcript facts from `Transcript::facts()` once per transcript, translates reference and alternate codons through `codon::translate_codon_with_table` under the transcript's codon table, and pushes one `TranscriptConsequence` per transcript into `transcript_consequences`, each carrying a `ConsequenceList` of Sequence Ontology (SO) terms whose `rank()` and `impact()` follow Ensembl Variation's table. Plugins add `plugin_data` entries to the variant or to a consequence row. The writers finish the record through `location()`, `uploaded_variation()`, `display_allele()` and `Consequence::so_term()`. Had the `ALT` been the symbolic `<DEL>`, `routes_to_structural` would have sent the line to `sv_record::structural_variants`, which gives the record one SO term and class through `sv_routing::sv_type_abbreviation` and `sv_so_term_and_class`, spans each allele from `POS` + 1 through `sv_span::structural_start` and `structural_end`, builds a breakpoint record's own breakend and its mates, and expands a `<CNV:TR>` allele's length into `tr_alt_bases`.

## Modules

| Module           | Role                                                                                                                                                                                                                  |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `lib.rs`         | Module list and the workspace constants: the Ensembl VEP release the binary reports, the default upstream and downstream distance, buffer size and cache region size                                                   |
| `allele.rs`      | Shared-prefix and suffix trimming to VEP-style minimised alleles (`"-"` for an empty side), transition versus transversion classification, IUPAC (International Union of Pure and Applied Chemistry) ambiguity codes |
| `assembly.rs`    | `Assembly`, the reference genome parsed leniently from `--assembly` (`GRCh37`, `hg19` and `37` are one value) with its Ensembl and UCSC (University of California, Santa Cruz) names                                 |
| `cache_info.rs`  | `CacheInfo`, the metadata a cache directory's `info.txt` or `info.json` carries: species, assembly, version, prediction availability, variation columns, source versions                                              |
| `codon.rs`       | Standard-code codon translation with the vertebrate mitochondrial table's four exceptions, complement and reverse complement, three-letter amino acid names                                                            |
| `consequence.rs` | The 41 SO consequence terms with rank, impact, accession and text; `TranscriptConsequence`, the row the effects engine produces; `LofteeContext`, the transcript structure a consequence-filtering plugin reads       |
| `coordinate.rs`  | `Strand`, `GenomicPosition`, `GenomicRange`, the coordinate ceiling `MAX_COORD`, Perl-faithful `POS` parsing and `normalize_chromosome`                                                                                |
| `prediction.rs`  | Decoder for the gzip-compressed SIFT (Sorting Intolerant From Tolerant) and PolyPhen prediction matrices embedded in the cache, and the Perl-style formatting of a score                                               |
| `sv_record.rs`   | `structural_variants`, which turns one structural VCF record into the variants it writes rows for: class, span, breakends and mates, tandem-repeat expansion, allele numbering                                        |
| `sv_routing.rs`  | Whether a line is structural, how its `ALT` list is read, the type abbreviation Perl VEP derives from the joined ALTs and the SO term and class that abbreviation names                                                  |
| `sv_span.rs`     | The start and end of a symbolic record from `POS`, `SVLEN` and `END` under the breakpoint, tandem and general rules, with the errors under which Perl VEP skips the record                                                   |
| `synonyms.rs`    | `ChromosomeSynonyms`, a bidirectional name map with built-in `chr` defaults and a tab-separated file loader                                                                                                            |
| `transcript.rs`  | `Transcript` and its parts: `Exon`, `Intron`, `Translation`, `TranscriptVEFC`, the cDNA-to-genome `ExonCoordMapper`, protein features, prediction matrices, and the lazily computed per-transcript facts and ordinals |
| `variant.rs`     | `InputVariant`, one allele of one input record with its annotation results; `VariantClass` and its SO term; `ColocatedVariant`; the output helpers that print a variant as Perl VEP names it                            |
| `variation.rs`   | `CachedVariation`, one known variant from the variation cache: name, alleles, minor allele frequency, clinical significance, population frequencies                                                                   |
| `vcf_info.rs`    | `INFO` column lookups that reproduce the hash Perl VEP's parser builds: the last occurrence wins, a bare key is a flag, keys match exactly                                                                               |

## Public API

The items a caller reaches for, by module. The record types (`InputVariant`, `TranscriptConsequence`, `Transcript`, `CachedVariation`, `CacheInfo`) derive `serde::Serialize` and `serde::Deserialize`, and `InputVariant` and `TranscriptConsequence` are the JSON (JavaScript Object Notation) documents a dynamic-library plugin receives ([docs/plugins.md](../../docs/plugins.md#what-run-receives)).

| Item                                                             | Module        | What it gives a caller                                                                                                                                   |
| ---------------------------------------------------------------- | ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `InputVariant`                                                   | `variant`     | One allele of one input record: coordinates, alleles, class, `transcript_consequences`, `colocated_variants`, `plugin_data`, the structural-variant fields |
| `InputVariant::new`, `InputVariant::without_allele`              | `variant`     | Build a variant from a chromosome, span and allele pair; build the no-allele record of a line whose every `ALT` is `*`                                      |
| `location`, `uploaded_variation`, `variation_name`, `display_allele` | `variant` | The Location, `Uploaded_variation`, JSON `id` and Allele columns as Perl VEP prints them                                                                   |
| `VariantClass`, `classify_variant`                               | `variant`     | The sequence and structural classes with `so_term()`, `is_structural()` and `is_span_type()`                                                              |
| `ColocatedVariant`                                               | `variant`     | A known variant matched at the input position                                                                                                              |
| `Consequence`, `ConsequenceList`, `Impact`                       | `consequence` | The 41 terms with `rank()`, `impact()`, `so_term()`, `so_accession()`, `from_so_term()` and `most_severe()`; `Ord` follows rank                          |
| `TranscriptConsequence`, `FeatureType`                           | `consequence` | One row of annotation against a transcript, regulatory feature or motif feature, with `feature_overlap()` for the structural-variant overlap columns     |
| `LofteeContext`                                                  | `consequence` | Exons, introns, coding bounds and `cds_*_NF` flags attached to a row only when a consequence-filtering plugin is active                                    |
| `Transcript`, `TranscriptVEFC`, `ExonCoordMapper`                | `transcript`  | The cached transcript model; `facts()` and `feature_ordinals()` compute once per transcript and are shared by every row against it                        |
| `PredictionMatrix::lookup`, `decode_prediction`, `format_prediction` | `transcript`, `prediction` | A SIFT or PolyPhen prediction for a protein position and amino acid, printed as Perl prints the score (`0.01`, not `0.010`)              |
| `Strand`, `GenomicRange`, `parse_pos`, `normalize_chromosome`    | `coordinate`  | Strand as `1` or `-1`, one-based inclusive ranges, `POS` as Perl numifies it, `chr21` to `21` and `chrM` to `MT`                                            |
| `translate_codon_with_table`, `reverse_complement`               | `codon`       | Amino acid for a codon under the standard or mitochondrial table; strand conversion of a sequence                                                          |
| `trim_alleles`, `classify_ts_tv`, `ambiguity_code`               | `allele`      | Minimised alleles with the bases trimmed from the start; transition or transversion; the IUPAC code for a base set                                        |
| `Assembly::parse`                                                | `assembly`    | `Grch37` or `Grch38` from any common spelling, `None` for anything else                                                                                    |
| `structural_variants`, `tandem_repeat_alt_bases`, `parse_ci_field` | `sv_record` | The variants of one structural record; a `<CNV:TR>` allele's expanded length; a `CIPOS` or `CIEND` pair                                                   |
| `routes_to_structural`, `written_alts`, `sv_type_abbreviation`, `sv_so_term_and_class`, `perl_true` | `sv_routing` | The structural-path decision, the `ALT` list as Perl reads it, the record's type and SO term, and Perl's string truth (`""` and `0` are false, every other string is true), applied to an `SVTYPE` value and to a Perl-style plugin switch such as LoFTEE's `run_splice_predictions` |
| `structural_start`, `structural_end`, `SpanClass`, `SpanError`   | `sv_span`     | A symbolic record's span and the reason Perl VEP would skip it                                                                                              |
| `info_field`, `info_value`, `info_value_written`                 | `vcf_info`    | An `INFO` key as Perl's hash holds it: flag, written value or absent                                                                                          |
| `CacheInfo`, `CachedVariation`, `ChromosomeSynonyms`             | `cache_info`, `variation`, `synonyms` | Cache metadata, a known variant from the variation cache, a chromosome name map                                                   |

The full item list with doc comments is the rustdoc, built below.

## Run it

The crate follows the workspace's standard commands ([the contributing guide](../../CONTRIBUTING.md)). Scoped to this crate:

```bash
cargo build -p vep-core
cargo test -p vep-core
cargo doc --no-deps -p vep-core
```

`cargo test -p vep-core` runs the modules' unit tests and `tests/prediction_fixtures.rs`, which decodes two prediction matrices taken from Ensembl's public variation database and checks their per-category counts and sampled cells against the packing rule of Ensembl's `ProteinFunctionPredictionMatrix.pm`.

## Related

- Depends on vep-core: [vep-io](../vep-io/README.md) (readers build `InputVariant`, writers print it), [vep-effects](../vep-effects/README.md) (fills `transcript_consequences` against `Transcript`), [vep-cli](../vep-cli/README.md) (loads the cache into `Transcript`, `CachedVariation` and `CacheInfo`, drives the pipeline), [vep-builtins](../vep-builtins/README.md) and [vep-plugin](../vep-plugin/README.md) (read and extend the rows), [vep-fasta](../vep-fasta/README.md), [vep-cache-builder](../vep-cache-builder/README.md) and [vep-cache-converter](../vep-cache-converter/README.md).
- vep-core depends on no other workspace crate.
- [docs/plugins.md](../../docs/plugins.md): the plugin interface, whose JSON documents are the serde form of `InputVariant` and `TranscriptConsequence`.
- [docs/cache-setup.md](../../docs/cache-setup.md#sift-and-polyphen-predictions): where the prediction matrices come from and how the decoder prints them.
- [docs/intended-divergences.md](../../docs/intended-divergences.md): the records Perl VEP drops before annotation, which `sv_routing` and `sv_record` carry under the `vep_skip` mark, and the symbolic `<CNV:TR>` reading whose expanded length `tr_alt_bases` holds.
- [docs/output-formats.md](../../docs/output-formats.md): the columns the output helpers of `InputVariant` and `Consequence` print.
- [Root README](../../README.md): what vep-rs is, the plugin matrix and the documentation index.
