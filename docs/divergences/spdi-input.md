# SPDI notation resolved from the cache and the FASTA

A difference in which inputs resolve, on an input format no measured dataset uses; nothing is set aside.

Ensembl VEP reads SPDI notation (`--format spdi`) only with a database connection (`Parser/SPDI.pm` line 96, `Cannot use SPDI format in offline mode`). Each line goes to its variation database adaptor (`VariationFeatureAdaptor.pm` lines 2477-2611). The adaptor resolves the sequence name through the core database's slice adaptor (line 2533) and places the variant at the 0-based position plus one (line 2524). It then checks the deleted bases against the genome or reads them from it (lines 2545-2596).

vep-rs resolves the same notations from the cache and the reference FASTA: the sequence name through the cache's chromosomes and their synonyms, the reference bases from `--fasta`, which the format therefore needs. The shapes, positions, refusals and warning texts are the adaptor's, and the `tests/golden/116/GRCh38-input-spdi` corpus compares its rows with Ensembl VEP 116.2's database-mode output. Three readings differ.

## A version-less accession

When no sequence region or synonym matches a name exactly, Ensembl VEP's slice adaptor queries the synonym table again with the name as a prefix (the core API's `SliceAdaptor.pm` lines 2871-2884). So `NC_000021` resolves through the synonym `NC_000021.9` to chromosome 21, and the variation API's own tests (`modules/t/variationFeatureAdaptor.t`) assert that `NC_000013:32954017:G:0` is a valid deletion against a database whose synonyms carry `NC_000013.10` alone. The adaptor strips a version only from an `LRG` or `NT` name before the lookup (lines 2499-2502).

vep-rs matches a name against the cache's chromosomes and their synonyms exactly, regardless of ASCII case, never as a prefix. A version-less accession names nothing in the cache's table (`21 NC_000021.9`), so the record is skipped with the adaptor's own text for a name no slice resolves:

```
NC_000021:10538884:A:C
```

```
WARNING: line 1 skipped (NC_000021:10538884:A:C): Unable to parse SPDI
notation 'NC_000021:10538884:A:C': Sequence name NC_000021 not valid
```

Ensembl VEP annotates it as it annotates `NC_000021.9:10538884:A:C`, which both engines write as nine rows on ENSG00000274391 at `21:10538885`, Allele `C`, each `intron_variant` (two with `NMD_transcript_variant` beside it), named by the notation as written.

An HGVS genomic notation reads the same way: `fetch_by_hgvs_notation` fetches its slice through the same adaptor (`VariationFeatureAdaptor.pm` line 1889), so Ensembl VEP annotates `NC_000021:g.30282169C>A` as it annotates `NC_000021.9:g.30282169C>A`, a `stop_lost` at `21:30282169`, named as written. vep-rs matches the whole name, in any letter case, against the FASTA's records and their aliases, never a prefix, and skips the notation:

```
WARNING: line 1 skipped (NC_000021:g.30282169C>A): Unable to parse HGVS
notation 'NC_000021:g.30282169C>A': The reference FASTA has no sequence named
NC_000021
```

## A span past the end of the chromosome

Ensembl VEP fetches the whole-sequence slice before it reads the shape (line 2533) and tests no coordinate against the sequence's length. An insertion reads no sequence at all (lines 2545-2549). A deletion by length reads a slice the slice adaptor builds past the end with no bounds test (`SliceAdaptor.pm` lines 430-438), which under `--fasta` is filled with `N` (`Utils/FastaSequence.pm` lines 619-639, one `N` for a start past the end). The record checks (`Parser.pm` lines 457-649) test no sequence length either, so the record is annotated at a coordinate the chromosome lacks. Chromosome 21 of GRCh38 ends at 46709983; run on the equivalent VCF records, Ensembl VEP 116.2 writes one intergenic row for each:

```
21	46709984	past_end_ins	N	NA
21	46709983	past_end_del	NN	N
```

| Record | Location | Allele | Consequence |
| --- | --- | --- | --- |
| `past_end_ins` | `21:46709984-46709985` | `A` | `intergenic_variant` |
| `past_end_del` | `21:46709984` | `-` | `intergenic_variant` |

vep-rs reads the sequence's length before it reads the shape and refuses a span the FASTA does not hold, skipping the record with `Span <start>-<end> lies past the end of sequence <name> (<length> bases)`. So `NC_000021.9:46709984::A`, an insertion after a base the chromosome lacks, and `NC_000021.9:46709983:1:`, a deletion of the base after the last, are each skipped with that warning and write no row. An insertion after the last base, `NC_000021.9:46709983::A`, is annotated by both engines as `21:46709983-46709984 A intergenic_variant`.

### An HGVS genomic notation past the end

An HGVS genomic notation (`g.`, `m.`) whose span runs past the last base is read the same way by Ensembl VEP. `fetch_by_hgvs_notation` takes the span's bases from the slice (`VariationFeatureAdaptor.pm` lines 2060-2062), which the core sequence adaptor extends past the end with no bounds test and pads with `N` (`SequenceAdaptor.pm` lines 148-152 and 222-233). So `21:g.46709984del` is the deletion `N/-` at 46709984 and `21:g.46709983_46709984insA` the insertion `-/A` after the last base, and each is annotated as an intergenic row.

vep-rs reads the genome before it reads the shape and skips each with `Unable to parse HGVS notation '<notation>': The reference FASTA has no bases at 21:46709984-46709984` (for the insertion, `at 21:46709983-46709984`). The last base itself is inside the sequence: `21:g.46709983del` and `21:g.46709982_46709983insA` annotate as intergenic rows on both engines.

## A deletion of length 0

The adaptor reads `<sequence>:<position>:0:` as a deletion (line 2552) and ends the span at the position (line 2563). The reference allele is the empty slice's sequence (line 2564): the core API's `Slice.pm` line 693 returns an empty string for a slice whose start is its end plus one. The adaptor keeps the feature `/-` from the position plus one to the position (line 2566), a feature with no allele to annotate. So Ensembl VEP 116.2 writes no row for it in the default and tab formats, carries its VCF line with an empty REF column and no consequences, and writes a JSON object with `allele_string` `/-` and `most_severe_consequence` `?`.

vep-rs refuses the notation and skips the record, naming the empty span:

```
21:10538884:0:
```

```
WARNING: line 1 skipped (21:10538884:0:): Unable to parse SPDI notation
'21:10538884:0:': Deleted sequence length (0) leaves nothing to delete: span
10538885-10538884 is empty
```

A length written `00` is read the same way by both engines. The notation with `0` for both parts, `21:10538884:0:0`, is refused by both with `SPDI notation not supported.` (lines 2540-2542).

To reproduce: vep-rs with `--format spdi --fasta` over the GRCh38 cache and the records above, and with `--format hgvs --fasta` for the notations. Ensembl VEP 116.2 reads both formats with a database connection alone; it reads the two VCF records with the GRCh38 cache and [the command of the closing section](reproducing.md#reproducing-these-examples), adding `--fasta`.
