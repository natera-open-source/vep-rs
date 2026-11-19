# Cache conversion

An Ensembl VEP Storable cache converts to the JSON layout vep-rs annotates from; [Cache setup](cache-setup.md) covers building one natively with `vep-cache-builder`.

The converter is `scripts/data/storable_to_json.pl`. It runs inside the `ensemblorg/ensembl-vep` container ([scripts/README.md](../scripts/README.md) gives the invocation), and its output goes straight to `--json_cache`.

The converter reads an Ensembl VEP Storable cache directory and writes the same JSON layout for every chromosome directory the cache holds. Those are the assembled chromosomes and the patch, alternate-locus, unplaced-scaffold and `LRG_` directories alike, as Ensembl VEP itself lists them. Beside the transcripts it writes an `info.json` derived from the cache's `info.txt` and copies the cache's `chr_synonyms.txt`. The `info.json` holds `contig_order`: the converted contigs as autosomes by number, then X, Y and MT, then every other contig by name.

Every field the cache dumper keeps on a transcript is carried, flattened or dropped under a rule the parity audit names (`scripts/concordance/cache_parity_audit.pl`):

- Carried:
  - the gene, transcript and protein identifiers (UniProt isoforms included, and on a merged cache the source)
  - the attributes
  - the exons and introns
  - the translation with its sequence edits
  - the translateable sequence and peptide
  - the codon table
  - the exon coordinate mapper's pairs
  - the SIFT and PolyPhen prediction matrices
- Flattened: the UTR objects, to their sequence strings, and each protein feature's analysis, to the label Ensembl VEP's `DOMAINS` column prints.
- Dropped:
  - the Perl object references (the slice, the gene object, the exon coordinate mapper's trees)
  - the per-exon and per-intron strand, which is the transcript's
  - a sequence edit's description
  - the dump-time flags

A cache written with the Sereal serialiser (`serialiser_type sereal` in its `info.txt`, or `.sereal` shards) is refused with an error naming the serialiser. Only Storable caches convert. A run in which any shard fails to read or write stops with an error before `info.json` is written.

Ensembl VEP's `convert_cache.pl` rewrites its own cache's region variation files as one tabix-indexed `all_vars.gz` per chromosome and, with `--sereal`, its transcript shards as Sereal. It has no counterpart and none is needed: the converter copies a cache's variation files in either form. `vep-cache-converter` is a different tool, which copies a flat JSON transcript cache into the per-chromosome runtime layout.

## Variation files

The variation files are copied as they are into `variations/<chr>/`. They are each chromosome's `all_vars.gz` with its tabix index or, for a cache dumped without the tabix conversion, its `<start>-<end>_var.gz` region files. `all_vars.gz` is the bgzip file the Ensembl VEP dumper writes. It has the chromosome first, `.` for an empty cell and the columns of the cache's `variation_cols` line.

vep-rs opens each chromosome's `all_vars.gz` on first use and reads it through its index. It parses a region file whole the first time a variant falls in it. `--check_existing` and the frequency flags therefore read the converted cache as Ensembl VEP reads the original. A cache whose `variations/` directory holds no such file warns once that the co-location columns stay empty.

A conversion made by a converter older than 0.5.0 carries no `variations/<chr>/all_vars.gz`. Convert the cache again, or copy each chromosome's `all_vars.gz` and its `.csi` from the Ensembl VEP cache into `variations/<chr>/`.

A region file writes an empty cell for a value the record lacks. Ensembl VEP reads an empty `clin_sig_allele` or `clinical_impact` cell inside a line as a value. `CLIN_SIG` then lists no term for the record under the default `--clin_sig_allele 1`, and `--clinvar_somatic_classification` writes its key with an empty value. vep-rs reads the two cells the same way.

## Regulatory features

The cache's regulatory build, the `<start>-<end>_reg.gz` files beside the transcript files, converts into `regulatory/<chr>/<start>-<end>.json`. That is one shard per region, holding the region's regulatory features and transcription-factor binding motifs with their binding matrices. `vep-cache-converter` copies the shards, as `.json` or `.json.gz`, and both tools set the `regulatory` marker of `info.json` only once the tree has been written. A conversion made without the tree carries none: `--regulatory` then warns once and writes no regulatory row ([Regulatory and motif features](regulatory.md#the-cache)).

## Other species

Any species with an Ensembl or Ensembl Genomes cache converts the same way. [Other species](cache-species.md) gives the conversion and a run for each genome the golden corpora cover: mouse, the fruit fly, Arabidopsis and yeast. It also gives the genetic codes and chromosome names of each cache.

## RefSeq and merged caches

The RefSeq and merged caches Ensembl ships (`homo_sapiens_refseq` and `homo_sapiens_merged` under the cache root) convert the same way. The converter writes `source_type` into `info.json` from the species directory's suffix: `refseq`, `merged` or `ensembl`, a key of its own with no counterpart in Ensembl VEP's `info.txt`. A run whose `--refseq` or `--merged` disagrees with the declared set is refused before any input is read, where Ensembl VEP filters whatever models the cache holds by the flags alone; a cache converted without the key is checked against nothing.

Under `--json_cache` the flags name the set of the cache given, so `--refseq` on a converted Ensembl cache without the `source_type` key keeps only the transcripts that pass its test of a RefSeq id, as Ensembl VEP does. Of a cache of `ENST` ids those are the transcripts on `MT`, where the test passes an id made of capitals and digits, so a record on any other chromosome gets an intergenic row.

A merged cache's transcripts carry `source_cache`, `Ensembl` or `RefSeq`, exactly where the Storable has it: on the region files both caches held. A region file one cache alone held carries no key, and its transcripts print no `SOURCE` in either engine.

The `bam` line of a RefSeq or merged cache's `info.txt`, the alignment file its models were corrected against, is carried too. A cache with it is a BAM-edited cache. A run on one prints `INFO: BAM-edited cache detected, enabling --use_transcript_ref; use --use_given_ref to override this`, reads each row's reference allele from the transcript and adds the `GIVEN_REF`, `USED_REF` and `BAM_EDIT` columns, as Ensembl VEP does. `--use_given_ref` keeps the input's allele and drops `GIVEN_REF` and `USED_REF`; `BAM_EDIT` stays on an edited transcript's rows in the default and JSON formats, with no header line, as Ensembl VEP writes it.

A model the alignment corrected carries `bam_edit_status` `ok` beside its `_rna_edit` attributes and its edited cDNA (`spliced_seq`). Its exon coordinate mapper's pairs fold the insertions and deletions in, read from the Storable mapper object rather than rebuilt from the exons. `BAM_EDIT` prints `OK` on its rows and its reference allele is cut from the edited cDNA, coding or not. A model the correction failed on, status `failed`, keeps its own attributes, which can include an `_rna_edit`, loses only the alignment's edits and has no cDNA cached: `BAM_EDIT` prints `FAILED` and its rows read the genome in vep-rs. No editing happens at run time: `--bam`, which edits the models against an alignment file as the cache is read, is refused.

### A transcript model the cache carries at two loci under `--refseq` or `--merged`

One behaviour of Ensembl VEP 116.2 on a RefSeq or merged cache is not reproduced. The format comparator sets it aside under the class named below, as [Intended divergences from Ensembl VEP](intended-divergences.md) lists it; its [How to read the counts](intended-divergences.md#how-to-read-the-counts) defines the terms. The class carries its tier under that page's rule: tier 1 is a shape the chromosome 21 population cuts carry with no flag, tier 2 one that needs a flag or plugin on ordinary data, tier 3 one ordinary data does not carry. The counts are over the cuts' 15,646,481 records, the 111 golden corpora and their 392 flag sets.

**The defect in one sentence.** Under `--refseq` or `--merged`, Ensembl VEP keeps one transcript per stable identifier among the models a buffer loads, so a model the cache carries at two loci keeps only the copy listed first whenever the buffer's regions carry both.

**Why this is a defect.** The method's own description says it "uniquifies" the list by `dbID`, the key that tells two features apart, and its first pass does that. A second pass for a RefSeq or merged source then keeps one transcript per stable identifier, the `ensembl`-sourced copy, else the lowest `dbID`: a rule for one model loaded from two sources at one locus. Run over every transcript of the regions the buffer loads, it also drops a distinct feature, with its own `dbID`, strand and exons, that merely shares the identifier, and which copy survives depends on what else the buffer holds.

**What vep-rs does instead.** vep-rs keeps every transcript the cache distinguishes by source and `dbID`, so both loci of such a model are annotated whatever the buffer holds.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2):

- `Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm`, lines 254-264 and 277-299, `merge_features`: the description and the `dbID` pass; lines 313-343: the pass per stable identifier, line 330 the `ensembl` preference, line 337 the lowest `dbID`.
- `Bio/EnsEMBL/VEP/AnnotationSource.pm`, line 139, `get_all_features_by_InputBuffer`: the method runs over each buffer's features.

**What correcting it would change.** The second pass would collapse only copies at one locus, and a model mapped to two places would keep both.

**The class it produces.** `merged_same_stable_id_collapse`, which the format comparator applies when a comparison case names it. It accepts a transcript row vep-rs alone writes under `--refseq` or `--merged`, its `SOURCE` `RefSeq` or absent, on a record Ensembl VEP wrote a transcript row for, of a model Ensembl VEP wrote no row for on that record. A record whose every model is a lost copy, another allele's row of a kept model and a duplicated row stay differences, and the JSON record's `most_severe_consequence` is left uncompared.

The chromosome 21 of the GRCh38 merged cache carries 325 stable identifiers at more than one locus, every one a `compmerge.*` model that only `--all_refseq` admits, and the chromosome 21 GRCh38 RefSeq cache carries the same 325; the GRCh37 chromosome 21 caches carry none. The `GRCh38-merged` corpus's 76 records under `--merged --all_refseq` lose 128 rows on 59 identifiers in Ensembl VEP and none in vep-rs, with no other row differing; no committed flag set runs that pair. The `GRCh38-refseq` corpus's `--all_refseq` run keeps every row only by record placement: no record of that corpus loads the regions of its two dual-locus models' other copies.

**How often the shape occurs.** Tier 2: the shape is a transcript whose stable identifier its cache carries at two or more loci. The count is over the caches committed with the golden corpora, not over the full caches of the paragraph above: of their 77,675 transcripts, 38 have the shape (489 per million), all in one corpus's cache. The buffer that loads both copies keeps the one listed first. One flag set runs `--all_refseq`, the GRCh38 RefSeq corpus's `--refseq --all_refseq`.

#### `compmerge.264.Heart.chr21` at two loci of chromosome 21

**In short.** The `GRCh38-merged` golden corpus's cache holds `compmerge.264.Heart.chr21` twice, both copies from the RefSeq side: `dbID` 1733463 on the reverse strand at 32601667-32612845 and `dbID` 1631814 on the forward strand at 39451905-39515506. Under `--merged --all_refseq` over the corpus's 76 records, one buffer, vep-rs writes `non_coding_transcript_exon_variant` rows on that identifier for `rs1036535345` at `21:32612111` and for `rs745798146` at `21:39462549`, with `STRAND=-1` and `STRAND=1`. Ensembl VEP 116.2's second pass keeps the copy its buffer lists first, both carrying the source `ensembl`, and writes no row of the identifier for the records at the other locus; under `--merged` alone neither engine admits the model.

A model whose first-listed copy overlaps no record loses every row: `compmerge.109.Heart.chr21` sits at 21:16194573-16607134 and at 21:17513092-17593914, no record lies in the first span, and Ensembl VEP writes no row for `21_17535482_G/T` on the second.

## The published conversions

The published concordance figures were measured on such a conversion ([scripts/README.md](../scripts/README.md#reproducing-the-published-concordance)), which holds the transcript set identical between the two engines. The caches those figures were measured on carry 195,232 transcripts (GRCh37) and 505,231 (GRCh38), the canonical chromosomes only. The GRCh37 conversion is of the cache-version-113 Storable that Ensembl VEP release 115 ships for that assembly. It differs from the cache-version-115 build in 310 of those transcripts (the `storable_113_*` columns of `manuscript/data/transcript_set_parity.csv`).

The conversion carries each translation's sequence edits, and vep-rs applies them the same way. The edits (`seq_edits`) are the initial-methionine, selenocysteine and amino-acid-substitution edits Ensembl VEP applies to the reference peptide at annotation time.

A cache converted without that key, the published conversions among them, still loads. On such a cache, a variant at codon 1 of a transcript with an alternative initiation codon is called from the cached peptide rather than from the edit. Converting the cache again restores Ensembl VEP's call.

A natively built cache carries no edit list either. Its peptide is the protein FASTA's sequence, with the edits already applied, and the builder writes none.

## Where next

- [Cache setup](cache-setup.md): building a cache natively, naming it, the reference FASTA and the chromosome names.
- [Other species](cache-species.md): the conversion command and a run for mouse, the fruit fly, Arabidopsis and yeast.
- [Regulatory and motif features](regulatory.md): the rows and fields a run writes from the cache's regulatory tree.
- [Intended divergences](intended-divergences.md): where the two engines differ, and why.
