# RefSeq and merged cache columns

A converted RefSeq or merged cache carries the attributes Ensembl VEP writes about each RefSeq model's alignment to the genome.

`--refseq` and `--merged` add the columns that print them. `--refseq` adds `REFSEQ_MATCH` and `REFSEQ_OFFSET` after `RefSeq`; `--merged` adds `REFSEQ_MATCH`, `SOURCE` and `REFSEQ_OFFSET`. Both sit before `GIVEN_REF`, `USED_REF` and `BAM_EDIT`, each name listed once when both flags are on. The JSON keys are `refseq_match`, `source`, `refseq_offset` and `bam_edit`; the Parquet columns carry the column names.

## What `REFSEQ_MATCH`, `SOURCE` and `BAM_EDIT` print

`REFSEQ_MATCH` lists the model's `rseq_*` match codes in the order the cache holds them. The codes are `rseq_mrna_match`, `rseq_5p_mismatch`, `rseq_cds_mismatch` and the others of the Ensembl core attribute set. They are joined with `,` in the default and tab formats and `&` in the VCF, and are an array in JSON. A transcript with no such code prints nothing. The cache also carries the `enst_refseq_compare` attribute of a model compared with an Ensembl transcript, which neither engine prints.

`SOURCE` names the cache a merged cache's transcript came from, `Ensembl` or `RefSeq`; a single-source cache stamps no transcript, so the column stays empty there.

`BAM_EDIT` prints `OK` or `FAILED` for a model the cache's alignment file touched ([the transcript's reference allele](../output-formats.md#the-transcripts-reference-allele)).

## When `REFSEQ_OFFSET` is set

`REFSEQ_OFFSET` is set under `--hgvs` on an `NM_` or `XM_` model whose alignment inserted or deleted bases upstream of the allele. Its value is the sum of those edits' length changes, positive for inserted bases, over the edits starting at or before the allele's cDNA start. It prints only when it is not 0, so a model whose edits cancel out, an allele upstream of every edit and an intronic allele, which has no cDNA start, carry none.

## HGVSc on an edited model

The HGVSc of a model whose alignment inserted or deleted bases is numbered against the edited model, the transcript the identifier names. An exonic position takes its edited coordinate and an intronic distance counts from the edited exon end. A genomic base the alignment removed from the model is described from the nearest modelled base, as an intronic base is. A single-base substitution in a coding exon takes the CDS coordinate the edited mapper gives. A model whose correction failed (`BAM_EDIT=FAILED`) keeps none of the alignment's edits, and its HGVSc is counted along the genome's exons ([the position past such a model's own edit](../divergences/refseq-failed-edit-model-hgvsc-offset.md)).

The reference and alternate alleles are spelled relative to the transcript, the reference being the allele the row was analysed against. So under the auto-enabled `--use_transcript_ref` a substitution names the transcript's base where the genome differs.

An alternate allele equal to that base gets no HGVSc, as any allele equal to its reference gets none. Ensembl VEP writes the base changed to itself there (`n.2023T>T`), the documented class [`edited_model_no_change_hgvsc`](../divergences/edited-model-no-change-hgvsc.md).

A shifted insertion's inserted bases are spelled as the shift rotates them, with or without `--shift_3prime 1`. Ensembl VEP spells them as written except on a coding row under that flag. The difference is an intended divergence ([the inserted bases of a shifted insertion's HGVSc on an edited model](../divergences/refseq-edited-model-shifted-insertion-bases.md)).

An insertion or a deletion on such a model is shifted 3' along the genome before its notation is written, and along the model's own sequence only on bases the genome lacks. A run that writes such a notation counts those alleles in one `INFO:` line at its end, withheld under `--quiet`.

Ensembl VEP's source holds a second shift along the edited sequence (`TranscriptVariationAllele.pm` lines 211-256). Away from a transcript's last bases and its edited bases it returns before that shift (lines 182-187), copying a shift it already has ([the HGVS_OFFSET of an insertion or deletion on an edited model](../divergences/refseq-edited-model-hgvs-offset.md)). Its GenBank protein name (lines 1750-1756) needs the transcript's analysis record, which a cache does not carry, so both engines name the cached protein identifier. An LRG name takes no version in either engine.

Ensembl VEP numbers these positions differently. Its cache dumper strips the transcript's edit-aware mapper, so it walks the exons in their unedited numbering and corrects an exonic position on an `NM_` or `XM_` model by the offset above. An intronic allele, which has no offset, is anchored on the unedited exon end, and an Ensembl model with such an edit is never corrected. The difference is an intended divergence ([intended divergences](../intended-divergences.md)); the comparator masks the class.

## A variant inside bases the alignment deleted

Genomic bases the alignment deleted from the model carry no transcript coordinate. A variant wholly inside them takes the term of the transcript region whose span holds it, with no cDNA, CDS or protein position:

- `coding_sequence_variant` in the coding sequence;
- the UTR term in a UTR;
- the exon term on a model with no coding sequence.

The run writes one warning per such model, `Bases deleted by the alignment of <transcript> carry no transcript coordinate; a variant inside them takes the term of the transcript region whose span holds it`. Ensembl VEP writes `intergenic_variant` on a coding model's rows, the default on a mapper gap. A non-coding model's exon term reads the exon's genomic span on both engines.

The difference is an intended divergence ([A variant inside the bases a RefSeq model's alignment deleted](../divergences/refseq-edited-model-deleted-span-term.md)); the comparator masks the class `refseq_edited_model_deleted_span_term`.

## The transcript-assembly mismatch warning

An alternate allele equal to the bases of the model's first edit, in a coding exon, is the variant Ensembl VEP reads as invalid. The run then writes `WARNING: Transcript-assembly mismatch in <Uploaded_variation>` to standard error and the warnings file, once per record, for a row the output prints. A marked row that `--pick` or another record-level flag drops draws nothing, as do `--most_severe` and `--summary`, which print no transcript row, and a run under `--use_given_ref`. Ensembl VEP under `--fork` writes the line behind the child process's id, a documented intended divergence, [`fork_warning_pid_prefix`](../divergences/fork-warning-pid-prefix.md).

In a `--json` run Ensembl VEP writes the line with no record name, once for the run or for each forked child. vep-rs names each record in every format, a second documented intended divergence, [`json_mismatch_warning_names_no_record`](../divergences/json-mismatch-warning-names-no-record.md).

## Where next

- [The transcript's reference allele](../output-formats.md#the-transcripts-reference-allele): how `GIVEN_REF`, `USED_REF` and `BAM_EDIT` are set.
- [An HGVSc naming a change of a base to itself](../divergences/edited-model-no-change-hgvsc.md): the record and what each engine writes.
- [A variant inside the deleted bases](../divergences/refseq-edited-model-deleted-span-term.md): the term each engine writes on such a variant.
- [Intended divergences](../intended-divergences.md): where the two engines differ, and why.
