# The stop bases of a RefSeq mitochondrial model whose stop codon an edit completes

**The defect in one sentence.** On a transcript model whose `_rna_edit` appends bases after an exon's last genomic base, Ensembl VEP 116.3 pairs the genomic bases the edit replaces with the last bases of the replacement, not the first. A record at those bases is read one or two positions too far along the coding sequence. `MT 9990` on the RefSeq model `COX3` is the stop codon's first base, c.784; it is written at `CDS_position` 786, the codon's third base, as `stop_retained_variant`, where vep-rs writes `CDS_position` 784 and `stop_lost`.

**Why this is a defect.** Six mitochondrial RefSeq models (`ND1`, `ND2`, `COX3`, `ND3`, `ND4`, `CYTB`) end one or two bases into their stop codon: the genome holds the codon's first base or two, and the rest is added after the transcript's last genomic base. The cache records that as an `_rna_edit` replacing the exon's last one or two bases by the whole codon. On each model the reference sequence at the replaced positions equals the codon's first bases and no other placement of them, so a genomic stop base is the codon's first or second base. Reading it as the second or third applies the record to the wrong base, and the term follows the wrong codon.

**What vep-rs does instead.** With a reference sequence, vep-rs reads the genome's bases at the replaced positions when it loads the transcript. Where the completed stop codon begins with them and no other placement of them inside the codon is an identity, it pairs them with the codon's first bases. The codon's remaining bases stay in the coding sequence and have no genomic position. A run without a reference sequence, and an edit of any other shape, keeps the pairs the cache holds.

**Where it is in Ensembl VEP's source** (`ensembl` at release/116):

- `Bio/EnsEMBL/TranscriptMapper.pm`, `_load_mapper`, line 179: for an edit that changes the cDNA's length, the pair before the edit ends just before the edit's first base.
- The same routine, lines 209-211: for an edit that lengthens the cDNA, the next pair's cDNA start moves past the added length while its genomic start stays at the edit's first genomic base, so the replaced genomic bases meet the replacement's last bases.
- The same routine, line 169: the edit's bases are not compared with the genome.

**What correcting it would change.** `_load_mapper` would place the added bases after the replaced ones when the replacement begins with the genome's bases. A record at a genomic stop base of the six models, or between the last sense base and the stop codon, would carry the positions, codons, residues and terms it carries on a transcript whose stop codon is whole on the genome.

**The class it produces.** `refseq_edited_model_appended_edit_pairs`: differing `cDNA_position`, `CDS_position`, `Codons`, `Amino_acids` and consequence terms on a record's row on one of the six models, and a differing `HGVSc` under `--hgvs`. The concordance datasets are annotated on the Ensembl caches, which hold no such edit, so the entry is documented on this page and the comparator carries no rule for it.

## 1. A substitution of the first stop base of `COX3`

**In short.** `MT 9990 T A` on the RefSeq GRCh38 cache under `--refseq` changes the one stop base the genome holds for `COX3`. Ensembl VEP writes `CDS_position` 786 and `stop_retained_variant`; vep-rs writes `CDS_position` 784 and `stop_lost`.

### 1.1 What a reader needs to know first

Some mitochondrial transcripts end before their stop codon does. The genome holds the codon's first base, or its first two, and the bases added to the transcript's 3' end complete it. A RefSeq cache holds such a model with its whole coding sequence, the completed stop codon included, and with one `_rna_edit` attribute that says which of the exon's last bases the completed codon replaces.

A mapper pairs each base of the model with a base of the genome. The completed codon has three bases and the genome has one or two of them, so one or two of the codon's bases have no genomic base. Which ones is the question this page answers from the reference sequence.

### 1.2 The record

```text
MT	9990	stop_first_base	T	A	.	.	.
```

`COX3` has one exon, `MT:9207-9990`, 784 bases. Its `_rna_edit` replaces cDNA 784 by three bases, and its coding sequence has 786 bases: 262 codons, the last the stop codon at c.784-786.

### 1.3 What each engine writes

| Field | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- |
| `Consequence` | `stop_retained_variant` | `stop_lost` |
| `cDNA_position` | 786 | 784 |
| `CDS_position` | 786 | 784 |
| `Protein_position` | 262 | 262 |
| `Amino_acids` | `Y/*` | `*/K` |
| `Codons` | `taT/taA` | `Taa/Aaa` |

### 1.4 Why Ensembl VEP writes that

The mapper of `COX3` holds two pairs: cDNA 1-783 with `MT:9207-9989`, and cDNA 786 with `MT:9990`. `_load_mapper` ends the first pair before the edit and starts the second two bases later on the cDNA, at the same genomic base. `MT:9990` is therefore read as c.786, the stop codon's third base.

### 1.5 Why vep-rs's row is the correct one

The reference base at `MT:9990` is the completed codon's first base and differs from its second and third. The genome's base is therefore c.784, and c.785-786 are the bases added after the transcript's last genomic base. The record changes the first base of the stop codon, which leaves no stop codon under the vertebrate mitochondrial code: `stop_lost` at c.784.

The same substitution on the first stop base of `ATP6`, whose stop codon is whole on the genome, is `stop_lost` at its own position in both engines.

### 1.6 Try it yourself

```bash
printf 'MT\t9990\tstop_first_base\tT\tA\t.\t.\t.\n' > record.vcf
vep -i record.vcf -o vep_rs.txt --offline --refseq \
  --json_cache refseq_GRCh38_json_cache --assembly GRCh38 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --force_overwrite --no_stats
grep COX3 vep_rs.txt | cut -f 7-12
```

The same record through Ensembl VEP 116.3 with `--cache --refseq` on the RefSeq GRCh38 cache of release 116 gives the other column.

### 1.7 How the comparator treats this class

It does not: the class has no comparator rule. The concordance datasets are annotated on the Ensembl caches, whose models of the same six genes carry no edit and pair each stop base with its own position in both engines.

## The other rows of the six models

On 672 records that touch the seven genomic stop bases of the six models, or the bases beside them, vep-rs's terms are the terms the same changes get on six mitochondrial models whose stop codon is whole on the genome.

| Kind of record | Rows | Rows differing from Ensembl VEP |
| --- | --- | --- |
| Substitution of a stop base | 21 | 21 |
| Deletion of the stop bases | 8 | 8 |
| Deletion from the last sense base into the stop codon | 7 | 7 |
| Insertion between the last sense base and the stop codon | 504 | 504 |
| Insertion between the two stop bases of `ND1` | 84 | 84 |
| Substitution or deletion of the last sense base | 24 | 0 |
| Insertion after the transcript's last base | 24 | 0 |

The consequence terms differ on 575 of the 624 differing rows; on the other 49 the terms agree and the positions differ. Under `--hgvs` the `HGVSp` follows the same pairs: Ensembl VEP prints one where vep-rs prints none on 271 of the 672 rows, and both print one in differing cells on 211; vep-rs reads each record at the place the reference sequence gives. vep-rs prints an `HGVSp` where Ensembl VEP prints none on 3 of the 672 rows, each a deletion of one or two stop bases the genome holds: 2 on `ND1` by this entry's pairs and 1 on `CYTB` by the class [`hgvs_shift_sequence_edge_window`](hgvs-shift-near-sequence-edge.md).
