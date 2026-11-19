# `--hgvsg_use_accession` in an offline run

**The defect in one sentence.** Under `--hgvsg --hgvsg_use_accession`, an offline run of Ensembl VEP 116.2 looks for a record's RefSeq accession among the synonyms of the chromosome name exactly as the input writes it. On a name the run resolves that has no `NC_` synonym of its own, the run dies inside the Ensembl API with `Can't call method "db" on an undefined value`, and writes no row for the record or for any record after it. vep-rs writes the chromosome's accession for every spelling it resolves, and every record's rows.

**Why this is a defect.** The flag asks for the accession of the record's chromosome, which does not depend on how the input spells the chromosome. Ensembl VEP reads `Chr21`, `CHR21` and `NC_000021.9` as chromosome 21 and annotates each against it. Its synonyms table pairs `21` with `NC_000021.9`, and with a database connection it prints that accession for all three names. Offline, the lookup takes one step in the table from the name as written. It finds the accession for `21`, `chr21` and `CM000683.2` and for no other name of the chromosome, and the accession itself is among the names it misses.

The miss is not handled. The lookup passes an undefined name on, the HGVS code then asks the chromosome's slice for its synonyms, and the slice of an offline run has no database to ask. The run ends on a Perl error that names neither the flag nor the record.

**What vep-rs does instead.** vep-rs reads the accession under the name as written and, where that name has none, under the cache's name for the chromosome the record resolves to, the route it takes for `--ga4gh_vrs`. It writes `NC_000021.9` for every spelling of chromosome 21 it resolves, the name Ensembl VEP writes with a database connection, and it writes every record's rows. [Output control](../cli/output-control.md) describes the flag.

A scaffold or a patch that the FASTA holds takes the RefSeq name the table holds for it, and a sequence the table gives no RefSeq name takes the cache's name for it. Ensembl VEP's connected run writes those names on a human scaffold, a human patch and a yeast chromosome.

The table carries no label for the source of a name, and the connected run takes a RefSeq name only from the source its database labels `RefSeq_genomic`. Where the database files a sequence's RefSeq names under another label, as the fly database does, the connected run writes the sequence's own name and vep-rs a RefSeq accession from the table, the versioned form where the table lists both. Under `--ga4gh_vrs` vep-rs writes its object there on a sequence's `NC_` accession, where Ensembl VEP writes none, offline or with a database.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/` unless a path says otherwise):

- `Config.pm`, line 199: the flag, under the comment `force HGVSg to return on chromosome accession instead of input chr name`.
- `BaseVEP.pm`, lines 637-665: the cache's `chr_synonyms.txt` is read into a table keyed by name. Lines 653-654 store the pair of each line under both of its names, and nothing joins two pairs.
- `OutputFactory.pm`, lines 1368-1372: without a database, the accession is the first name matching `NC_` among the synonyms of `$vf->{chr}`, the name as the input writes it. Line 1372 passes it to `hgvs_genomic`, undefined when there is none.
- `Bio/EnsEMBL/Variation/VariationFeature.pm` of `ensembl-variation` at release/116, lines 1912-1916: given no name, `hgvs_genomic` asks the slice for its `RefSeq_genomic` synonyms.
- `Bio/EnsEMBL/Slice.pm` of `ensembl` at release/116, line 2985: `get_all_synonyms` calls `$self->adaptor->db`, and the slice of an offline run has no adaptor.

**What correcting it would change.** A lookup that followed the written name to the chromosome the run annotates the record against would find `NC_000021.9` for every name of chromosome 21. An offline run would then write the rows its connected run writes, which are vep-rs's rows.

**The class it produces.** `hgvsg_use_accession_ends_the_run_offline`: every row of an offline run under `--hgvsg --hgvsg_use_accession` from the first record the run resolves whose chromosome name has no `NC_` synonym of its own. The reference, Ensembl VEP's side of a comparison, ends at that record, with exit status 255 under `--fork` and 25 without it, where the candidate, vep-rs's side, writes the record's rows and every later record's. No comparator codes the class. The one committed flag set that passes both flags, `accession` of the `GRCh38-identifiers` golden corpus, names each of its 47 records `21`, and the engines agree on every row of it.

**How often the shape occurs.** Tier 3, see [Intended divergences from Ensembl VEP](../intended-divergences.md): the shape is a chromosome name outside the three names of each chromosome that Ensembl's human synonyms table pairs with its accession, in a run under both flags. No record of the cuts carries one, since each cut names the chromosome `21` or `chr21`. Of the 392 committed flag sets, 1 passes both flags, and 10 golden corpora hold a record under such a name.

## 1. One substitution on a chromosome written `Chr21`

**In short.** The record `Chr21 10538885 A C` is a substitution in an intron on chromosome 21, under a spelling both engines resolve and Ensembl's synonyms table does not hold. Under `--hgvsg --hgvsg_use_accession`, an offline run of Ensembl VEP 116.2 dies with `Can't call method "db" on an undefined value`, its header written and no row; its exit status is 255 under `--fork` and 25 without it. vep-rs writes the record's nine rows with `HGVSg=NC_000021.9:g.10538885A>C`, the rows Ensembl VEP writes with a database connection.

### 1.1 What a reader needs to know first

`--hgvsg` adds the `HGVSg` field, the variant in HGVS notation on its chromosome. The notation opens with the chromosome's name as the input writes it, as in `21:g.10538885A>C`. `--hgvsg_use_accession` puts the chromosome's RefSeq accession in that place, `NC_000021.9` for chromosome 21 of GRCh38.

Each cache carries a table of chromosome synonyms, `chr_synonyms.txt`, one pair of names a line. Ensembl's table for GRCh38 gives chromosome 21 four names: `21`, `chr21`, `CM000683.2` and `NC_000021.9`. A run reads a record under any of them, and under a name the table lacks whose `chr` prefix, in any letter case, comes off to leave the cache's name ([Chromosome names](../cache-setup.md#chromosome-names)).

An Ensembl VEP run reads its annotation from a local cache, from a database, or from both, and `--offline` forbids every database connection.

### 1.2 The record

The record is `intron_snv_1` of the `GRCh38-identifiers` golden corpus with its chromosome written `Chr21`. It is read against that corpus's cache and reference sequence (`tests/golden/116/GRCh38-identifiers/`):

```
Chr21	10538885	intron_snv_1	A	C	.	.	.
```

Both engines run offline under the same two flags:

```
--hgvsg --hgvsg_use_accession
```

### 1.3 What each engine writes

Ensembl VEP 116.2 writes its header and no row and writes this line to standard error:

```
Can't call method "db" on an undefined value at /opt/vep/src/ensembl-vep/Bio/EnsEMBL/Slice.pm line 2985, <__ANONIO__> line 3.
```

Run with `--fork 1`, as the reproduction runs it, the line `Died in forked process` and a process number follow, and the exit status is 255. Run without `--fork` it writes the line alone and exits with status 25.

vep-rs writes the record's nine rows, one for each of nine transcripts, each with the accession. The first two:

```
intron_snv_1	Chr21:10538885	C	ENSG00000274391	ENST00000427445	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C
intron_snv_1	Chr21:10538885	C	ENSG00000274391	ENST00000612746	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1;HGVSg=NC_000021.9:g.10538885A>C
```

Run with `--cache` and without `--offline`, so that it may connect to Ensembl's public database, Ensembl VEP writes vep-rs's nine rows byte for byte.

The name decides the outcome. Offline, both engines write `HGVSg=NC_000021.9:g.10538885A>C` on the nine rows of the record named `21`, `chr21` or `CM000683.2`. Ensembl VEP dies the same way on `CHR21` and on `NC_000021.9`, where vep-rs writes the accession. In a longer input Ensembl VEP writes the rows of the records above the first such name and ends at it. Under `--hgvsg` alone both engines write `HGVSg=Chr21:g.10538885A>C`.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP maps the record's name to the cache's chromosome when it reads the cache, and the record keeps the name as written. Offline, the accession lookup reads the synonyms table under that written name and takes the first synonym matching `NC_`. The table holds the pair of each line under both of its names and never joins two pairs. Only a name that shares a line with the accession finds it: `21`, `chr21` and `CM000683.2`. `Chr21` stands on no line, and the three lines of `NC_000021.9` pair it with no second accession.

The lookup calls the HGVS code with what it found, a name or nothing. Given nothing, that code asks the record's slice for its `RefSeq_genomic` synonyms, the route a connected run takes for every record under the flag. The slice of an offline run has no database adaptor, and the call on it is the error of 1.3.

Source (`ensembl-vep` at release/116.2, `Bio/EnsEMBL/VEP/` unless a path says otherwise):

- `BaseVEP.pm`, lines 689-738: a name is mapped to the cache's chromosome through a synonym or, at lines 722-726, with its `chr` prefix removed in any letter case.
- `BaseVEP.pm`, lines 653-654: the table stores each pair under both names.
- `OutputFactory.pm`, lines 1368-1372: the offline lookup and the call that passes its result on.
- `Bio/EnsEMBL/Variation/VariationFeature.pm` of `ensembl-variation` at release/116, lines 1912-1916, and `Bio/EnsEMBL/Slice.pm` of `ensembl` at release/116, line 2985: the request for the slice's synonyms and the call that fails.

### 1.5 Why vep-rs's row is the correct one

Position 10538885 of the record lies on chromosome 21 of GRCh38 whatever the input calls the chromosome. Both engines annotate it against the transcripts of chromosome 21, and Ensembl VEP writes the nine rows itself once `--hgvsg_use_accession` is left off. `NC_000021.9` is that chromosome's RefSeq accession in Ensembl's own synonyms table, and Ensembl VEP prints it for this record whenever it can reach its database. The flag asks for that name. An exit that loses the record and every record after it answers nothing.

### 1.6 Try it yourself

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --fasta <corpus>/reference.fa.gz --hgvsg --hgvsg_use_accession
```

vep-rs writes its rows of 1.3. Ensembl VEP 116.2 takes the same two flags with a FASTA file and the cache options of [Reproducing these examples](reproducing.md), and ends with the line of 1.3. With `--cache` and a connection to Ensembl's public database in place of `--offline`, it writes vep-rs's rows.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/hgvsg_use_accession_ends_the_run_offline.toml` runs the record of 1.2 on both engines offline and records the location, transcript and `HGVSg` of each row: no row and exit status 255 for Ensembl VEP, the nine rows with the accession for vep-rs.

No comparator sets the class aside, since no committed run meets it: the `accession` flag set of the `GRCh38-identifiers` golden corpus names its records `21`, and vep-rs's rows there equal Ensembl VEP's.
