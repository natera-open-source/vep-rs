# A name the synonyms file gives to several sequences

**The defect in one sentence.** Where a cache's synonyms file gives one name to two or more sequences the cache holds, Ensembl VEP 116.3 reads a record under that name on whichever of them comes first in Perl's hash order. That order changes between runs of the same command. A breakend record in the same buffer whose mate lies on another of those sequences can move the record to the mate's sequence, as section 1.4 shows. vep-rs skips the record with a warning that names the sequences, since the file does not say which of them the name stands for.

**Why this is a defect.** A record's rows are those of the sequence it is read on. Ensembl VEP reads the record on the first of the name's synonyms the cache holds, in the order Perl's `keys` returns them. Perl's documentation of `keys` says that order differs between runs: the hash seed is random for each process unless `PERL_HASH_SEED` is set. One command over one cache and one input therefore writes the rows of one sequence in one run and of another in the next, with no warning in either. Rows that change with a hash seed are not a function of the input.

**What vep-rs does instead.** vep-rs resolves a name the same way in every run. A name the cache holds is itself, and a name with one synonym the cache holds is that synonym. A name with two or more is left as written, so no cache lookup matches it. The record is skipped with the warning `Chromosome <name> is a synonym of <n> chromosomes of the cache (<first>, <second>, ...) and is read as none of them`. The warning names the first two in the file's order, and `, ...` follows them only where there are more than two.

A record under such a name passes through the same checks as a record on a chromosome the cache lacks, in the same order. Under `--dont_skip`, and for a structural record whose span is over `--max_sv_size`, it is kept as written and gets its `intergenic_variant` row.

The entry covers the records of a VCF, an Ensembl-format and a region input, whose names the cache resolves. An SPDI notation is resolved by its own reader, and a `--custom` file and a GFF3 or GTF read beside a cache key their own sequences. Each of those still takes the first of several sequences in the file's order. An HGVS notation and a regulatory GFF resolve such a name by rules of their own, and a record that one of these files lists is admitted, not skipped.

The reference FASTA has a reader of its own too, which reaches a record through its aliases and takes the first of several by name. A record under such a name that is skipped is not read against a sequence. One that is kept as written, or a `<CNV:TR>` record, which is expanded before the chromosome check, can still be read against one of them.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, under `Bio/EnsEMBL/VEP/`; the same lines stand at release/116.2):

- `CacheDir.pm`, line 254: the cache's `chr_synonyms.txt` is the synonyms file where the command gives no `--synonyms`. `Runner.pm`, line 131, reads it.
- `BaseVEP.pm`, lines 646-656: each line of the file is split on whitespace, and the first field and each later field are recorded both ways as keys of a hash of hashes.
- `Parser/VCF.pm`, lines 339-350: a substitution keeps the name as written (line 348). `Parser.pm`, lines 770-802: the name passes the check, the valid names being widened by their synonyms (line 782).
- `AnnotationSource.pm`, line 179: the cache source asks `get_source_chr_name` for the name it reads region files under. `AnnotationSource/Cache.pm`, lines 90-92: the valid names are the cache's.
- `BaseVEP.pm`, lines 711-716: the loop walks `keys` of the name's synonyms and takes the first the cache holds. Lines 695-697 and 734 keep the mapped name for the rest of the run.
- `InputBuffer.pm`, lines 348-357, reached from `AnnotationType/Transcript.pm`, line 108, through `InputBuffer.pm`, lines 286-345: the record matches a loaded transcript of any held sequence the name is a synonym of. `AnnotationType/Transcript.pm`, lines 122-143, keeps a record on the sequence of the first transcript it is matched with, a breakend record through a mate too (`InputBuffer.pm`, lines 335-342), and keeps a breakend record on the sequences its mates name as well.
- `InputBuffer.pm`, lines 195-200 and 224-234: a buffer ends where the written name changes. `AnnotationSource.pm`, lines 250-257, lists region files of a further sequence only for the mate of a breakend record. In a buffer with no such mate on another of the name's sequences, the region files chosen at `AnnotationSource.pm` line 179 therefore decide every row of a record that is no breakend. With one, lines 119-131 take the listed region files held in memory first, line 134 reads the others, and line 143 filters each sequence's transcripts, as section 1.4 shows.
- `BaseVEP.pm`, lines 596-621, and `Pipeline/DumpVEP/CreateDumpJobs.pm`, lines 200-209: the dump that writes the file, described in section 1.1.

**What correcting it would change.** A loop that walked the name's synonyms in one fixed order would write the same rows in every run of a command. It would still read the record on a sequence the file does not single out. The engines would agree on such a record only if Ensembl VEP refused the name too, or if the file gave the name to one held sequence.

**The class it produces.** `synonym_of_several_sequences_by_hash_order`: every row of a record whose name the cache holds no sequence by and the synonyms file gives to two or more sequences the cache holds. Where no breakend record of the record's buffer has its mate on another of them, Ensembl VEP 116.3 writes, for a record that is no breakend, the rows of one of those sequences, chosen by hash order. Section 1.4 gives the other case. vep-rs writes no row and one warning. A comparison that is given the Ensembl VEP cache sets Ensembl VEP's rows aside, as section 1.7 describes.

**How often the shape occurs.** Tier 3: no record of the chromosome 21 cuts is written under such a name. No committed cache's synonyms file, and no committed `--synonyms` file, gives one name to two of a cache's sequences. The difference needs no flag, so each of the 511 committed flag sets meets its requirements.

Real data carries the shape. The table counts the names of `chr_synonyms.txt` in six ready-made Ensembl caches: four of Ensembl Metazoa, cache version 63, and the human caches of release 116.

| Cache | Assembly | Names in the file | Given to two or more held sequences |
| --- | --- | --- | --- |
| *Anopheles coluzzii* | AcolN3 | 518 | 351 |
| *Anopheles funestus* | idAnoFuneDA-416_04 | 1,322 | 581 |
| *Anopheles gambiae* | GCA943734735v2 | 576 | 166 |
| *Drosophila miranda* | D.miranda_PacBio2.1 | 527 | 212 |
| Human | GRCh38 | 3,373 | 0 |
| Human | GRCh37 | 856 | 0 |

One name of the *Anopheles coluzzii* file is given to 105 held sequences. On a record under that name, the identical Ensembl VEP 116.3 command, run four times with no seed set, wrote four different row sets. The human caches hold no such name on either assembly, so no record of the human datasets reaches the difference.

## 1. One name given to chromosomes 21 and X

**In short.** A synonyms file gives the name `shared_name` to chromosomes 21 and X, and one substitution is written under that name. Ensembl VEP 116.3 writes the 2 rows of chromosome X under `PERL_HASH_SEED=0` and the 13 rows of chromosome 21 under `PERL_HASH_SEED=11`, each with `PERL_PERTURB_KEYS=0`. vep-rs writes no row in any run, and one warning that names both chromosomes.

### 1.1 What a reader needs to know first

Each line of a cache's `chr_synonyms.txt` holds names that stand for each other: the first name and each later one, in both directions. `--synonyms <file>` replaces the cache's file ([The synonyms file](../cache-setup.md#the-synonyms-file)). A record under a name the cache holds no sequence by is read on the sequence that name is a synonym of.

Ensembl VEP's dump writes the file from the database's synonyms and links names with one level of indirection (`BaseVEP.pm` lines 596-621). For each database synonym, lines 610-619 link the new name to every name already linked to its partner. The dump writes each linked pair once, in the order `keys` returns them (`Pipeline/DumpVEP/CreateDumpJobs.pm` lines 200-209).

Where one synonym is shared by many sequences, the other names of those sequences are linked to each other. In the *Anopheles coluzzii* cache the name `X` is paired with 114 held sequences, the chromosome and 113 scaffolds.

The file holds unordered pairs alone. A pair that comes from the database and a pair the indirection added cannot be told apart in it. For a name paired with several held sequences, the file therefore does not say which of them the name stands for.

### 1.2 The record

The name `shared_name` is invented. vep-rs reads the cache of the committed corpus `tests/golden/116/GRCh38-synonyms`, which holds chromosomes 1, 21 and X. Ensembl VEP reads the GRCh38 cache of release 116, of which the corpus's cache is a slice. Both engines are given a synonyms file of two tab-separated lines with `--synonyms`:

```
21	shared_name
X	shared_name
```

The input is one record: the corpus's one-base substitution `rs575418327`, at position 28885369 of chromosome 21, written under the invented name.

```
shared_name	28885369	shared_name_record	C	T	.	.	.
```

### 1.3 What each engine writes

vep-rs writes no row for the record in any run and exits 0. Its warnings file holds one line:

```
WARNING: line 1 skipped (shared_name 28885369 shared_name_record C T . . .): Chromosome shared_name is a synonym of 2 chromosomes of the cache (21, X) and is read as none of them
```

Ensembl VEP 116.3 under `PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0` writes 2 rows, those of chromosome X at that coordinate:

```
shared_name_record	shared_name:28885369	T	ENSG00000169306	ENST00000378993	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
shared_name_record	shared_name:28885369	T	ENSG00000169306	ENST00001135211	Transcript	intron_variant	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
```

Ensembl VEP 116.3 under `PERL_HASH_SEED=11 PERL_PERTURB_KEYS=0`, on the same command, writes 13 rows, those of chromosome 21:

```
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00000303775	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00000351429	Transcript	5_prime_UTR_variant	5	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000215005	ENST00000447985	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=1911;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00000460212	Transcript	5_prime_UTR_variant,NMD_transcript_variant	3	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000307311	ENST00000825033	Transcript	non_coding_transcript_exon_variant	13	-	-	-	-	-	IMPACT=MODIFIER;STRAND=1
shared_name_record	shared_name:28885369	T	ENSG00000307311	ENST00000825034	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=56;STRAND=1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00000873666	Transcript	5_prime_UTR_variant	5	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00000922530	Transcript	5_prime_UTR_variant	5	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00001043897	Transcript	5_prime_UTR_variant,NMD_transcript_variant	16	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00001095638	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=18;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00001095639	Transcript	5_prime_UTR_variant,NMD_transcript_variant	7	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00001122419	Transcript	upstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2;STRAND=-1
shared_name_record	shared_name:28885369	T	ENSG00000156239	ENST00001131612	Transcript	5_prime_UTR_variant	5	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
```

From the Allele column on, the 13 rows equal the corpus's committed Ensembl VEP rows for the same position named `chr21`, the record `rs575418327` of `expected/default.txt.gz`. Neither run of Ensembl VEP writes a warning.

### 1.4 Why Ensembl VEP writes that

The two lines of the file give `shared_name` a hash of two keys, `21` and `X`. The record passes the parser's check, since the valid names are widened by their synonyms. The cache source then asks which name it reads region files under. The loop walks the keys of that hash and takes the first the cache holds. The cache holds both, so the first key the hash yields is taken, and it is kept for the rest of the run.

Which key a hash yields first follows the process's hash seed: with `PERL_PERTURB_KEYS=0`, `X` under seed 0 and `21` under seed 11. The region files of that chromosome are loaded, and the record matches a loaded transcript of any held sequence its name is a synonym of. The transcripts of the chosen chromosome at position 28885369 therefore give every row. Without `PERL_HASH_SEED`, Perl draws the seed at random for each process. Source: the lines listed under "Where it is in Ensembl VEP's source".

A breakend record in the record's buffer can change this where its mate lies on another of the name's sequences. With the line below written before the record, Ensembl VEP 116.3 writes the 13 rows of chromosome 21 for the record under both seeds. With the mate at `X:28885400` it writes the 2 rows of chromosome X under both. In both inputs the breakend record itself gets rows on the transcripts of its mate's chromosome alone.

```
shared_name	5000000	mate_on_21	C	C[21:28885400[	.	.	SVTYPE=BND
```

The cache source lists the megabases around the mate too, and the record stays on the sequence of the first loaded transcript it matches. A breakend record is also matched through its mates, and stays on that sequence and on those its mates name. Listed files the previous buffer left in memory are walked first. The others follow in the input's order: a record's own megabases, then its mates', each span widened by the 5,000 bases of the default `--distance`. Here the mate's file comes before the record's megabase of the chosen chromosome. With the breakend at position 28885300 or 28885500 instead, the chosen chromosome's file comes first, and the seed decides the record's rows as in section 1.3.

Of each sequence, Ensembl VEP keeps the transcripts within that distance of the span the buffer covers on it, which for the mate's chromosome can be the mate alone. With the breakend at position 5000000 and its mate at `21:28875000`, the record gets 10 of the 13 rows of chromosome 21 under seed 0 and all 13 under seed 11. With the mate at `21:28100000`, no transcript kept around the mate reaches the record, and the seed decides.

A record written `21` or `X` is read in another buffer, so it ends the buffer of the records written before it, and otherwise acts only through the files it leaves in memory. With `21 28885400` written before the breakend at position 28885300 with its mate at `21:28885400`, the record gets the 13 rows of chromosome 21 under seed 0. vep-rs skips the record in each of these inputs, and each breakend record with it.

### 1.5 Why skipping the record is the correct outcome

The right rows for the record are those of the sequence its name stands for. The synonyms file does not say which that is: it pairs `shared_name` with `21` and with `X`, and marks neither pair. An engine that picks one of several sequences, by hash order or by the file's order, is right for one of them at most and says nothing when it is wrong.

Ensembl VEP's pick also changes from run to run, so its rows are not a function of its input: the 2 rows of chromosome X in one run, the 13 of chromosome 21 in another. One draw in several lands on the right sequence, and the run that does so reads the same as the runs that do not.

vep-rs does not find the right sequence either. It writes no row and says why, the same on every run. The warning names the sequences, so a user can write the record under the cache's own name, or give a `--synonyms` file that pairs the name with one sequence.

### 1.6 Try it yourself

Write the record of section 1.2 to `record.vcf` under a VCF header, as [Reproducing these examples](reproducing.md) describes. vep-rs, from the repository's root, the first line writing the synonyms file:

```
printf '21\tshared_name\nX\tshared_name\n' > synonyms.txt
vep -i record.vcf -o vep-rs.txt --offline \
  --json_cache tests/golden/116/GRCh38-synonyms/json_cache \
  --assembly GRCh38 --synonyms synonyms.txt --no_stats --force
```

vep-rs writes the output file with no data row, and `vep-rs.txt_warnings.txt` holds the warning. Ensembl VEP 116.3, over the GRCh38 cache of release 116:

```
docker run --rm -e PERL_HASH_SEED=0 -e PERL_PERTURB_KEYS=0 \
  -v "$PWD":/work -v /path/to/vep_cache:/cache \
  ensemblorg/ensembl-vep:release_116.3 \
  vep -i /work/record.vcf -o /work/vep.txt --offline --cache \
  --dir_cache /cache --cache_version 116 \
  --species homo_sapiens --assembly GRCh38 --format vcf \
  --synonyms /work/synonyms.txt --no_stats --force_overwrite
```

This run writes the 2 rows of chromosome X. With `-e PERL_HASH_SEED=11` in place of `-e PERL_HASH_SEED=0` it writes the 13 rows of chromosome 21. With neither variable set, the rows are those of one chromosome or the other, by the seed Perl draws.

### 1.7 How the comparators treat this class

`vep-outputs-compare` sets the class's rows aside when it is given `--ensembl-cache <dir>`, the Ensembl VEP cache directory the reference run read. It takes the sequences that cache holds from its directories, as Ensembl VEP lists them (`CacheDir.pm` line 431), and reads its `chr_synonyms.txt`. It sets aside the reference's rows at a Location whose sequence name that cache does not hold and whose synonyms name two or more sequences it holds, where vep-rs wrote no row at that Location. It does so only for a reference output whose header states a command line that names no `--synonyms`, since only such a run read the cache's own file.

A record vep-rs wrote a row for stays a difference, as under `--dont_skip`, and so do rows on a sequence the cache holds. Without the option, under the 115.2 reference, and for a cache whose file gives no such name, nothing is set aside and no report changes. A reference run that was given `--synonyms`, as in the example above, keeps its rows open, and so does one whose output states no command line. The rule reads no warnings file and no flag of vep-rs's run. The format comparison, `vep-formats-compare`, holds no rule for the class.

`scripts/genomes/check_genome.sh` names the Ensembl cache it was given or extracted to the comparator, and a run given a converted cache alone names none. The reproduction `scripts/divergences/reproductions/synonym_of_several_sequences_by_hash_order.toml` holds the record, the synonyms file and both engines' rows, Ensembl VEP's under `PERL_HASH_SEED=0`.
