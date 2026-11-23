# A start site on or one base above a record's first base under `--nearest`

**The defect in one sentence.** Under `--nearest`, Ensembl VEP 116.3's search above a record begins two bases above the record's first base and never returns a protein-coding start site on the first base or one base above it. Such a site is named only where the search at or below the record's end returns it: the record must reach the site and cover no higher one. On any other record where the searches return a site but none nearer, Ensembl VEP names a farther site in its place, or only one of two sites at one distance. vep-rs searches from the first base as well, so the lowest site at or above it is a candidate on every record.

**Why this is a defect.** `NEAREST` names the protein-coding transcription start site nearest the record. No site is nearer to a record of one base than a site one base above it, unless a site lies on the record itself. Ensembl VEP names that site for a record two bases below it (1.6) and for records farther below it and above it (1.3). For the record one base below it alone, it names the nearest of the sites its two searches return, which can lie thousands of bases away; where neither search returns one, it returns an empty list (`TranscriptTree.pm`, line 239) and names no site.

A site on a record's first base is at distance zero, and no site is nearer. Ensembl VEP names such a site for a record of one base (2.6), and it names a site on a record's last base. Where the record covers a higher site as well, the site on the first base is not named, though no site Ensembl VEP names is nearer. The higher site lies on the record's last base in 2.3 and one base below it in the second record of 2.6.

**What vep-rs does instead.** vep-rs keeps Ensembl VEP's two candidates and adds a third: the lowest site at or above the record's start, which Ensembl VEP's search above the start passes over where it lies on the first base or one base above it. The nearest of the three is named, and every one of the three at that distance. The distances, the ties and the identifiers printed are Ensembl VEP's. A record with no site on its first base or one base above it gets the same sites in both engines.

Nearest means nearest among the candidates. As in Ensembl VEP, only the sites the searches return are candidates, however near an end of the record another site lies. In both engines a site below a record that covers a site is never a candidate. The lowest site above a record of two or more positions is a candidate exactly where the record covers no site two bases or more above its first base, the place Ensembl VEP's search above the start begins.

**Where it is in Ensembl VEP's source** (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `TranscriptTree.pm`, line 165: a start site at coordinate `t` is stored in the interval tree as the half-open interval from `t - 1` to `t`. Lines 172 and 204 make the same conversion where they ask the tree for overlaps.
- `TranscriptTree.pm`, lines 233-236: the two candidates are `fetch_nearest_down($e)` and `fetch_nearest_up($s)`, the record's end and start as given. The start is not converted.
- `Set::IntervalTree`, the library behind the tree, version 0.12 in the release's image. Its documentation says that all intervals are half-open and that `fetch_nearest_up` searches for the closest interval that does not contain the query. A site at `s + 1` is the interval from `s` to `s + 1`, which contains `s`. A site at `s` is the interval from `s - 1` to `s`, which ends at the query and lies below it.
- `TranscriptTree.pm`, lines 242-253: the nearer of the candidates found is kept, both on a tie.

**What correcting it would change.** Asked with `$s - 2`, the search above the start would return the lowest site at or above the start. The interval of a site at `s` starts at `s - 1` and that of a site at `s + 1` at `s`, both above that query. Ensembl VEP would then name the sites vep-rs names on every record of the three shapes below. Asked with `$s - 1`, the conversion lines 172 and 204 make, it would reach a site at `s + 1` and still pass over a site at `s`.

**The class it produces.** `nearest_start_site_one_base_above`: a `NEAREST` value on a row whose record has a protein-coding start site on its first base or one base above it. Ensembl VEP's value names the nearest of the sites its two searches return, and vep-rs's the nearest of those sites and of the lowest site at or above the first base. The two values agree where Ensembl VEP's search at or below the end returns that lowest site, or where one of its searches returns a nearer site. The format comparator accepts the two values under the class's name where each names its own set of sites and the two sets differ.

**How often the shape occurs.** Tier 2, in the terms of [the index](../intended-divergences.md#differences-outside-the-classes), which defines the tiers, the cuts, the corpora and the flag sets. The engines name different sites for three shapes of record, counted against the protein-coding start sites of the release 116 caches:

| Shape of the record | Records of the cuts |
| --- | --- |
| One base, one base below a site, no site on it (section 1) | 464 |
| Two or more positions, a site on the first base, a higher site on a later base (section 2) | 10 |
| Four or more positions, a site on the second base, none on the first or the last, a higher site between the second base and the last | 7 |

The two engines name different sites on every record of these shapes. Together 481 of the 15,646,481 records of the cuts carry a shape, 30.7 per million. Of the 511 committed flag sets, 3 pass `--nearest`, and 34 of the 152 corpora hold a record of one of the shapes.

In the third shape Ensembl VEP's search above the start begins on the record's third base, so both its questions find a higher site. vep-rs names the site on the second base, one base from the start, or that site and the highest one where the highest lies one base below the last base.

On one more shape the engines agree: a record of three or more positions that covers a site on its second base alone, with a site one base above its last base. Ensembl VEP's search above the start returns that outer site, and both sites are named, tied at one base. vep-rs keeps Ensembl VEP's candidate, so it names both too.

## 1. A substitution one base below the start site of KCNE2

**In short.** `21 35736322 . T C`, a 1000 Genomes Phase 3 record on GRCh37, lies one base below the start site of ENST00000290310, a protein-coding transcript of KCNE2. Under `--nearest symbol`, Ensembl VEP 116.3 writes `NEAREST=SMIM11`, a gene whose start site lies 11,477 bases above the record. The same row gives the record's distance from the KCNE2 transcript as 1. vep-rs writes `NEAREST=KCNE2`.

### 1.1 What a reader needs to know first

`--nearest transcript|gene|symbol` adds one value to every row of a record: the identifier of the protein-coding transcript whose start site lies nearest the record among the candidates ([Default and tab output](../output/vep-and-tab.md#nearest)). A transcript's start site is its start on the forward strand and its end on the reverse. A candidate is measured to the nearer end of the record. A candidate inside a longer record is measured the same way, so it is at distance zero only on an end. Two candidates at the smallest distance are both printed.

The value is the record's, whatever a row's own transcript is: every row of a record carries the same value, an intergenic row included.

### 1.2 The record

```
21	35736322	.	T	C	.	.	.
```

It is read against the release 116 GRCh37 cache with one flag:

```
--nearest symbol
```

The protein-coding start sites of the cache around the record:

| Site | Transcript | Symbol | From the record |
| --- | --- | --- | --- |
| 35445870 | ENST00000608209 | SLC5A3 | 290,452 below |
| 35736323 | ENST00000290310 | KCNE2 | 1 above |
| 35747799 | ENST00000399292, ENST00000399299 | SMIM11 | 11,477 above |

### 1.3 What each engine writes

The record has one row in both engines, on the KCNE2 transcript. Its Location, Allele, Gene, Feature and Consequence, then its Extra column, in Ensembl VEP 116.3:

```
21:35736322  C  ENSG00000159197  ENST00000290310  upstream_gene_variant
IMPACT=MODIFIER;DISTANCE=1;STRAND=1;NEAREST=SMIM11
```

vep-rs:

```
21:35736322  C  ENSG00000159197  ENST00000290310  upstream_gene_variant
IMPACT=MODIFIER;DISTANCE=1;STRAND=1;NEAREST=KCNE2
```

Every other column agrees. The dataset's records at 21:35736193 and 21:35736342, 130 bases below the site and 19 above it, print `NEAREST=KCNE2` in both engines.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP keeps the start sites in an interval tree whose intervals are half-open. It stores the KCNE2 site, 35736323, as the interval from 35736322 to 35736323. The interval holds the coordinate 35736322 and stops before 35736323.

For each record it asks the tree twice: for the nearest interval at or below the record's end, and for the nearest interval above the record's start. Both questions carry the record's own coordinate, 35736322. The first finds the SLC5A3 site, 290,452 bases below: the KCNE2 interval ends above the record, so it is no answer to that question. The second skips an interval that holds its query, so it passes over the KCNE2 site and finds the SMIM11 site.

The nearer of the two is kept: SMIM11, at 11,477 bases. Where Ensembl VEP asks the same tree for overlaps, it lowers the start by one first, the conversion it made when it stored the site. The second question is not lowered. A record two bases below the site asks with 35736321, which the KCNE2 interval does not hold, and gets KCNE2.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `TranscriptTree.pm`, line 165: the interval a site is stored as.
- `TranscriptTree.pm`, lines 233-236: the two questions.
- `TranscriptTree.pm`, lines 242-253: the choice between the two answers.

### 1.5 Why vep-rs's row is the correct one

The record is the one base at 35736322. ENST00000290310 lies on the forward strand and starts at 35736323, so its start site is one base above the record. Ensembl VEP's own row prints that distance, `DISTANCE=1`, on its upstream row for the transcript.

The SMIM11 site is at 35747799, 11,477 bases above the record, and the SLC5A3 site is 290,452 bases below it. No site lies on the record. One base is the smallest of the three distances, so the nearest protein-coding start site is KCNE2's, the value vep-rs writes.

### 1.6 Try it yourself

Two made-up records lie one base and two bases below the start site of ENST00000286808, a transcript of CLDN17. It lies on the reverse strand, so its start site is its end, 30166805. The cache of the `GRCh38-identifiers` golden corpus holds the transcript:

```
21	30166804	one_below_site	A	G	.	.	.
21	30166803	two_below_site	C	T	.	.	.
```

vep-rs's command:

```
vep -i records.vcf -o out.txt --offline --species homo_sapiens \
  --json_cache tests/golden/116/GRCh38-identifiers/json_cache \
  --assembly GRCh38 --nearest symbol
```

vep-rs writes `NEAREST=CLDN17` on the row of both records. Ensembl VEP 116.3 takes the same flag with the cache options of [Reproducing these examples](reproducing.md). Its cache directory must be writable: the first `--nearest` run on a cache writes its table of start sites there. It writes `NEAREST=CLDN17` for the second record and `NEAREST=CLDN8` for the first. The CLDN8 site, 30216097, is the next site above in the release 116 GRCh38 cache, 49,293 bases from the record.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/nearest_start_site_one_base_above.toml` runs the two records of 1.6 and the three of 2.6 on both engines under `--nearest transcript`. It records the `NEAREST` value of each record's row on the transcript its site belongs to. For the records of 1.6 Ensembl VEP's values are `ENST00000399899`, the CLDN8 transcript, and `ENST00000286808`. vep-rs's are `ENST00000286808` twice.

The 34 corpora that hold a record of one of the shapes run no flag set with `--nearest`, so no committed corpus documents a row under the class. The format comparator, `tools/vep-formats-compare`, reads the start sites of vep-rs's cache when it is given the cache with `--nearest-sites`. For a row's location it works out the sites each engine keeps: Ensembl VEP from its two questions, the upper one starting two bases above the first base, and vep-rs from those two and a third that starts on the first base.

An identifier names a site when a transcript with that id, gene or symbol starts there, and a value's identifiers are read in any order. The comparator accepts the two values under the class where the engines keep different sites, Ensembl VEP's identifiers name exactly the sites it keeps, and vep-rs's exactly the sites vep-rs keeps. Any other difference in `NEAREST` stays a difference, and so does every such pair without the cache.

## 2. An insertion between two start sites of TMPRSS2

**In short.** `21 41508115 . C CCCTCCG`, a 1000 Genomes high-coverage record on GRCh38, inserts six bases between two neighbouring positions. Each is the start site of a protein-coding transcript of TMPRSS2: ENST00000966232 below the insertion and ENST00000966231 above it. Under `--nearest transcript`, Ensembl VEP 116.3 writes `NEAREST=ENST00000966231`. vep-rs writes `NEAREST=ENST00000966231,ENST00000966232`.

### 2.1 What a reader needs to know first

An insertion lies between two positions, and both engines measure it over those two: Ensembl VEP gives it a start one above its end and swaps them for the search (`TranscriptTree.pm`, line 228). A site on either position is at distance zero. A deletion or a substitution of several bases is measured over its first and its last base in the same way. What 1.1 says of the value holds here too.

### 2.2 The record

```
21	41508115	.	C	CCCTCCG	.	.	.
```

It is read against the release 116 GRCh38 cache with one flag:

```
--nearest transcript
```

The protein-coding start sites of the cache around the record, all of TMPRSS2 on the reverse strand:

| Site | Transcript | From the record |
| --- | --- | --- |
| 41508105 | ENST00001063019 | 10 below its first position |
| 41508115 | ENST00000966232 | on its first position |
| 41508116 | ENST00000966231 | on its second position |
| 41508120 | ENST00001063018 | 4 above its second position |

### 2.3 What each engine writes

The record has 38 rows in both engines, each on a transcript of TMPRSS2. The rows on the two transcripts, their Location, Allele, Feature and Consequence, then their Extra column, in Ensembl VEP 116.3:

```
21:41508115-41508116  CCTCCG  ENST00000966231  5_prime_UTR_variant
IMPACT=MODIFIER;STRAND=-1;NEAREST=ENST00000966231
21:41508115-41508116  CCTCCG  ENST00000966232  upstream_gene_variant
IMPACT=MODIFIER;DISTANCE=0;STRAND=-1;NEAREST=ENST00000966231
```

vep-rs:

```
21:41508115-41508116  CCTCCG  ENST00000966231  5_prime_UTR_variant
IMPACT=MODIFIER;STRAND=-1;NEAREST=ENST00000966231,ENST00000966232
21:41508115-41508116  CCTCCG  ENST00000966232  upstream_gene_variant
IMPACT=MODIFIER;DISTANCE=0;STRAND=-1;NEAREST=ENST00000966231,ENST00000966232
```

Every other column agrees, and the 36 other rows of the record carry the same value in each engine.

### 2.4 Why Ensembl VEP writes that

The record's two positions are 41508115 and 41508116. The question at or below the end, 41508116, finds the site of ENST00000966231 on that position. The question above the start carries 41508115 and reaches no site below 41508117. The interval of the site at 41508116 holds the query. The interval of the site at 41508115 ends at it, so it lies below the query and is no answer to a question about what lies above. The question finds the site at 41508120.

The nearer of the two is kept: the site of ENST00000966231, at distance zero. The site of ENST00000966232, also at distance zero, was never a candidate. The question below the end reaches a site on the record's last position, and the question above the start cannot reach one on its first. A site on the first position is named only where the record covers no higher site, since the question below the end then finds it.

Source (`ensembl-vep` at release/116.3, `Bio/EnsEMBL/VEP/`):

- `TranscriptTree.pm`, line 165: the interval a site is stored as.
- `TranscriptTree.pm`, line 228: an insertion's start and end are swapped.
- `TranscriptTree.pm`, lines 233-236: the two questions.
- `TranscriptTree.pm`, lines 242-253: the choice between the two answers.

### 2.5 Why vep-rs's row is the correct one

The six bases are inserted between 41508115 and 41508116. ENST00000966232 lies on the reverse strand and ends at 41508115, so its start site is the base directly below the insertion. ENST00000966231 ends at 41508116, the base directly above it. Each site is on one of the record's two positions, and the record is symmetric about them: no measure makes one nearer than the other.

Ensembl VEP's own rows place the record the same way. Its upstream row for ENST00000966232 prints `DISTANCE=0`, and its row for ENST00000966231 lies in that transcript's 5' UTR. Where two sites it has found lie at one distance, Ensembl VEP prints both (`TranscriptTree.pm`, lines 250-253). Both sites lie at distance zero here, so the nearest protein-coding start sites are both, the value vep-rs writes. vep-rs prints the site found at or below the end first.

### 2.6 Try it yourself

Three made-up records lie at two start sites of SAMSN1 on the reverse strand: ENST00000959353's at 14546313 and ENST00000400566's at 14546315. The first deletes the bases from one site to the other. The second deletes one base more, so the higher site lies one base below its last base. The third is a substitution on the lower site. The cache of the `GRCh38-identifiers` golden corpus holds both transcripts:

```
21	14546312	first_base_and_last_base	AACT	A	.	.	.
21	14546312	first_base_and_inside	AACTG	A	.	.	.
21	14546313	on_the_site	A	C	.	.	.
```

Run with the command of 1.6 and `--nearest transcript`, vep-rs writes `NEAREST=ENST00000400566,ENST00000959353` on the rows of the first record. It writes `NEAREST=ENST00000959353` on the rows of the second and of the third. Ensembl VEP 116.3 writes `NEAREST=ENST00000400566` for the first two records and `NEAREST=ENST00000959353` for the third. Its question above the start of each deletion reaches no site below 14546315 and finds the higher site there.

### 2.7 How the comparator treats this class

The reproduction of 1.7 runs these three records beside the two of 1.6. For them Ensembl VEP's values are `ENST00000400566`, `ENST00000400566` and `ENST00000959353`. vep-rs's are `ENST00000400566,ENST00000959353`, `ENST00000959353` and `ENST00000959353`.

The comparator takes the rows of the first two records, and the rows of the record of 2.2, under the class as 1.7 describes. Ensembl VEP's two questions keep the higher site alone. vep-rs's third question finds the site on the first base, which is kept, and the higher site is kept beside it where it lies on the last base. The third record's rows agree.
