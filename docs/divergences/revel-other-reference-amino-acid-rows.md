# The REVEL score of a table row of another reference amino acid

**The defect in one sentence.** Ensembl VEP's REVEL plugin reads each table row's reference amino acid and tests nothing with it, so a missense row takes the score of a table row that names another reference amino acid, the score of a substitution the row does not report. vep-rs writes a table row's score only where both of its amino-acid cells equal the output row's, and no score where the table holds no such row.

**Why this is a defect.** The table's publisher describes a row's two amino-acid columns as "aaref: Reference amino acid" and "aaalt: Alternate amino acid". It describes the transcript column as the "Ensemble transcript ID(s) corresponding to the given amino acid substitution", included "to indicate which transcript(s) the amino acid substitution applies to" ([REVEL downloads](https://sites.google.com/site/revelgenomics/downloads)). A row therefore scores one substitution, from its reference residue to its alternate residue, on the transcripts it lists. Where the transcript model reads another reference residue at the codon, the output row's `Amino_acids` reports another substitution than the table row scores. In the example below the codon is `cTc`, leucine, and the table row that lists the transcript scores serine to proline.

**What vep-rs does instead.** vep-rs writes a table row's score only where both of its amino-acid cells equal the output row's: the reference residue and the alternate residue of `Amino_acids`, compared as written. The published table spells its residues with the twenty one-letter codes alone, as the output row does, so no other sign is read as the same residue. A table row of another reference residue is passed over and the walk goes on in file order.

The output row takes the score of the first table row that names both of its residues and lists its transcript. Under `no_match` it takes the first table row that names both residues, whatever transcripts that row lists. It takes no score where the table holds no such row.

The plugin's parameters and fields are described under [REVEL](../plugins/annotation-store.md#revel).

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116):

- `REVEL.pm`, lines 210 and 229, `parse_data`: a row's cells are assigned to a list of variables, `$refaa` among them. The hash the function returns holds the alleles, the positions, `altaa`, the transcript ids and the score; `$refaa` is given no key, and no other line names it.
- `REVEL.pm`, line 193, `run`: for a row whose transcript ids are filled, the score is returned where one id is the transcript's stable id and `$_->{altaa} eq $tva->peptide`, the alternate amino acid alone. Line 197 makes the same test without the transcript, for a row with no ids and for every row under `no_match`.
- `REVEL.pm`, lines 170-199: the table rows of the variant's span are walked in file order, and the first row that passes is the one returned.

**What correcting it would change.** A missense row would carry a REVEL score only from a table row that scores its own substitution on its transcript. The rows that lose a score stand on transcripts the table lists under another reading of the coding sequence than the cache's transcript model. On most of them the model pads an incomplete coding start with the first exon's phase and the table reads the transcript without the pad, so the table's rows for the transcript name the residues of other codons.

**The class it produces.** `revel_other_reference_amino_acid_rows`: a `REVEL` field (the JSON `revel` member) where Ensembl VEP's value is the score of the table row its walk reaches, a row whose reference amino acid is not the first residue of the row's `Amino_acids`. vep-rs writes nothing there, or the score of the first table row of the row's own allele that names both residues and, where transcripts are matched, lists the row's transcript. The class reads the table the run's plugin read and sets nothing aside without it, under the test of section 1.7.

A comparison run without the class reports a disagreement of this shape as a difference between the engines and counts it against vep-rs. The concordance figures compare the engines' consequence terms and read no plugin field, so this difference does not move them.

**How often the shape occurs.** The index, [Intended divergences from Ensembl VEP](../intended-divergences.md), places the entry in tier 2 of its three tiers: ordinary data carries the shape, and the difference needs a plugin. The index counts how often a shape occurs over five sets of chromosome 21 records, of gnomAD, of 1000 Genomes and of ClinVar, 15,646,481 records in all. Of them 273 (17.4 per million) give a missense row whose transcript the table lists under another reference residue beside the row's alternate residue. None of the committed test corpora under `tests/golden/116/` holds such a record. Of the corpora's committed runs, which the index calls flag sets, 6 load the REVEL plugin.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A ClinVar substitution on a transcript whose coding start is incomplete

**In short.** The ClinVar record `21 43059244 471355 A G` on the reverse-strand transcript `ENST00000430013`, against the published GRCh38 REVEL table's six rows at the position. The transcript model reads the codon `cTc`, leucine, changed to `cCc`, proline: `Amino_acids` is `L/P`. The table's `A G` row that lists `ENST00000430013` scores serine to proline, `0.033`; its `A G` row for leucine to proline scores `0.419` and lists the gene's complete transcripts. Ensembl VEP 116.3 writes `REVEL=0.033` on the `ENST00000430013` row; vep-rs writes no `REVEL` there. Both engines write `REVEL=0.419` on the rows of the four listed transcripts the cache holds.

### 1.1 What a reader needs to know first

The published REVEL table has one row per position, alternate base and amino-acid substitution. Each row names its reference amino acid, its alternate amino acid, its score and the transcripts the substitution applies to. The plugin matches a table row to an output row by the row's transcript and its alternate amino acid.

`ENST00000430013` carries the `cds_start_NF` flag: its coding start is incomplete, and its sequence begins inside a codon. Ensembl VEP pads such a coding sequence by the first exon's phase, here one base, so its codons fall where they do on the gene's complete transcripts. The table lists the transcript under a reading without the pad, one base out of step, so its rows for the transcript name the residues of other codons.

### 1.2 The record

ClinVar's record 471355 on GRCh38, a one-base substitution:

```
21	43059244	471355	A	G	.	.	.
```

The published table's rows at the position, its third column the GRCh38 position:

```
21	44479354	43059244	A	C	L	R	0.395	ENST00000398158;ENST00000398165;ENST00000359624;ENST00000352178;ENST00000398168;ENST00000539520;ENST00000544202
21	44479354	43059244	A	C	S	A	0.032	ENST00000430013
21	44479354	43059244	A	G	L	P	0.419	ENST00000398158;ENST00000398165;ENST00000359624;ENST00000352178;ENST00000398168;ENST00000539520;ENST00000544202
21	44479354	43059244	A	G	S	P	0.033	ENST00000430013
21	44479354	43059244	A	T	L	H	0.361	ENST00000398158;ENST00000398165;ENST00000359624;ENST00000352178;ENST00000398168;ENST00000539520;ENST00000544202
21	44479354	43059244	A	T	S	T	0.044	ENST00000430013
```

The cache of the `GRCh38-plugin-revel` golden corpus (`tests/golden/116/GRCh38-plugin-revel/json_cache`) holds the transcripts. The corpus's table slice does not hold the position, so the reproduction writes its rows into a table of its own (1.6). Both engines run `--plugin REVEL,file=<table>`.

### 1.3 What each engine writes

The rows on `ENST00000430013` and on `ENST00000352178`, one of the listed complete transcripts, in the default format. Ensembl VEP 116.3:

```
471355	21:43059244	G	ENSG00000160200	ENST00000352178	Transcript	missense_variant	2002	1205	402	L/P	cTc/cCc	-	IMPACT=MODERATE;STRAND=-1;REVEL=0.419
471355	21:43059244	G	ENSG00000160200	ENST00000430013	Transcript	missense_variant	166	167	56	L/P	cTc/cCc	-	IMPACT=MODERATE;STRAND=-1;FLAGS=cds_start_NF;REVEL=0.033
```

vep-rs:

```
471355	21:43059244	G	ENSG00000160200	ENST00000352178	Transcript	missense_variant	2002	1205	402	L/P	cTc/cCc	-	IMPACT=MODERATE;STRAND=-1;REVEL=0.419
471355	21:43059244	G	ENSG00000160200	ENST00000430013	Transcript	missense_variant	166	167	56	L/P	cTc/cCc	-	IMPACT=MODERATE;STRAND=-1;FLAGS=cds_start_NF
```

The `REVEL` field by transcript:

| Feature | `Amino_acids` | Ensembl VEP 116.3 | vep-rs |
| --- | --- | --- | --- |
| `ENST00000352178` | `L/P` | `0.419` | `0.419` |
| `ENST00000430013` | `L/P` | `0.033` | none |

The rows on `ENST00000359624`, `ENST00000398158` and `ENST00000398165` carry `REVEL=0.419` in both engines, as the `ENST00000352178` row does. The record's other missense rows stand on transcripts no table row lists and carry no `REVEL` in either engine.

### 1.4 Why Ensembl VEP writes that

On the `ENST00000430013` row the loop at line 170 of `REVEL.pm` reads the table rows at the position in file order. The third row matches the record's alleles and names proline, the row's alternate amino acid, but lists other transcripts, so line 193 passes it over. The fourth row matches the alleles, lists `ENST00000430013` and names proline, so line 193 returns its score, `0.033`. Its reference amino acid, serine, was read into `$refaa` at line 210 and is tested by no line.

### 1.5 Why vep-rs's row is the correct one

The output row reports leucine to proline: its codon is `cTc`, the codon the gene's complete transcripts read at the position. The table scores that substitution `0.419`, in a row that does not list `ENST00000430013`. The row that lists the transcript scores serine to proline, a substitution the row does not report, so `0.033` beside `L/P` states a score the table does not give for the row's change. The table row for leucine to proline does not apply to the transcript by the table's own transcript list. vep-rs therefore writes no score on the row, and keeps `0.419` on the transcripts the leucine row lists.

### 1.6 Try it yourself

Write the two `A G` rows of 1.2 into a table, each with one of its transcripts, then annotate the record with each engine:

```
printf '%s\n' '#chr hg19_pos grch38_pos ref alt aaref aaalt REVEL Ensembl_transcriptid' \
  '21 44479354 43059244 A G L P 0.419 ENST00000352178' \
  '21 44479354 43059244 A G S P 0.033 ENST00000430013' \
  | tr ' ' '\t' | bgzip -c > revel_rows.tsv.gz && tabix -s 1 -b 3 -e 3 revel_rows.tsv.gz
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --assembly GRCh38 --plugin REVEL,file=revel_rows.tsv.gz
```

Ensembl VEP's `ENST00000430013` row carries `REVEL=0.033` and vep-rs's carries no `REVEL`; the `ENST00000352178` rows agree on `0.419`. With `no_match=1` appended to the plugin's parameters, both engines read the first row whose alternate amino acid is proline, the leucine row, and agree on `0.419` on every missense row.

### 1.7 How the comparator treats this class

The format comparator (`tools/vep-formats-compare`) applies `revel_other_reference_amino_acid_rows` when a comparison case names it and `--revel-table` hands it the table the run's plugin read, with `--revel-no-match` for a run that set `no_match`. Without the table the class sets nothing aside. A `REVEL` field is set aside, counted once per row, where both outputs spell the row's `Amino_acids` alike and Ensembl VEP's value is the score of the table row its walk reaches. That is the first row, in file order, that matches an alternate of the record, lists the row's transcript where transcripts are matched, and names the row's alternate amino acid. It must name another reference amino acid than the row's.

vep-rs's field must be what its own walk gives: the score of the first table row of the row's own allele that names both residues and, where transcripts are matched, lists the row's transcript, or nothing where there is none.

To find the row Ensembl VEP's walk reaches, the comparator needs every alternate allele of the row's input record. In the VCF and JSON formats it reads them from the output's own record: the `REF` and `ALT` columns of the row's line, or the `allele_string` of the row's object. In the default and tab formats it reads them from the row's `Uploaded_variation` column, the record's name, where the record has no identifier and both engines name it by its position and alleles. For a record with an identifier it reads them, in those two formats, from the input file handed to it by `--input`, from the one line that carries the identifier.

Where none of these gives the record's alleles, a row at one position with a one-base allele is still read. Ensembl VEP's walk can then have reached the first passing table row of the row's own allele or, before it in the file, the first passing row of another allele. The field is set aside only where every one of those rows that carries Ensembl VEP's score names another reference amino acid.

Every other difference in the `REVEL` field stays a difference. That covers an Ensembl VEP score that no such table row explains, an Ensembl VEP score from a table row that names the output row's own reference and alternate amino acids, and any vep-rs value that its own walk does not give.

A value that [`revel_sibling_allele_rows`](../plugins/annotation-store.md#9-revel-the-score-of-a-sibling-allele-on-a-multi-allelic-record) also accepts is counted under that class, which the comparator applies before this one. The golden harness, the test of vep-rs against the committed corpora, compares the `REVEL` fields of the two REVEL corpora, `GRCh37-plugin-revel` and `GRCh38-plugin-revel`, verbatim: no committed row has the shape, so the class applies to none of their rows.
