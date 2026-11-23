# A MaxEntScan matrix line that is no number

**The defect in one sentence.** Ensembl VEP's MaxEntScan plugin multiplies a matrix line that is no number as whatever number Perl converts its opening characters to, `0.5x` as 0.5 and an empty line as 0, and scores every row from it. vep-rs refuses the matrix file when the plugin loads, naming the file, the line and its text, and the run goes on without the plugin.

**Why this is a defect.** The plugin multiplies every line it reads and takes the logarithm of the product, so a line stands for a number. A line that is not one still yields a score, printed with three decimals like every other, and nothing in the output tells it from a score over a number. Perl itself marks the conversion. The module runs under `use warnings`, and each multiplication of such a line raises `Argument "0.5x" isn't numeric in multiplication (*)`. That warning names the module's own line 860, never the matrix file or the line in it.

**What vep-rs does instead.** vep-rs reads each line as the number Perl's arithmetic makes of its text wherever Perl makes one without that warning. Those are a whole number, a decimal, an exponent form, blanks around any of them, the ten bytes `0 but true`, and an infinity or a not-a-number in the spellings Perl takes whole. Any other line is refused as the plugin loads, `<file> line <n> is not a number: '<text>'`. The run goes on without the plugin, as it does for any plugin that fails to load, and every row is written without a `MaxEntScan_` or `MES-` field.

Perl raises the warning for some texts in one place of a product and not in another. `nanx` as the first factor beside a decimal is converted in silence, and beside a whole number, or as any later factor, under the warning. vep-rs refuses such a line wherever it stands.

**Where it is in Ensembl VEP's source** (`MaxEntScan.pm` of the `VEP_plugins` repository at release/116):

- Line 75: `use warnings`.
- Lines 754-759, `score5_makescorematrix`, and 832-837, `score3_makemaxentscores`: each line is kept as text after `chomp` and one `s/\s//`, with no test of what it holds.
- Line 689, `$a * $d`, and line 860, the product of nine entries: the kept text is multiplied, so Perl converts the number that opens it, 0 where none does, and warns.

**What correcting it would change.** A test of each kept line as the tables are read would refuse the file at start-up, where the module already refuses a directory that does not exist (line 100) and a file it cannot open (lines 737, 753 and 831). Every line that is a number would be read as before.

**The class it produces.** `maxentscan_matrix_line_not_a_number`: the rows of a run whose matrix directory holds a line Perl converts only under its `isn't numeric` warning. Ensembl VEP writes the plugin's fields from the converted number, and vep-rs writes none. The class is documented only: no committed run reads such a line, so no manifest names it and no comparator applies it.

**How often the shape occurs.** Tier 3: the shape is a line of a matrix file. Of the 511 committed flag sets, 6 pass the plugin, each over the made-up matrices `scripts/fixtures/maxentscan_fixture_matrices.py` writes, whose every line is a number.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. A line that reads `0.5x`

**In short.** Line 1,024 of `splicemodels/me2x3acc1` in a made-up matrix directory reads `0.5x`. The record `21 14219563 G A` lies in the acceptor window of the first intron of ENST00000400577, whose score reads that line. Ensembl VEP 116.3 writes `MaxEntScan_ref=-0.860`, the score with the line taken as 0.5, and two `isn't numeric` lines in its warnings file. vep-rs writes the row without a `MaxEntScan_` field and one warning that names the file and the line.

### 1.1 What a reader needs to know first

The acceptor model scores a 23-mer, the last 20 bases of an intron and the first 3 of the exon after it. Its score multiplies nine entries, one from each of the files `me2x3acc1` to `me2x3acc9` (`MaxEntScan.pm` lines 851-860). A file holds one entry a line, and the bases of the window choose the line. A substitution in the window is scored twice, on the reference sequence and on the changed one, as `MaxEntScan_ref` and `MaxEntScan_alt`.

The directory here is made up. `scripts/fixtures/maxentscan_fixture_matrices.py --seed 1 --plain-lines` writes it, every number drawn from the seed, and one line of it is then rewritten by hand.

### 1.2 The record

One record of the `tests/golden/116/GRCh38-plugin-maxentscan-fixture` corpus, the first exonic base after intron 1 of ENST00000400577, a forward-strand transcript:

```
21	14219563	acceptor_exon_base	G	A	.	.	.
```

Both engines run `--fasta <fasta> --plugin MaxEntScan,<synthetic_matrices>` with the corpus's cache and reference. Line 1,024 of `splicemodels/me2x3acc1`, which the reference window and the changed window both read, holds `0.5x`.

### 1.3 What each engine writes

Ensembl VEP 116.3's row on ENST00000400577:

```
acceptor_exon_base	21:14219563	A	ENSG00000185272	ENST00000400577	Transcript	missense_variant,splice_region_variant	127	97	33	A/T	Gcg/Acg	-	IMPACT=MODERATE;STRAND=1;MaxEntScan_alt=-2.629;MaxEntScan_diff=1.769;MaxEntScan_ref=-0.860
```

Its warnings file, the module's path shortened; the number in the first line differs from run to run:

```
WARNING: 12 : Argument "0.5x" isn't numeric in multiplication (*) at MaxEntScan.pm line 860, <__ANONIO__> line 3.
Argument "0.5x" isn't numeric in multiplication (*) at MaxEntScan.pm line 860, <__ANONIO__> line 3.
```

vep-rs's row and the line of its warnings file:

```
acceptor_exon_base	21:14219563	A	ENSG00000185272	ENST00000400577	Transcript	missense_variant,splice_region_variant	127	97	33	A/T	Gcg/Acg	-	IMPACT=MODERATE;STRAND=1
```

```
WARNING: Failed to instantiate plugin MaxEntScan: <synthetic_matrices>/splicemodels/me2x3acc1 line 1024 is not a number: '0.5x'
```

Both engines exit 0. With the line as the generator wrote it, both write `MaxEntScan_ref=-1.563`.

### 1.4 Why Ensembl VEP writes that

The module keeps the line as the text `0.5x`. At line 860 the text is the first factor of the product, and Perl's conversion reads the longest number that opens it, 0.5, and raises its warning for the rest. The product, the logarithm and the three decimals follow as for any entry, so the row carries a score of a table whose line 1,024 holds 0.5. Ensembl VEP copies the Perl warning into its warnings file, once as raised and once under its own `WARNING:` prefix.

Source (`MaxEntScan.pm` at release/116):

- Lines 832-837: the line is kept as text.
- Line 860: the product that converts it.
- Line 245: the score printed with three decimals.

### 1.5 Why vep-rs's row is the correct one

The file states no number on that line. `0.5x` is not 0.5: the text may be a number cut short or joined to the next, or no number at all. A score needs nine numbers, and with one of them unread no score follows from the file. `-0.860` is the score of another table, one the user never supplied.

The conversion is also not one rule. An empty line converts to 0, on which the logarithm of line 796 or the division of line 860 dies and the row loses every field. `high` converts to 0 too, and `1e5.5` to 100,000. vep-rs treats each alike: the line is named and no score is written from it.

### 1.6 Try it yourself

```
python3 scripts/fixtures/maxentscan_fixture_matrices.py --seed 1 \
  --plain-lines --out synthetic_matrices
awk 'NR == 1024 { $0 = "0.5x" } { print }' \
  synthetic_matrices/splicemodels/me2x3acc1 > me2x3acc1
mv me2x3acc1 synthetic_matrices/splicemodels/me2x3acc1
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --species homo_sapiens --assembly GRCh38 \
  --fasta <corpus>/reference.fa.gz \
  --plugin MaxEntScan,synthetic_matrices
```

vep-rs writes the row without a score and the warning in `out.txt_warnings.txt`. Ensembl VEP with the same directory writes `MaxEntScan_ref=-0.860`; the commands of [Reproducing these examples](reproducing.md#the-commands) give the engines and their caches.

### 1.7 How the comparator treats this class

The reproduction `scripts/divergences/reproductions/maxentscan_matrix_line_not_a_number.toml` writes the directory, rewrites the line and holds Ensembl VEP's `MaxEntScan_ref` on the row, `-0.860`, beside vep-rs's empty field. No committed flag set reads a matrix line that is no number, so no manifest names the class and the class is documented only.
