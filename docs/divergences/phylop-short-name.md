# The key of a `--phyloP` source

**The defect in one sentence.** Ensembl VEP's `--phyloP <N>` option expands into a `--custom` source whose short name is `phlyoP<N>way`, two letters transposed, so every header line, column and key the source adds is misspelt while the file it reads is named `phyloP<N>way`. vep-rs gives the same source the short name `phyloP<N>way`, and every value agrees.

**Why this is a defect.** The key names the UCSC track the source reads, phyloP, the phylogenetic p-value score, and the file's own name spells it so: `hg38.phyloP100way.bw`. The sibling `--phastCons <N>` option spells its key `phastCons<N>way`, as its file is spelt. The transposition occurs once in Ensembl VEP's source, in the option set that writes the key, and in the test that pins it, and nowhere else. A user who filters on the column must type the misspelling, and a header that names the file spells the track two ways on one line.

**What vep-rs does instead.** vep-rs expands `--phyloP <N>` into the same source with the short name `phyloP<N>way`, the name `--custom <file>,phyloP<N>way,bigwig,exact` would give it. Every value agrees; the key alone differs.

**Where it is in Ensembl VEP's source** (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/Config.pm`, line 555: the `phyloP` entry of the option set writes `%ucsc_data_root%%ucsc_assembly%/phyloP%phyloP%way/%ucsc_assembly%.phyloP%phyloP%way.bw,phlyoP%phyloP%way,bigwig,exact`, the path spelt `phyloP` twice and the short name `phlyoP`.
- `t/Config.t`, lines 99-100: the module's test asserts the expansion with `phlyoP7way` and `phlyoP100way` as the short names.

**What correcting it would change.** The class vanishes: both engines print `phyloP<N>way`, and a `filter_vep` expression written against the misspelt key needs the corrected one there, as it does on vep-rs's output.

**The class it produces.** `conservation_short_name_spelling`, applied by the golden harness and by the format comparator to the reference text before the comparison. The concordance figures score consequence sets and never see the key, so nothing is set aside in them.

Tuple, pair, set aside, charged and the other scoring terms are defined in [How to read the counts](../intended-divergences.md#how-to-read-the-counts).

## 1. The misspelt `phyloP` key

**In short.** `--phyloP 100` adds one bigWig source to a run. Ensembl VEP names it `phlyoP100way`, from the option set's own text, and that key heads its header line, its tab column, its `Extra` item, its CSQ field, its `##INFO` line and its JSON member. vep-rs names it `phyloP100way`, the spelling of the track, of the file the source reads and of every other place Ensembl VEP writes the word. The values under the two keys are identical.

### 1.1 What a reader needs to know first

`--phyloP <N>` and `--phastCons <N>` are shortcuts for a `--custom` bigWig source of `type=exact`. Each expands into one source whose file lies in a UCSC tree, `<root>/<assembly>/phyloP<N>way/<assembly>.phyloP<N>way.bw`, where `--ucsc_assembly` names `<assembly>` and `--ucsc_data_root` names `<root>`. vep-rs reads the tree from a local directory alone, as [bigWig sources and conservation shortcuts](../output/custom-bigwig.md) describes.

A custom source's short name is its key in every format. In the default format it heads the `## <key> : <file>` header line and the `<key>=<value>` item of `Extra`; in the tab format it names a column. In the VCF it names a field of the CSQ list and its `##INFO=<ID=<key>` line, and in the JSON it names a member of the record's `custom_annotations` object. A score the source holds at the variant's position is the value under the key.

### 1.2 The record

The record is `21 5028999 . A G` from `variants.vcf` of the `GRCh38-custom` golden corpus (`tests/golden/116/GRCh38-custom/`). The corpus's `custom/ucsc/` tree holds a synthetic `hg38.phyloP100way.bw` and `hg38.phastCons100way.bw`, written by `scripts/golden/make_custom_fixtures.py` with one one-base item per record of the corpus:

```
21	5028999	.	A	G	.	.	.
```

Both engines run the `phylop-phastcons` flag set of the corpus's manifest:

```
--phyloP 100 --phastCons 100 --ucsc_assembly hg38 --ucsc_data_root <corpus>/custom/ucsc/
```

### 1.3 What each engine writes

The two header lines the shortcuts add, then the record's row, in the default format. Ensembl VEP:

```
## phlyoP100way : [PATH]/hg38.phyloP100way.bw
## phastCons100way : [PATH]/hg38.phastCons100way.bw
21_5028999_A/G	21:5028999	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;phastCons100way=0.998000025749207;phlyoP100way=-0.980000019073486
```

vep-rs:

```
## phyloP100way : [PATH]/hg38.phyloP100way.bw
## phastCons100way : [PATH]/hg38.phastCons100way.bw
21_5028999_A/G	21:5028999	G	-	-	-	intergenic_variant	-	-	-	-	-	-	IMPACT=MODIFIER;phastCons100way=0.998000025749207;phyloP100way=-0.980000019073486
```

The two rows agree in every value. The key differs in the header line and in the `Extra` item, and the same key difference appears in the tab format's column line, in the VCF's `##INFO` line and CSQ list, and in the JSON record:

```
##INFO=<ID=phlyoP100way,Number=.,Type=String,Description="[PATH]/hg38.phyloP100way.bw">
"custom_annotations":{"phastCons100way":[{"name":0.998000025749207}],"phlyoP100way":[{"name":-0.980000019073486}]}
```

against vep-rs's `##INFO=<ID=phyloP100way,...>` and `"phyloP100way":[{"name":-0.980000019073486}]`. The `phastCons100way` key agrees between the engines everywhere.

### 1.4 Why Ensembl VEP writes that

Ensembl VEP turns `--phyloP <N>` into a `--custom` string through a table of option sets. The `phyloP` entry of that table builds the source's file path from `%ucsc_data_root%`, `%ucsc_assembly%` and the width, spelling `phyloP` in both path segments, and then writes the short name as `phlyoP%phyloP%way`. The short name becomes the source's key, so the key is misspelt wherever the source prints it, while the path on the same line is not.

The module's own test asserts this expansion with the misspelt names, so the text is pinned rather than incidental, and the sibling `phastCons` entry on the next lines spells its short name as its path.

Source (`ensembl-vep`):

- `Bio/EnsEMBL/VEP/Config.pm`, lines 552-557: the `phyloP` option set.
- `Bio/EnsEMBL/VEP/Config.pm`, lines 559-564: the `phastCons` option set, spelt consistently.
- `t/Config.t`, lines 99-100: the expansions asserted with `phlyoP7way` and `phlyoP100way`.

### 1.5 Why vep-rs's row is the correct one

The key names what the source holds, the UCSC phyloP track, and the file the option itself resolves to is spelt `phyloP<N>way` in both its directory and its name. The option is `--phyloP`. Ensembl VEP spells the word `phyloP` in the option's name, in the path it builds and in every message about it; the short name is the one place it does not. The `--phastCons` sibling, built by the same table, spells its key as its file.

A filter expression is where the spelling matters. On Ensembl VEP's output the field is `phlyoP100way`, so `filter_vep -f "phlyoP100way lt -0.9"` selects the rows scoring below -0.9; on vep-rs's output the field is `phyloP100way`, so an expression written against the misspelt column must use the corrected name there. [`vep filter`](../cli/filter.md) resolves a field by its exact name, its case variants or a unique prefix of one key, and `phlyoP100way` is none of these for `phyloP100way`.

### 1.6 Try it yourself

Run the record of 1.2 through both engines, the corpus's `custom/ucsc/` tree as the data root. vep-rs's command is below; Ensembl VEP takes the same three shortcut flags with the cache and FASTA options of [Reproducing these examples](reproducing.md):

```
vep -i record.vcf -o out.txt --offline --json_cache <corpus>/json_cache \
  --fasta <corpus>/reference.fa.gz --species homo_sapiens --assembly GRCh38 \
  --phyloP 100 --ucsc_assembly hg38 --ucsc_data_root <corpus>/custom/ucsc/
```

| Output | Ensembl VEP | vep-rs |
| --- | --- | --- |
| Header line | `## phlyoP100way : [PATH]/hg38.phyloP100way.bw` | `## phyloP100way : [PATH]/hg38.phyloP100way.bw` |
| `Extra` item | `phlyoP100way=-0.980000019073486` | `phyloP100way=-0.980000019073486` |
| Tab column | `phlyoP100way` | `phyloP100way` |
| JSON member | `"phlyoP100way":[{"name":-0.980000019073486}]` | `"phyloP100way":[{"name":-0.980000019073486}]` |

Naming the same file as a plain `--custom` source with the short name `phyloP100way` makes Ensembl VEP print the corrected key with the same values, which is the source vep-rs's shortcut builds:

```
--custom file=<corpus>/custom/ucsc/hg38/phyloP100way/hg38.phyloP100way.bw,short_name=phyloP100way,format=bigwig,type=exact
```

### 1.7 How the comparator treats this class

The golden harness (`crates/vep-cli/tests/common/mod.rs`, `masked_reference_text`) and the format comparator (`scripts/concordance/compare_vep_formats.py`, `reference_text_masked`) apply the class `conservation_short_name_spelling` when a corpus flag set or a comparison case names it in `masks`. Under the class every `phlyoP<N>way` in the reference text, a digit run between the prefix and `way`, is read as `phyloP<N>way` before the header lines are compared and the rows are parsed. The comparator counts each rewritten occurrence under the class in its report. Nothing else changes: a value under the key compares as any other value, and the file name, spelt correctly on both sides, is untouched.

A comparison case's floor on the source's column is read under the renamed key whichever spelling names it, so the default, tab and VCF formats count the annotated rows on both sides. The JSON format carries the source once per record under `custom_annotations`, so the floor is not applicable there and the text formats of the case check it.

The `phylop-phastcons` flag set of the `GRCh38-custom` corpus names the class, and so does every comparison case that runs `--phyloP`. The class applies to no other key, and a run that does not name it compares a `phlyoP<N>way` key as written.
