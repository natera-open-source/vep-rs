# Coming from Ensembl VEP

This page is for a reader who runs Ensembl VEP today and wants to know what changes when the same command line runs on vep-rs.

## What stays the same

vep-rs takes the same input files and flag names as Ensembl VEP and writes the same output layouts. The consequence terms are Ensembl VEP's own Sequence Ontology vocabulary. The VEP default, `--vcf`, `--json` and `--tab` outputs match Ensembl VEP's layouts. vep-rs adds `--parquet`.

`--plugin <Name>,<params>` and `--dir_plugins` keep their roles. A plugin name resolves to one of the built-in plugins first, and otherwise to a dynamic library in `--dir_plugins`. `--fasta` takes an uncompressed, `faidx`-indexed genome and gives `--hgvs` the reference context for the 3' shifting of indels.

## What differs

**The transcript cache.** `--json_cache` names a directory built by `vep-cache-builder` or converted from an Ensembl VEP cache ([Cache setup](../cache-setup.md)). `--cache`, `--offline`, `--dir_cache`, `--cache_version` and `--species` are accepted but do not locate a cache, and vep-rs never contacts a database. `--database`, `--host`, `--port`, `--user` and `--password` do not exist and are rejected as unknown flags. The transcript set is the cache's contents, so `--refseq`, `--merged` and `--gencode_basic` select nothing.

**Flags accepted without effect, and flags refused.** Some Ensembl VEP flags are accepted but change nothing, among them `--regulatory`, `--custom`, `--coding_only`, `--pick_order`, `--domains`, `--mirna` and the statistics flags, so no summary file and no regulatory, custom or domain annotation is written. Three flags, `--regulatory_gff`, `--extended_promoters` and `--custom_suppress_filter`, are refused with an error that says what each would do, and [Compatibility with Ensembl VEP](../cli/compatibility.md) lists every flag in both groups.

**Consequence terms and input formats.** The six regulatory-feature consequence terms are never assigned, because no regulatory-feature cache is loaded. `--format` takes `vcf`, `ensembl`, `hgvs`, `region` and `guess`. HGVS input is genomic (`:g.`) only, and the `id`, `spdi` and `caid` formats are not read.

**Where the outputs differ.** Each case where the two engines' outputs are known to differ is traced to the mechanism in Ensembl VEP's source and reproducible from the record. The [Intended divergences](../intended-divergences.md) pages list those in the consequence terms; the HGVSp of a shifted insertion is under [`--hgvs`](../cli/output-control.md#opt_hgvs).

## The commands to swap in

Three commands cover the move. The first runs once per release and assembly; the other two are your existing command before and after.

**1. Build the cache.** Build the transcript cache, or convert the Ensembl VEP cache you already have with `scripts/data/storable_to_json.pl`. The downloads of the GTF and the genome FASTA are in the [Quick start](quick-start.md#2-get-a-transcript-cache).

```bash
vep-cache-builder \
  --species homo_sapiens --assembly GRCh38 --release 116 \
  --gtf Homo_sapiens.GRCh38.116.gtf.gz \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --output-dir $HOME/.vep/json_cache/homo_sapiens/116_GRCh38
```

**2. Your Ensembl VEP command.** As it runs today, for example:

```bash
vep --cache --offline --dir_cache $HOME/.vep --assembly GRCh38 \
  -i input.vcf -o output.txt
```

**3. The same command on vep-rs.** With `--json_cache` added; the other flags stay and are accepted:

```bash
vep --cache --offline --dir_cache $HOME/.vep --assembly GRCh38 \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  -i input.vcf -o output.txt
```

`--json_cache` names the directory holding `info.json` and `transcripts/`, so `116_GRCh38` itself and not its parent; [Cache setup](../cache-setup.md) covers the GRCh37 command, variation data and the SIFT and PolyPhen matrices.
