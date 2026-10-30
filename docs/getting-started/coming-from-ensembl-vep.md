# Coming from Ensembl VEP

This page is for a reader who runs Ensembl VEP and wants to know what changes when the same command line runs on vep-rs.

## What stays the same

vep-rs takes the same input files and flag names as Ensembl VEP and writes the same output layouts. The consequence terms are Ensembl VEP's own Sequence Ontology vocabulary. The VEP default, `--vcf`, `--json` and `--tab` outputs match Ensembl VEP's layouts. vep-rs adds `--output_format parquet`.

`--plugin <Name>,<params>` and `--dir_plugins` keep their roles. A plugin name resolves to one of the built-in plugins first, and otherwise to a dynamic library in `--dir_plugins`. `--fasta` takes a `samtools faidx`-indexed genome, plain or compressed with `bgzip` (Ensembl's `dna_index` files), and gives `--hgvs` the reference context for the 3' shifting of indels; `--hgvs` without a FASTA is refused, as Ensembl VEP refuses it offline.

## What differs

**The transcript cache.** The cache is a directory built by `vep-cache-builder` or converted from an Ensembl VEP cache ([Cache setup](../cache-setup.md)). It is found under `--dir_cache`, `--species` and `--cache_version` in Ensembl VEP's `<dir_cache>/<species>/<version>_<assembly>` layout, or named directly with `--json_cache`; a run with no cache stops with an error, and vep-rs never contacts a database. `--cache` and `--offline` are accepted: vep-rs always annotates from the local cache, and a directory given to `--cache` is the cache root when `--dir_cache` is absent. `--refseq` and `--merged` select the `_refseq` or `_merged` species directory. `--refseq` also resolves an HGVS gene name among RefSeq transcripts first, under `--json_cache` too, where `--merged` changes nothing. `--gencode_basic` drops every transcript without the `gencode_basic` attribute as the cache loads.

**Flags accepted without effect, and flags refused.** `--custom` and the statistics flags are accepted without effect: no custom annotation or summary file is written. Database flags are refused, as are regulatory GFF sources and extended promoters; [Compatibility with Ensembl VEP](../cli/compatibility.md) lists every flag. `--regulatory`, alone or through `--everything`, annotates regulatory and motif features from a converted Ensembl VEP cache ([Regulatory and motif features](../regulatory.md)).

**Consequence terms and input formats.** The six regulatory and motif consequence terms need `--regulatory` and a cache converted with its regulatory features. `--format` takes `ensembl`, `vcf`, `hgvs`, `id`, `spdi`, `region` and `guess`; `guess` detects the first data line. HGVS genomic, transcript and protein input resolves against `--fasta` and the transcript cache; region and SPDI input also need `--fasta`. Identifier lists resolve against cache variation files and are refused without them; CAid input is refused, as in offline Ensembl VEP.

**Where the outputs differ.** Each case where the two engines' outputs are known to differ is traced to the mechanism in Ensembl VEP's source and reproducible from the record; the [Intended divergences](../intended-divergences.md) pages list them all.

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

**2. Your Ensembl VEP command.** As it runs, for example:

```bash
vep --cache --offline --dir_cache $HOME/.vep --assembly GRCh38 \
  -i input.vcf -o output.txt
```

**3. The same command on vep-rs.** A cache built or converted into Ensembl VEP's own layout, `$HOME/.vep/homo_sapiens/116_GRCh38`, is found by `--dir_cache`, `--species` and `--cache_version`, and the command runs unchanged. A cache elsewhere, such as the one step 1 builds, is named with `--json_cache`; the other flags stay and are accepted:

```bash
vep --cache --offline --dir_cache $HOME/.vep --assembly GRCh38 \
  --json_cache $HOME/.vep/json_cache/homo_sapiens/116_GRCh38 \
  -i input.vcf -o output.txt
```

`--json_cache` names the directory holding `info.json` and `transcripts/`, so `116_GRCh38` itself and not its parent; [Cache setup](../cache-setup.md) covers the GRCh37 command and variation data, and [SIFT and PolyPhen predictions](../cache-predictions.md) covers the matrices.
