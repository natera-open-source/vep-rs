# SIFT and PolyPhen predictions

`--sift` and `--polyphen` take `p` (prediction), `s` (score) or `b` (both), reading matrices embedded in the cache. `vep-cache-builder --auto` reads per-transcript prediction matrix blobs from the release's variation table dumps, and `--include-predictions` fetches them from Ensembl's public MySQL server. The builder stores them base64-encoded and gzip-compressed beside each transcript.

## Build the cache and run

```bash
vep-cache-builder --release 116 --assembly GRCh38 --output-dir /path/to/cache \
  --gff3 ... --gtf ... --include-predictions --codon-table MT=2 \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GRCh37: --assembly GRCh37 --mysql-port 3337
# --codon-table MT=2: the mitochondrion's genetic code (or pass --core-dumps)

vep -i input.vcf --sift b --polyphen b --json_cache /path/to/cache --offline
```

Each `...` stands for a file of the release: its Ensembl GFF3 after `--gff3`, its Ensembl GTF after `--gtf` ([Build a cache](cache-setup.md#build-a-cache)). `input.vcf` is the user's own file.

## What a run prints

At runtime the decoder in `crates/vep-core/src/prediction.rs` looks up `(protein_position, alt_amino_acid)` and prints the prediction as Ensembl VEP does: `deleterious(0.01)`, `probably_damaging(0.998)`. The score is `prob/1000`, so `0.01` rather than `0.010`; under `s` the bare score is printed, under `p` the bare prediction.

PolyPhen reads the HumVar matrix, as Ensembl VEP does by default. `--humdiv` reads the HumDiv matrix instead, and a cache without that matrix prints an empty column, with no error in either engine.

## Which changes get a prediction

A prediction is looked up for every single-residue substitution inside the coding sequence, the shape Ensembl VEP predicts on:

- a `missense_variant`
- a `start_lost` whose codon still translates (`M/T`)
- a same-length substitution of several bases within one codon

It is never looked up for:

- a synonymous change (whose peptide allele string is the one residue, `C`, not `C/C`)
- a stop gain or loss
- an unknown alternate residue
- an insertion
- a deletion

## When a tool is refused

The tools are checked against the cache's metadata before any record is read, as in Ensembl VEP. A cache whose `info.json` names neither a `source_versions` entry for the tool nor the tool's marker stops a run that asks for it. The line is Ensembl VEP's text: `ERROR: [cache] SIFT not available` or `ERROR: [cache] PolyPhen not available`.

Ensembl VEP reads its `info.txt`'s `source_sift` and `source_polyphen` lines alone (`AnnotationSource/Cache/Transcript.pm` line 155). So a converted cache keeps a marker only beside its version line, and the two engines stop on the same caches.

Under `--everything` or `--safe` the run instead drops the tool and writes no column for it. It prints `INFO: disabling SIFT` on standard error (withheld under `--quiet`).

A natively built cache names no tool version. Its `sift` and `polyphen` markers, written only when the build embedded a matrix of that tool, count as availability instead. A native cache of a species with predictions whose build read no matrix is therefore the one case the flags refuse that Ensembl VEP's own caches never present: a build with `--output-dir` and no `--include-predictions`, or a GRCh37 installation without the switch.

## How the matrices are fetched and stored

An `--auto` installation reads the matrices from the table dumps the release publishes of the variation database: `attrib_type`, `attrib`, `translation_md5` and `protein_function_predictions`. It reads them for any species whose database holds matrices, and checks each file against its `CHECKSUMS` line like every other input.

The GRCh37 tree publishes the matrices table without the `translation_md5` table that keys it. So a GRCh37 installation holds no matrices, and its closing warning names `--include-predictions`. Beside `--auto` that switch fetches the matrices from Ensembl's MySQL server instead, for human on GRCh37 (`--mysql-port 3337`) and GRCh38 alone.

With `--include-predictions`, which serves human on GRCh37 and GRCh38 alone, the build fetches the matrices from:

- the server `ensembldb.ensembl.org`, anonymous read-only
- port 3306 for the GRCh38 databases and `--mysql-port 3337` for the GRCh37 ones
- the database `homo_sapiens_variation_<release>_<n>`, `<n>` being `38` for GRCh38 and `37` for GRCh37
- the tables `translation_md5`, `protein_function_predictions`, `attrib` and `attrib_type`, keyed on the MD5 of the full peptide

A matrix stays gzip-compressed as the cache stores it and is inflated at the first lookup against its transcript. So a run without `--sift` or `--polyphen` inflates none.

## A converted cache

A converted cache carries the SIFT, PolyPhen HumVar and PolyPhen HumDiv matrices the Storable holds for each translated transcript. They are base64-encoded under the same `predictions_data` key and in the same binary layout the native builder writes. So `--sift` and `--polyphen` read either kind of cache alike.

`info.json` keeps the `sift` and `polyphen` markers only when the conversion wrote a matrix of that tool and `info.txt` named the tool's `source_<tool>` version. A conversion that kept a marker without its version line runs the tool here where Ensembl VEP refuses it. Converting the cache again drops the marker.

The `dbNSFP` plugin is a separate route to third-party predictions, read from the columns a user names in their own file ([Plugins](plugins.md)).
