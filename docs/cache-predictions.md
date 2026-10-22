# SIFT and PolyPhen predictions

`--sift` and `--polyphen` take `p` (prediction), `s` (score) or `b` (both), reading matrices embedded in the cache. `vep-cache-builder --include-predictions` fetches per-transcript prediction matrix blobs from Ensembl's public MySQL server and stores them base64-encoded and gzip-compressed beside each transcript. The server is `ensembldb.ensembl.org`, anonymous read-only; port 3306 serves the GRCh38 databases and `--mysql-port 3337` the GRCh37 ones. The database is `homo_sapiens_variation_<release>_<assembly>`, the tables `translation_md5`, `protein_function_predictions` and `attrib`, keyed on the MD5 of the full peptide.

At runtime the decoder in `crates/vep-core/src/prediction.rs` looks up `(protein_position, alt_amino_acid)` and prints the prediction as Ensembl VEP does: `deleterious(0.01)`, `probably_damaging(0.998)`. The score is `prob/1000`, so `0.01` rather than `0.010`; under `s` the bare score is printed, under `p` the bare prediction.

A prediction is looked up for every single-residue substitution inside the coding sequence, the shape Ensembl VEP predicts on: a `missense_variant`, a `start_lost` whose codon still translates (`M/T`), a same-length substitution of several bases within one codon. It is never looked up for a synonymous change (whose peptide allele string is the one residue, `C`, not `C/C`), a stop gain or loss, an unknown alternate residue, an insertion or a deletion.

PolyPhen reads the HumVar matrix, as Ensembl VEP does by default; `--humdiv` reads the HumDiv matrix instead, and a cache without that matrix prints an empty column, with no error in either engine.

The tools are checked against the cache's metadata before any record is read, as in Ensembl VEP. A cache whose `info.json` names neither a `source_versions` entry for the tool nor the tool's marker stops a run that asks for it, with Ensembl VEP's text: `ERROR: [cache] SIFT not available` or `ERROR: [cache] PolyPhen not available`. Ensembl VEP reads its `info.txt`'s `source_sift` and `source_polyphen` lines alone (`AnnotationSource/Cache/Transcript.pm` line 155), so a converted cache keeps a marker only beside its version line and the two engines stop on the same caches. Under `--everything` or `--safe` the run instead drops the tool with `INFO: disabling SIFT` on standard error (withheld under `--quiet`) and writes no column for it.

A natively built cache names no tool version; its `sift` and `polyphen` markers, written only when the build embedded a matrix of that tool, count as availability instead. A native cache built without `--include-predictions` is therefore the one case the flags refuse that Ensembl VEP's own caches never present.

```bash
vep-cache-builder --release 116 --assembly GRCh38 --output-dir /path/to/cache \
  --gff3 ... --gtf ... --include-predictions \
  --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GRCh37: --assembly GRCh37 --mysql-port 3337

vep -i input.vcf --sift b --polyphen b --json_cache /path/to/cache --offline
```

A converted cache carries the SIFT, PolyPhen HumVar and PolyPhen HumDiv matrices the Storable holds for each translated transcript. They are base64-encoded under the same `predictions_data` key and in the same binary layout the native builder writes, so `--sift` and `--polyphen` read either kind of cache alike. `info.json` keeps the `sift` and `polyphen` markers only when the conversion wrote a matrix of that tool and `info.txt` named the tool's `source_<tool>` version. A conversion that kept a marker without its version line runs the tool here where Ensembl VEP refuses it; converting the cache again drops the marker. The `dbNSFP` plugin is a separate, unsupported route to third-party SIFT and PolyPhen scores ([Plugins](plugins.md)).
