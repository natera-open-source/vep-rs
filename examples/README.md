# Examples

Two small human VCFs to annotate before any download: each is 175 records cut from the golden corpus of its assembly under `tests/golden/116/`, so the pruned cache that ships beside that corpus covers every record and the run needs no cache build.

| File | Assembly | Records | Cache |
| --- | --- | --- | --- |
| `homo_sapiens_GRCh37.vcf` | GRCh37 | 175 | `tests/golden/116/GRCh37/json_cache` |
| `homo_sapiens_GRCh38.vcf` | GRCh38 | 175 | `tests/golden/116/GRCh38/json_cache` |

The records are the corpora's own: public ClinVar, gnomAD and 1000 Genomes variants with the corpus record's identifier in the `ID` column, each with explicit `REF` and `ALT` alleles of at most 50 bases together, the `QUAL`, `FILTER` and `INFO` columns empty. Each file carries single-nucleotide variants, insertions and deletions over the autosomes and X (23 chromosomes in the GRCh37 file, 22 in the GRCh38 file). A record whose span exceeds `--max_sv_size`, a symbolic allele and a tandem repeat are left to the corpora, so every example record annotates with no warning.

From the repository root, after `cargo build --release`:

```bash
target/release/vep \
  --json_cache tests/golden/116/GRCh37/json_cache \
  --assembly GRCh37 \
  -i examples/homo_sapiens_GRCh37.vcf -o output.txt
```

The same command with `GRCh38` in both places annotates the other file. Each record's rows are the rows the corpus run writes for it: `crates/vep-cli/tests/examples.rs` compares the two on every `cargo test`, and the release and container smoke tests annotate `homo_sapiens_GRCh37.vcf` beside the corpus.

A cache for real data, and the flags the quick start continues with, are on [Quick start](../docs/getting-started/quick-start.md). Examples of the other species the golden corpora cover are on [Other species](../docs/cache-species.md).
