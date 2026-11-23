# A gene-model file that writes `chrM`, read beside a cache

Every file here is made up. `chrm_gene_synthetic.gtf.gz` is a GTF of one
gene, `MITOG1`, with one transcript, `MITOT1`, of one exon at positions 101 to
200 of the forward strand, of the biotype `Mt_tRNA`. Its three lines write the
sequence `chrM`. The file is sorted, compressed with `bgzip` and indexed with
`tabix -p gff`. It copies no published gene model.

The file is read with `--gtf` beside a GRCh37 cache whose synonyms link `MT`
and `chrM`, over the made-up reference of 16,569 bases in
`../mitochondrial-frame-check/reference_synthetic.fa`. For a record on `MT`
inside the exon vep-rs writes the row of `MITOT1` beside the cache's rows.
Ensembl VEP 116.3 writes the cache's rows and no row of `MITOT1`, and writes
that row when the same three lines name the sequence `MT`.

`crates/vep-cli/tests/gene_model_chrm_beside_a_cache.rs` runs the record of
`scripts/divergences/reproductions/file_beside_cache_chrm_transcripts_dropped.toml`
on these files, and
[the entry's page](../../../../../docs/divergences/file-beside-cache-chrm-transcripts-dropped.md)
explains the difference.
