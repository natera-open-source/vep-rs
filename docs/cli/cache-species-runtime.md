# Cache, species and runtime

These flags name the annotation cache and the reference FASTA, the species and assembly, and the batch size, thread counts and distances of a run. The cache is found as Ensembl VEP finds it, under `--dir_cache`, `--species` and `--cache_version`, or named directly with `--json_cache`. Building or converting a cache is described under [Cache setup](../cache-setup.md).

## Cache and offline

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--cache` | String | none | Annotate from the local cache, which vep-rs always does; a directory value is the cache root, as `--dir_cache`, when that flag is absent |
| `--offline` | bool | false | Accepted without effect: vep-rs always runs offline and never contacts a database |
| `--dir` | String | `~/.vep` | Top-level VEP data directory: the cache root when neither `--dir_cache` nor a `--cache` directory is given |
| `--dir_cache` | String | none | Cache root: the cache read is \<dir\_cache\>/\<species\>[\_refseq\|\_merged]/\<cache\_version\>\_\<assembly\>, as Ensembl VEP lays it out (defaults to `--dir`) |
| `--cache_version` | u32 | none | Cache version: the prefix of the cache directory under the species directory (defaults to the VEP release) |
| `--json_cache`, `--full_cache_dir` | String | none | The cache directory itself, the one holding info.json and transcripts/ (built by vep-cache-builder or converted from an Ensembl VEP cache); it is read in place of the `--dir_cache` layout, and only `--assembly` is checked against its info.json |
| `--show_cache_info` | bool | false | Print the cache's source versions as sorted key\<TAB\>value lines and exit |
| `--fasta`, `--fa` | String | none | FASTA file for reference sequence lookup: plain with its .fai, or bgzipped with the .fai and .gzi that samtools faidx writes; a gzip file that is not BGZF is refused and no index is built. A directory stands for every indexed .fa and .fa.gz in it, read together. Without the flag, the first indexed .fa or .fa.gz by name in the cache directory is used |
| `--fasta_dir` | String | none | Directory holding FASTA files of several species or assemblies: the .fa or .fa.gz whose name carries `--species` and the assembly is used, a toplevel file before any other, a bgzipped file before a plain one and the first by name of two tied on both (where Ensembl VEP keeps the directory listing's order), in place of `--fasta` and of the cache directory's FASTA; when none qualifies the run warns and falls back to `--fasta` or the cache directory's FASTA |
| `--no_fasta` | bool | false | Do not use a FASTA file found in the cache directory |
| `--synonyms` | String | none | Chromosome synonyms file: on each line, whitespace-separated names of one chromosome, the first the reference name; an input chromosome the cache does not name is looked up through its synonyms. Replaces the cache's own chr\_synonyms.txt, which is read when this flag is absent |
| `--shift_3prime` | `0`, `1` | none | Move each insertion or deletion to its most 3' position on the transcript, over a 1,000-base window of the reference, before its consequence terms, positions, codons and amino acids are computed (1), or annotate it where the input placed it (0); needs a FASTA and stops the run without one. The HGVS notations shift under `--shift_hgvs`, not this flag. The shift runs along the genome only: a transcript whose cached sequence is edited away from the genome is not shifted again along its own sequence |
| `--shift_hgvs` | `0`, `1` | none | Move each insertion or deletion to its most 3' position on the transcript before its HGVSc and HGVSp are built (1, the default), or describe it where the input placed it (0), which also sets `--shift_3prime` to 0 and cannot be combined with `--shift_3prime` or `--shift_genomic` |
| `--shift_genomic` | `0`, `1` | none | Under `--shift_3prime` 1, print the shifted position in Location and the shifted sequence in Allele on every row (1), intergenic rows moving towards higher coordinates, or keep the input position (0); without `--shift_3prime` 1 nothing moves |
| `--shift_length` | bool | false | Add the SHIFT\_LENGTH column: under `--shift_3prime` 1, the bases each insertion or deletion moved along the transcript, 0 on every other transcript row of a sequence variant and none on a structural record's rows |
| `--sereal` | bool | false | Accepted without effect: the cache's serialiser is read from its own metadata |

## Species and assembly

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--species`, `-s` | String | `homo_sapiens` | Species: names the cache directory under `--dir_cache`; the Latin name only, with no alias such as human |
| `--assembly`, `-a` | String | none | Genome assembly (e.g., GRCh38): chooses among several caches of one version and must equal the cache's own name exactly; hg19 and hg38 are read as GRCh37 and GRCh38 |
| `--grch37` | bool | false | Annotate against GRCh37: sets `--assembly` GRCh37, over an explicit `--assembly`, as Ensembl VEP's option set does |

## Runtime

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--buffer_size` | usize | 5000 | Number of variants to process per batch |
| `--max_sv_size` | i64 | 10000000 | A record whose end minus start exceeds this many bases is oversize, a sequence record included; -1 lifts the limit, as in Ensembl VEP 116. With a FASTA, a \<CNV:TR\> record whose repeat sequence is longer than this stays a symbolic tandem repeat rather than expanding to its alleles (under -1 every one does, and 0 reads as 5000), as in Ensembl VEP 116. Ensembl VEP 116 drops an oversize record before annotation and writes nothing for it in any format; vep-rs annotates it: rows in the default, tab and parquet outputs, its VCF line without consequences, no JSON object (docs/output-formats.md and docs/intended-divergences.md describe the divergence). A region, SPDI or HGVS transcript record whose span exceeds the limit is skipped with Ensembl VEP's warning, as every oversize record is in Ensembl VEP 116; the skip is made as the line is parsed, so it holds under `--dont_skip` too, where Ensembl VEP keeps the feature (`--dont_skip` names the records it cannot keep). An HGVS genomic notation's end minus its start must be smaller than this value, 0 reading as 5000 and -1 as -1 itself, so under -1 a notation written in genome order is skipped, as Ensembl VEP 116 skips it |
| `--fork`, `--threads` | usize | 0 | Number of parallel annotation worker threads (Ensembl VEP forks processes; vep-rs runs threads and is parallel by default). Default 0 = auto: uses the number of available logical CPUs, capped at 32 (one reader thread feeds the workers, so extra workers wait at the batch barrier). Set to 1 to force single-threaded annotation, or any positive value to pin the worker count; user-supplied values are not capped |
| `--decompression_threads` | String | 0 | Number of BGZF decompression worker threads for a BGZF-compressed input (a plain gzip input is decompressed on one thread). Default 0: uses the number of available logical CPUs, capped at 10. `auto` chooses the count from the input's inflated bytes per record: a BGZF input with a .tbi or .csi index beside it whose records inflate to at least 8,000 bytes gets 32 workers on an aarch64 host, capped at the available logical CPUs; any other input, and every input on another architecture, gets the 0 rule. `auto` is opt-in, never the default. A positive value pins the worker count, e.g. `--decompression_threads=1` for single-threaded; a value above the number of available logical CPUs is lowered to that number |
| `--pipeline_mode` (hidden) | `batched`, `streaming` | `streaming` | How the batches move through the pipeline. `streaming` hands each batch to the worker pool as independent tasks that annotate and render their own records (and parse them, for a VCF input) and takes the next batch without waiting for them; the writer puts the rendered output back in input order. `batched` parses, annotates and renders each batch on the pool as a whole before the next one starts. A run with plugins, or one reading known variants from the cache's variation files or from the `--check_svs` file, uses `batched` whatever this says: plugins see a whole batch at once, and each of those sources prefetches one batch's regions at a time |
| `--pipeline_depth` (hidden) | u64 | 0 | Batches a pipeline stage may hold ahead of the next one, at most 1024; 0 is the built-in depth |
| `--render_task_records` (hidden) | u64 | 0 | Input records per pool task in both pipeline modes (per render task on the plugin path), at most 1,000,000; 0 is the built-in sizing |
| `--writer_strategy` (hidden) | `buffered`, `direct`, `vectored` | `buffered` | How rendered output reaches the file or standard output. `buffered` collects it in a buffer of `--output_buffer_size` bytes; `direct` writes each rendered chunk with one call; `vectored` writes the chunks that are ready together with one vectored call. Every strategy writes the same uncompressed bytes |
| `--output_buffer_size` (hidden) | u64 | 0 | Bytes of the `buffered` writer strategy's buffer, at most 1073741824 (1 GiB); 0 is the built-in size. The other strategies ignore it |
| `--resolved_config` (hidden) | bool | false | Print the resolved runtime knobs and the input's shape as name\<TAB\>value lines and exit without annotating |
| `--distance` | String (comma-separated) | 5000 | Upstream and downstream distance in bp: one integer, or the pair "up,down" (a downstream part written 0 reads as the upstream value) |

