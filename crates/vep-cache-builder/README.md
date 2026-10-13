# vep-cache-builder

`vep-cache-builder` is the binary that builds the JSON (JavaScript Object Notation)
transcript cache the `vep` binary reads through `--json_cache`, from an Ensembl
release's GFF3 (General Feature Format version 3) and protein FASTA, without Perl VEP
(Variant Effect Predictor). No workspace crate depends on it; it uses
`vep-fasta` for genome reads and `vep-core` for the reverse complement, and `vep-cli`
loads what it writes. The output is `info.json`, one transcript shard per chromosome and
1 Mb region under `transcripts/` and, with `--include-variations`, matching shards under
`variations/`.

## How it works

Follow one protein-coding transcript. `download::resolve_inputs` takes the GFF3 and
protein FASTA paths from `--gff3` and `--fasta`, or downloads them from the Ensembl
FTP (File Transfer Protocol) site into the system temporary directory under a
`vep_cache_builder_<species>_<assembly>_<release>` name that a later run reuses, and
hands both paths to `main`. `fasta::parse_protein_fasta` reads the FASTA once and maps
the unversioned transcript identifier in each header's `transcript:` token to its
peptide. `gff3::parse_gff3` reads the gene, transcript, exon and CDS (coding sequence)
lines into a `GeneRecord` holding `TranscriptRecord`s, each with its sorted `ExonRecord`s
and `CdsRecord`s: it strips the `chr` prefix and renames `chrM` to `MT`, drops the
`gene:` and `transcript:` identifier prefixes, and reads the canonical, MANE (Matched
Annotation from NCBI and EMBL-EBI), TSL (transcript support level), APPRIS (Annotating
Principal Splice Isoforms) and GENCODE (encyclopedia of genes and gene variants) `basic`
tags plus the HGNC (HUGO Gene Nomenclature Committee) identifier from `Dbxref`. Three
optional inputs join here: `predictions::fetch_predictions` hashes every peptide with
MD5 (Message-Digest algorithm 5) and fetches the matching SIFT and PolyPhen matrices
from Ensembl's public MySQL server in batches of 500 hashes, `gtf::parse_gtf_tags` reads
the `#!` header keys and the per-transcript `tag` values of the GTF (Gene Transfer
Format), and `vep_fasta::IndexedFasta` opens `--genome-fasta`. `builder::build_cache`
then turns the transcript into a `CacheTranscript`: `mapper::build_mapper` lays the
exons out 5' to 3' into cDNA-to-genome pairs and locates `cdna_coding_start`,
`cdna_coding_end` and the start phase from the CDS; the exon phases come from CDS
overlap, the introns from the
gaps between exons, the peptide from the FASTA map, the protein identifier from the
transcript identifier with `ENST` replaced by `ENSP`, the attributes from the `basic`
tag, the GTF's kept tags, TSL and APPRIS, and `translateable_seq` from the CDS segments
read out of the genome FASTA in transcript order, reverse-complemented on the minus
strand. The builder bins the transcripts by chromosome and by the 1 Mb region of their
start and writes each bin as a pretty-printed JSON array to
`transcripts/<chr>/<start>-<end>.json`; `builder::generate_info_json` closes with
`info.json`, whose `cache_version` is `--release`, whose `sift` and `polyphen` read `b`
when predictions were fetched, and whose `source_versions` come from the GTF header.

With `--include-variations`, `download::resolve_variation_inputs` scans
`--variation-vcf-dir` for `*-chr*.vcf.gz` files or downloads Ensembl's per-chromosome
VCFs (Variant Call Format) and their optional `.csi` indexes in parallel;
`variation::parse_variation_vcf` keeps each record whose `FILTER` is `PASS`, `.` or
empty and that has an alternate allele, reading the first `ID`, the `REF/ALT` allele
string, `MA`, `MAF`, `CLIN_SIG`, `PUBMED` and every population frequency `INFO` field as
an `ALT:freq` string; `variation::build_variation_cache` writes the same 1 Mb bins as
compact JSON under `variations/<chr>/`.

```text
<output-dir>/
  info.json
  transcripts/
    1/
      1-1000000.json
      1000001-2000000.json
  variations/
    1/
      1-1000000.json
```

The `variations/` tree exists only with `--include-variations`. A transcript shard
carries no `chr` field: the loader takes the chromosome from the directory name.

## Modules

| Module               | Role                                                                                                                                                                                         |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/main.rs`        | The `clap` argument struct and the run order: resolve inputs, parse FASTA and GFF3, fetch predictions, load the GTF and genome FASTA, build the transcript cache, write `info.json`, then the variation cache |
| `src/download.rs`    | Ensembl FTP URL construction (GRCh37 gene set capped at release 87), downloads into the temporary directory with reuse, the local variation directory scan and the parallel per-chromosome VCF download |
| `src/fasta.rs`       | Protein FASTA reader: unversioned transcript identifier to peptide, plain or gzipped                                                                                                         |
| `src/gff3.rs`        | GFF3 reader into `GeneRecord`, `TranscriptRecord`, `ExonRecord` and `CdsRecord`; chromosome normalisation, attribute URL-decoding, HGNC, TSL, APPRIS and MANE extraction                     |
| `src/gtf.rs`         | GTF `transcript` lines: the `#!` header keys and the five tags kept as transcript attributes (`basic` and `gencode_basic` both become `gencode_basic`)                                        |
| `src/mapper.rs`      | cDNA-to-genome mapper pairs, `cdna_coding_start`, `cdna_coding_end` and the start phase for one transcript                                                                                    |
| `src/builder.rs`     | The serialised cache types, exon phases, introns, `translateable_seq`, attribute codes, 1 Mb binning, the transcript shards and `info.json`                                                    |
| `src/predictions.rs` | SIFT and PolyPhen matrix fetch from Ensembl's MySQL server, keyed by the MD5 of the peptide and stored base64-encoded                                                                         |
| `src/variation.rs`   | Variation VCF reader (bgzipped or plain), the frequency field mapping to `ALT:freq` strings and the variation shards                                                                          |

## Command line

The binary has no subcommands. `--release`, `--assembly` and `--output-dir` are required;
every other flag is optional. [docs/cli-reference.md](../../docs/cli-reference.md) documents
the `vep` binary, not this one, so the flags are listed here.

| Flag                          | Role                                                                                                                                                      |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--release <N>`               | Ensembl release number; names the download paths and becomes `cache_version` in `info.json`                                                              |
| `--assembly <NAME>`           | `GRCh37` or `GRCh38`; selects the `grch37/` FTP tree and the prediction database                                                                          |
| `--species <NAME>`            | Ensembl species directory name, default `homo_sapiens`                                                                                                    |
| `--output-dir <DIR>`          | The cache directory to write                                                                                                                              |
| `--gff3 <FILE>`               | Local GFF3, plain or gzipped, instead of the download                                                                                                      |
| `--fasta <FILE>`              | Local protein FASTA instead of the download                                                                                                               |
| `--gtf <FILE>`                | Ensembl GTF of the same release: the `cds_start_NF`, `cds_end_NF` and `gencode_primary` attributes and the `source_versions` keys that the GFF3 lacks      |
| `--genome-fasta <FILE>`       | Indexed genome FASTA (`.fa` beside its `.fai`): `translateable_seq` for every coding transcript, the sequence behind codon-level consequences              |
| `--include-predictions`       | Fetch the SIFT and PolyPhen matrices; `homo_sapiens` on GRCh37 or GRCh38 only                                                                             |
| `--mysql-host <HOST>`         | Prediction server, default `ensembldb.ensembl.org`                                                                                                        |
| `--mysql-port <PORT>`         | Prediction server port, default `3306`                                                                                                                    |
| `--include-variations`        | Also build `variations/` from Ensembl's variation VCFs                                                                                                    |
| `--variation-vcf-dir <DIR>`   | Local directory of `*-chr*.vcf.gz` files instead of the download                                                                                          |
| `--chromosomes <A,B,...>`     | Chromosomes whose variation VCFs to download; default 1 to 22, X, Y and MT. The transcript build always covers every chromosome the GFF3 carries          |

The GRCh38 build command is in the [root README's quick start](../../README.md#quick-start);
[docs/cache-setup.md](../../docs/cache-setup.md) explains what each input decides and gives
the GRCh37 command. Progress is logged to standard error.

## Run it

```bash
cargo build --release -p vep-cache-builder
target/release/vep-cache-builder --help
cargo test -p vep-cache-builder
```

The unit tests run in-process with no network and no MySQL server: URL construction, GFF3
attribute decoding, mapper arithmetic, GTF tag selection, VCF parsing on temporary files
and region binning. On Linux the binary links OpenSSL (`reqwest` and `mysql` both use
`native-tls`), so the build needs the OpenSSL development package; macOS links
Security.framework. The workspace `rust-version` in `Cargo.toml` is set by this crate's
`mysql` dependency chain. `scripts/validation/validate_builder_cache.sh` builds the
binary, builds a one-region cache and compares it with captured fixtures
(`--fixture-mode`) or a Perl-derived JSON cache (`--baseline-cache`). The
[contributing guide](../../CONTRIBUTING.md) has the repository's standard build, test,
lint and format commands.

## Related

- [vep-fasta](../vep-fasta/README.md): `IndexedFasta`, the `.fa` plus `.fai` reader behind `--genome-fasta`.
- [vep-core](../vep-core/README.md): `codon::reverse_complement` for minus-strand CDS; `prediction.rs` decodes the matrices this crate embeds; `transcript.rs` holds the runtime `Transcript` the cache fills.
- [vep-cli](../vep-cli/README.md): `crates/vep-cli/src/json_cache.rs` loads the shards named by `--json_cache` and accepts both this crate's bare-array shard and the map shape a Perl-converted cache carries.
- [vep-cache-converter](../vep-cache-converter/README.md): the other route to a cache, validating a Perl-converted JSON cache and copying it into the runtime layout.
- [docs/cache-setup.md](../../docs/cache-setup.md): what each input decides, the GRCh37 command, predictions and variation data.
- [docs/cli-reference.md](../../docs/cli-reference.md): the `vep` flags that read what this crate writes: `--json_cache`, `--sift`, `--polyphen`, `--check_existing`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg` and `--max_af`.
- [Root README](../../README.md#crates): the index of every crate.
