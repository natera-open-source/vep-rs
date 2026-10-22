# vep-cache-builder

`vep-cache-builder` is the binary that builds the JSON (JavaScript Object Notation)
transcript cache the `vep` binary reads through `--json_cache`, from an Ensembl
release's GFF3 (General Feature Format version 3) and protein FASTA, without Ensembl VEP
(Variant Effect Predictor). No workspace crate depends on it; it uses
`vep-fasta` for genome reads and `vep-core` for the reverse complement, and `vep-cli`
loads what it writes. The output is `info.json`, one transcript shard per chromosome and
1 Mb region under `transcripts/` and, with `--include-variations`, one `all_vars.gz` with its
`.csi` per chromosome under `variations/`, the variation file Ensembl VEP's own cache carries,
then `CHECKSUMS.sha256`, one `sha256sum` line per file written, which `vep-cli` checks every
file it reads against.

## How it works

Follow one protein-coding transcript. `download::resolve_inputs` takes the GFF3 and
protein FASTA paths from `--gff3` and `--fasta`, or downloads them from the Ensembl
FTP (File Transfer Protocol) site into the system temporary directory under a
`vep_cache_builder_<species>_<assembly>_<release>` name that a later run reuses, and
hands both paths to `main`. `fasta::parse_protein_fasta` reads the FASTA once and maps
the unversioned transcript identifier in each header's `transcript:` token to its
peptide. `gff3::parse_gff3` reads the gene, transcript, exon and CDS (coding sequence)
lines into a `GeneRecord` holding `TranscriptRecord`s, each with its sorted `ExonRecord`s
and `CdsRecord`s. A gene is any line whose `ID` is `gene:<id>` and a transcript any line
whose `ID` is `transcript:<id>`, whatever the type column says (`gene`, `ncRNA_gene`,
`pseudogene`; `mRNA`, `lnc_RNA`, `miRNA` and the rest), so every gene and transcript
type reaches the cache. The reader strips the `chr` prefix and renames `chrM` to `MT`,
drops the `gene:` and `transcript:` identifier prefixes, and keeps what the Ensembl core
API writes on each line: the gene and transcript `version`, `ccdsid`,
`transcript_support_level`, the `tag` codes (`gencode_basic`, or `basic` on the GRCh37
release-87 file, `gencode_primary`, `Ensembl_canonical`, `MANE_Select`,
`MANE_Plus_Clinical`), the exon `ensembl_phase` and `ensembl_end_phase`, and the CDS
line's `protein_id` and translation `version`, which name the ENSP (Ensembl protein)
identifier. Two kinds of transcript never reach the cache, as Ensembl VEP's own dumper
drops them: a transcript of biotype `artifact` and a readthrough transcript (tagged
`readthrough_transcript` on the release 116 GRCh38 file, `readthrough_tra` by the core
API's attribute code).

Four optional inputs join here. `core_dumps::CoreDumps::load` reads the release's core
and variation MySQL table dumps from the `--core-dumps` directory (any table missing from
it is downloaded there first): the transcript attributes Ensembl VEP's dumper keeps (the
MANE accessions, APPRIS, TSL, the GENCODE labels, `cds_start_NF`, `cds_end_NF`, the
mature `miRNA` arms and the `ncRNA` structure), the CCDS and RefSeq xrefs of each
transcript, the UniProt, UniParc and UniProt isoform xrefs and the protein features of
each translation with their analysis display labels, the gene's display xref (the symbol,
its source database and, for HGNC, the HGNC identifier), the gene's canonical transcript,
the sequence regions' codon table, the `meta` source versions and, from the variation
database's `phenotype_feature` and `phenotype` tables, the count of significant gene
phenotype features Ensembl VEP's `--gene_phenotype` reads; the whole-genome GRCh38 dumps
(about 500 MB compressed, the phenotype features over half of it) read in about half a
minute and 3 GB of memory. `predictions::fetch_predictions`
hashes every peptide with MD5 (Message-Digest algorithm 5) and fetches the matching SIFT
and PolyPhen matrices from Ensembl's public MySQL server in batches of 500 hashes,
`gtf::parse_gtf_tags` reads the `#!` header keys and the per-transcript `tag` values of the
GTF (Gene Transfer Format), and `vep_fasta::IndexedFasta` opens `--genome-fasta`.
`builder::build_cache` then turns the transcript into a `CacheTranscript` the way Ensembl
VEP's own cache dumper reads the transcript off the core API. The exons are ordered 5' to
3' by `rank` and each carries the GFF3 line's `ensembl_phase` and `ensembl_end_phase`,
the exon's own Ensembl phases (an exon line without the pair gets them derived from the
CDS: -1 at an end in untranslated sequence, the CDS phase column with 1 and 2 swapped at
a coding 5' end, and `(phase + coding length) % 3` at a coding 3' end).
`mapper::build_mapper` lays the exons out into cDNA-to-genome pairs, locates
`cdna_coding_start` and `cdna_coding_end` from the CDS and records exon 1's phase as the
mapper's `start_phase`. With `--genome-fasta` the spliced cDNA is read exon by exon,
reverse-complemented on the minus strand, and cut into `five_prime_utr`,
`translateable_seq` (led by one `N` per base of the translation start exon's phase, so a
transcript whose coding sequence opens mid-codon keeps its frame) and `three_prime_utr`;
a transcript without a translation has the empty `translateable_seq`. `sorted_exons`
lists the exons by genomic start, the introns follow the exons in transcript order, the
translation's `start` and `end` are positions within its start and end exons, and
`codon_table` is the sequence region's attribute when the core dumps carry one, else 2
on `MT` and 1 elsewhere. The peptide comes from the FASTA map, the protein identifier and
translation version from the CDS line, `ccds` from `ccdsid` (or the dumper's `-`), and
the attributes from the core dumps when read and otherwise from the GFF3 `tag` list and
`transcript_support_level` (a MANE membership whose accession the GFF3 does not carry is
written with an empty value) plus the GTF's kept tags. With the core dumps the
transcript also carries the gene symbol's source and HGNC identifier, the gene phenotype
count, the RefSeq, SWISSPROT, TrEMBL, UniParc and UniProt isoform lists (`-` where there
are none), the protein features, the translation's sequence edits (`seq_edits`, the
`initial_met`, `_selenocysteine`, `amino_acid_sub` and `_stop_codon_rt` attributes as
`start`, `end` and `alt_seq`) and the sequence region's codon table; and, within each
1 Mb region, a transcript whose gene symbol matches that of a transcript carrying an HGNC
identifier takes that identifier, as Ensembl VEP's `merge_features` copies it at dump
time. The record carries the keys the dumper keeps and none of its own: no exon or
intron strand, no mapper pair count or map names. The builder bins the transcripts by chromosome and by the 1 Mb
region of their start and writes each bin as a pretty-printed JSON array to
`transcripts/<chr>/<start>-<end>.json`; `info::CacheInfoJson` writes `info.json`,
whose `cache_version` is `--release`, whose `sift` and `polyphen` read `b` only when a
matrix of that tool was written, whose `source_versions` come from the core `meta` table
(`assembly`, `gencode`, `genebuild`, as Ensembl VEP's cache info carries them) or, without
the dumps, from the GTF header, whose `contig_order` lists the chromosome directories
written (autosomes by number, X, Y, MT, then the rest by name) and whose `variation_cols`
and `var_type tabix` name the variation files' columns and layout only when the files were
written. With `--core-dumps`,
`info::write_chr_synonyms` writes `chr_synonyms.txt` beside it from the `seq_region` and
`seq_region_synonym` tables, the pairs Ensembl VEP's dump pipeline writes (each
sequence region's name against its synonyms, linked one step further through the
synonyms already recorded), so `NC_000021.9` or `CM000683.2` in an input resolves to the
cache's `21`. Last, `vep_core::cache_manifest::CacheManifest::write` lists every file of the
output directory in `CHECKSUMS.sha256`, the file `vep-cli` checks each cache file it reads
against.

Without `--core-dumps` the cache states what the GFF3 states and no more: no symbol
source or HGNC identifier, no gene phenotype count, no MANE accession, no APPRIS, no
mature miRNA arms, no RefSeq or UniProt lists, no protein features and no sequence edits, so the `SYMBOL_SOURCE`,
`HGNC_ID`, `GENE_PHENO`, `MANE_SELECT`, `APPRIS`, `miRNA`, `RefSeq`, `SWISSPROT`, `TREMBL`,
`UNIPARC`, `UNIPROT_ISOFORM` and `DOMAINS` columns stay empty on such a cache. The
Ensembl cache lists a translation's protein features in the database's row order, which
no dump carries; a native cache lists them by feature identifier, so `DOMAINS` may name
the same features in another order.

With `--include-variations`, `download::resolve_variation_inputs` scans
`--variation-vcf-dir` for the per-chromosome `*-chr*.vcf.gz` files and, under their Ensembl
names, the genome-wide `1000GENOMES-phase_3.vcf.gz`, `<species>_somatic.vcf.gz`,
`<species>_phenotype_associated.vcf.gz` and `<species>_clinically_associated.vcf.gz`, or
downloads them (Variant Call Format, each with its `.csi` index; a genome-wide file the
site does not publish, the 1000 Genomes file on GRCh37, is left out with a warning).
`variation::chromosome_records` reads one chromosome at a time: the per-chromosome file
whole and each genome-wide file through its index (or read through when it has none),
keeping each record whose `FILTER` is `PASS`, `.` or empty, adding the genome-wide
records the per-chromosome file lacks (a record is the same one when its name, position
and alleles are) with the somatic file's flagged `somatic`. Each record takes the form
the variation database holds, undoing what the export did to it: `end` is `POS` plus the
reference length minus one; a record whose class says the database holds an empty allele
(`TSA=deletion`, `TSA=insertion`) was exported behind an anchor base, so that base is
chopped from every allele, `start` moves one past `POS` and the emptied allele reads `-`,
and an insertion lies one base further right still (the export anchors it two bases before
the database start); a record marked `TSA=indel` has no empty allele in the database and
was anchored only if its alleles do not all share a first base, so it keeps the export's
alleles at `POS` (`A/AA`, and the dbSNP forms the export trimmed of their trailing bases,
`A/AGA` for the database's `AGA/AGAGA`) unless chopping leaves alleles that neither include
an empty one nor share a first base (`CC/CTTCT` to `C/TTCT`), the shape the export's anchor
leaves behind (a database record spelling one shared base before alleles that then differ
would be read as anchored; the sources' minimal spellings carry no such base); a COSMIC or
HGMD-PUBLIC variant
whose ALT is `.` or `<.>` reads `COSMIC_MUTATION` or `HGMD_MUTATION` on the span the
export's shape encodes, and every HGMD-PUBLIC variant reads `HGMD_MUTATION`. `clin_sig`
is the record's `CLIN_<term>` flags joined by `,` and `phenotype_or_disease` the
`E_Phenotype_or_Disease` evidence flag. The 1000 Genomes file's per-allele `AFR`, `AMR`,
`EAS`, `EUR` and `SAS` frequencies join each non-somatic record through the minimised
REF/ALT pair Ensembl VEP's own join matches on, printed as the dumper prints them
(`%.4g`), an allele at 0 in every population getting no entry; `AF`, which the file lacks,
is the allele's count over the panel's 5,008 chromosomes, each population's count recovered
from its frequency and sample size.
`clin_sig_allele`, `clin_sig_ref_allele` and `clinical_impact`, which Ensembl VEP joins
from the variation database's ClinVar phenotype features, are filled from the
`phenotype_feature` and `phenotype_feature_attrib` dumps of the `--variation-dumps`
directory (`variation::ClinVarAttribs`); no Ensembl VCF carries them.
`variation::write_chromosome` writes the chromosome's records sorted by start as
`variations/<chr>/all_vars.gz`, tab-delimited bgzip with the chromosome first and `.` for
an empty cell, and its `.csi` on the start column, the file and index the Ensembl VEP
dumper writes.

```text
<output-dir>/
  info.json
  chr_synonyms.txt      (with --core-dumps)
  transcripts/
    1/
      1-1000000.json
      1000001-2000000.json
  variations/
    1/
      all_vars.gz
      all_vars.gz.csi
```

The `variations/` tree exists only with `--include-variations`. A transcript shard
carries no `chr` field: the loader takes the chromosome from the directory name.

## Modules

| Module               | Role                                                                                                                                                                                         |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/main.rs`        | The `clap` argument struct and the run order: resolve inputs, resolve and read the core dumps, parse FASTA and GFF3, fetch predictions, load the genome FASTA and GTF, build the transcript cache, then the variation files, then `chr_synonyms.txt` and `info.json` |
| `src/download.rs`    | Ensembl FTP URL construction (GRCh37 gene set capped at release 87; the core and variation MySQL dump directories), downloads into the temporary directory with reuse, the `--core-dumps` directory fill, the local variation directory scan, the parallel per-chromosome VCF download and the genome-wide variation VCFs with their indexes |
| `src/fasta.rs`       | Protein FASTA reader: unversioned transcript identifier to peptide, plain or gzipped                                                                                                         |
| `src/gff3.rs`        | GFF3 reader into `GeneRecord`, `TranscriptRecord`, `ExonRecord` and `CdsRecord`; chromosome normalisation, attribute URL-decoding, the `tag` codes, `ccdsid`, `transcript_support_level`, the exon phases, the CDS `protein_id` and `version`, and the dumper's `artifact` and readthrough filters |
| `src/core_dumps.rs`  | The core and variation table dumps: MySQL escape decoding, the joins of `transcript_attrib`, `translation`, `protein_feature`, `xref`, `object_xref`, `gene`, `seq_region_attrib`, `meta`, `phenotype_feature` and `phenotype` into per-transcript and per-gene facts; the table reader `info.rs` uses for `seq_region_synonym` |
| `src/info.rs`        | `info.json` (species, assembly, version, the `sift` and `polyphen` markers, source versions, `variation_cols` and `var_type`, `contig_order`) and `chr_synonyms.txt` from the `seq_region_synonym` table |
| `src/gtf.rs`         | GTF `transcript` lines: the `#!` header keys and the five tags kept as transcript attributes (`basic` and `gencode_basic` both become `gencode_basic`)                                        |
| `src/mapper.rs`      | cDNA-to-genome mapper pairs, `cdna_coding_start` and `cdna_coding_end` for one transcript, with exon 1's phase as the start phase                                                              |
| `src/builder.rs`     | The serialised cache types, exon ordering and phases, introns, the spliced cDNA cut into UTRs and `translateable_seq`, the codon table, seq edits, attributes and xref lists from the GFF3 or the dumps, the per-region HGNC identifier copy, 1 Mb binning and the transcript shards |
| `src/predictions.rs` | SIFT and PolyPhen matrix fetch from Ensembl's MySQL server, keyed by the MD5 of the peptide and stored base64-encoded under the Ensembl matrix object's keys (`analysis`, `sub_analysis`, `translation_md5`) |
| `src/variation.rs`   | The variation VCF readers (per-chromosome whole, genome-wide through the index), each record in the variation database's form, the 1000 Genomes frequency join, the ClinVar columns from the variation dumps, and the `all_vars.gz` writer with its `.csi` |

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
| `--gtf <FILE>`                | Ensembl GTF of the same release: the `cds_start_NF` and `cds_end_NF` attributes and the `source_versions` keys that the GFF3 lacks, when no core dumps are read |
| `--core-dumps <DIR>`          | Directory of the release's core and variation MySQL table dumps (`<table>.txt.gz`, or `<table>.txt`); a table missing from it is downloaded there. Supplies every transcript attribute, xref, protein feature, symbol source, HGNC identifier, codon table, source version, chromosome synonym and gene phenotype count the GFF3 lacks |
| `--genome-fasta <FILE>`       | Indexed genome FASTA (`.fa` beside its `.fai`): `translateable_seq` for every coding transcript, the sequence behind codon-level consequences              |
| `--include-predictions`       | Fetch the SIFT and PolyPhen matrices; `homo_sapiens` on GRCh37 or GRCh38 only                                                                             |
| `--mysql-host <HOST>`         | Prediction server, default `ensembldb.ensembl.org`                                                                                                        |
| `--mysql-port <PORT>`         | Prediction server port, default `3306`                                                                                                                    |
| `--include-variations`        | Also write `variations/<chr>/all_vars.gz` and its `.csi` from Ensembl's variation VCFs: the per-chromosome files with the genome-wide somatic, phenotype-associated, clinically-associated and 1000 Genomes frequency files |
| `--variation-vcf-dir <DIR>`   | Local directory of the `*-chr*.vcf.gz` files and, under their Ensembl names, the genome-wide files with their `.csi` indexes, instead of the download        |
| `--variation-dumps <DIR>`     | Directory of the release's variation MySQL table dumps (`source`, `seq_region`, `attrib_type`, `phenotype`, `phenotype_feature`, `phenotype_feature_attrib`; nothing is downloaded): the ClinVar allele and somatic classification columns of the variation records, read with `--include-variations` for the chromosomes being built |
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
attribute decoding, the table dump joins on a hand-written slice, mapper arithmetic, GTF
tag selection, VCF parsing on temporary files and region binning. The integration test
`tests/native_cache_matches_ensembl_cache.rs` runs the binary on
`tests/fixtures/grch38_chr21/`, a release 116 slice of 27 chromosome 21 genes (the GFF3
lines, three GTF transcript lines, peptides and dump rows), and compares the identifier,
symbol, xref, attribute, protein feature and exon coordinate fields of the 63 transcripts
it writes with the same transcripts of the Ensembl VEP 116 GRCh38 cache. The integration test
`tests/variation_cache_from_dumps.rs` runs the binary with `--include-variations` and
`--variation-dumps` on indexed VCFs and dump rows written to a temporary directory, checks
the `all_vars.gz` lines and `info.json` it writes, and reads the records back through `vep`'s
variation source. On Linux the binary links OpenSSL (`reqwest` and `mysql` both use
`native-tls`), so the build needs the OpenSSL development package; macOS links
Security.framework. The workspace `rust-version` in `Cargo.toml` is set by the `noodles`
crates; this crate's `mysql` dependency chain is the next tier down.
`scripts/validation/validate_builder_cache.sh` builds the
binary, builds a one-region cache and compares it with captured fixtures
(`--fixture-mode`) or a Perl-derived JSON cache (`--baseline-cache`). The
[contributing guide](../../CONTRIBUTING.md) has the repository's standard build, test,
lint and format commands.

## Related

- [vep-fasta](../vep-fasta/README.md): `IndexedFasta`, the `.fa` plus `.fai` reader behind `--genome-fasta`.
- [vep-core](../vep-core/README.md): `codon::reverse_complement` for minus-strand CDS; `prediction.rs` decodes the matrices this crate embeds; `transcript.rs` holds the runtime `Transcript` the cache fills.
- [vep-cli](../vep-cli/README.md): `crates/vep-cli/src/json_cache.rs` loads the shards named by `--json_cache` and accepts both this crate's bare-array shard and the map shape a Perl-converted cache carries.
- [vep-cache-converter](../vep-cache-converter/README.md): the other route to a cache, validating a Perl-converted JSON cache and copying it into the runtime layout.
- [docs/cache-setup.md](../../docs/cache-setup.md): what each input decides, the GRCh37 command and variation data; [docs/cache-predictions.md](../../docs/cache-predictions.md): the SIFT and PolyPhen predictions.
- [docs/cli-reference.md](../../docs/cli-reference.md): the `vep` flags that read what this crate writes: `--json_cache`, `--sift`, `--polyphen`, `--check_existing`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg` and `--max_af`.
- [Root README](../../README.md#crates): the index of every crate.
