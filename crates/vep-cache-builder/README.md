# vep-cache-builder

`vep-cache-builder` is the binary that builds the JSON (JavaScript Object Notation)
transcript cache the `vep` binary reads through `--json_cache` or `--dir_cache`, from an
Ensembl release's GFF3 (General Feature Format version 3) and protein FASTA, or from any
GFF3 or GTF (Gene Transfer Format) and a genome FASTA, without Ensembl VEP (Variant Effect
Predictor). With `--auto <species> <assembly>` it installs a release in one command: every
input fetched from the Ensembl site and checked against its directory's `CHECKSUMS`, the
cache built in Ensembl VEP's layout and the indexed genome FASTA placed inside it. Such an
installation brings the known variants with their PubMed identifiers, synonyms and failed
flag, the SIFT and PolyPhen matrices and the regulatory and motif features with their
cell-type activity and matrix counts, wherever the release publishes a file of them for the
genome, and ends by naming each part and cell it left out. `vep-cli`
depends on the crate's library, whose `gene_model` module reads a gene-model file a run
names beside its cache; the crate uses `vep-fasta` for genome reads and `vep-core` for the
reverse complement, and `vep-cli` loads what it writes. The output is `info.json`, one
transcript shard per chromosome and 1 Mb region under `transcripts/`, with the variation
data (`--include-variations`, or an `--auto` installation) one `all_vars.gz` with its `.csi`
per chromosome under `variations/`, the variation file Ensembl VEP's own cache carries,
with the regulatory files (`--regulatory-gff`, `--emars-gff` and `--motif-gff`, or an
`--auto` installation) the feature shards and `binding_matrices.json` under `regulatory/`,
then `CHECKSUMS.sha256`, one `sha256sum` line per file written, which `vep-cli` checks every
file it reads against. Until the build ends that file is one line that lists no file, so
`vep` refuses by name a directory a build left early.

## How it works

Follow one protein-coding transcript. The transcript model is the GFF3 named by `--gff3`,
else the GTF named by `--gtf`, else the release's Ensembl GFF3, which
`download::resolve_gff3` downloads from the Ensembl FTP (File Transfer Protocol) site into
the build's work directory under a `vep_cache_builder_<species>_<assembly>_<release>`
name. The work directory is `<output-dir>/tmp/vep_cache_builder/<species>_<release>_<assembly>/`
(`download::WorkDir`); nothing is downloaded into the system temporary directory. The peptides come from the protein FASTA named by `--fasta`,
else downloaded the same way by `download::resolve_protein_fasta`; with `--genome-fasta`,
no `--fasta` and no `--include-predictions`, no protein FASTA is read at all. `main` hands the
resolved paths to `main::build` as a `BuildInputs`, and `build` removes
`<output-dir>/tmp/vep_cache_builder/` once the cache is built, the work directory and any
that a failed build of another release, assembly or species kept beside it, before the
checksum manifest lists the output directory's files.
A build that fails with a download in the work directory keeps it and names it at the end
of its error line, and the next run reuses the files it holds. Under `--auto`, `install::install`
fills the same `BuildInputs` instead: it refuses a release directory that already exists,
reads the `CHECKSUMS` file of each directory of the site (or of the `--site` tree) as
that directory's listing, places each file under `<dir_cache>/tmp/<species>_<release>_<assembly>/`
and checks its BSD `sum` (the 16-bit rotating checksum and 1 KiB block count the listing
carries) against the listing, keeping a downloaded file that still verifies and refusing a
mismatch by file name. The genome FASTA goes into the cache directory under construction:
the `dna_index/` file with its published `.fai` and `.gzi`, else the `dna/` file, which
`fasta_index::install_fasta` keeps when it is already BGZF or writes again in BGZF blocks,
then indexes with the two writers of `vep-fasta`, `index::write_fai` over the inflated
stream (`.fai`) and `index::write_gzi`, a walk of the block headers (`.gzi`). The build then runs as below with that FASTA as the genome,
the manifest lists the FASTA with the rest, the cache directory is renamed to
`<dir_cache>/<species>/<release>_<assembly>` and the work directory is removed.
A `dna_index/` file that does not open through its published indexes is installed as a
`dna/` file is, with a warning that carries the refusal.
Before the first fetch `install` also looks up the three parts an installation brings
beside the transcripts: the variation VCFs (`VariationVcfs::find`), the table dumps the
matrices are read from (`prediction_tables`) and the regulatory-feature and motif-feature
GFF3 of the genome's `annotation/` directory (`RegulatoryFiles::find`, which passes over
the EMAR file, Ensembl VEP's own cache holding no EMAR). It looks up the same way what two
of the parts read beside those files. For the known variants (`VariantCellSources::find`):
the tables of the variation database that hold the PubMed identifiers (`variation`,
`variation_citation`, `publication`), the synonyms (`variation`, `variation_synonym`,
`source`) and the failed flag (`variation`, `failed_variation`), and the failed-variants GVF
of the genome's GVF directory; and the `variation_feature` and `seq_region` tables
(`SequenceFeatureTables::find`), which the build of a species with one genome-wide variation
VCF asks for when a sequence of the cache has no record in that file. For the regulatory
tree (`RegulatoryCells::find`): the activity file of the `annotation/` directory and, in the
directory of the funcgen database named beside the core database, the tables of the
matrices' counts (`binding_matrix`, `binding_matrix_frequencies`), of their factors
(`binding_matrix_transcription_factor_complex`, `transcription_factor_complex`) and of the
regulatory build's version (`regulatory_build`). One of these the site does not publish is
logged in one line, `No <what> is read for the known variants: <reason>` or
`No <what> is read for the regulatory tree: <reason>`, and the cache is built without what it
carries. It fills `BuildInputs` with each
part the site publishes a file of and leaves out, without an error, a part it publishes
none of; a `homo_sapiens` listing that names some per-chromosome variation files and lacks
another, or lacks one `--chromosomes` names, is refused. Under `--include-predictions` the
matrices' source is Ensembl's MySQL server instead of the dumps. After the line that names
the installed directory, `close` writes what `left_out` found: with no part left out, one
line saying the cache holds each part of Ensembl VEP's own cache that a file of the release
carries; else one warning for each part left out,
`No <part> in the cache of <species> <assembly>: <reason>`; and a warning of the same form
for each cell of a part that was brought which the cache holds no value of: the gnomAD
frequencies of `homo_sapiens`, which no file of the release carries, and each cell the built
cache's `info.json` names under `left_out` (the PubMed identifiers, the variant synonyms and
the failed variants of the variation data; the cell-type activity and `HIGH_INF_POS` or
`MOTIF_SCORE_CHANGE` of the regulatory tree), with the reason its lookup kept, the cells of
one reason in one warning, naming a conversion of Ensembl VEP's cache as what carries them.
Between the warnings of the parts and of those cells, one of the same form names the
sequences of the cache that hold no known variant, the first ten by name and the rest by
their count: `sequences_without_known_variants` reads them from the built cache, each
sequence of `contig_order` whose variation file is absent or holds no record.
A build named one input at a time ends on the same warning, naming the build
(`sequences_warning_of_a_build`), between its own lines and the line of its cells.
The installation's last warnings are the lines every
build ends on for what its cache leaves out (`left_out.rs`): the columns of the core table
dumps a build did not read, the missing source of a GTF or GFF3 read as the transcript
model, the transcripts whose translation names no protein, and the 1000 Genomes `AF`
recovered from the continental frequencies. `info.json` lists the same parts by name under
`left_out`, with the cells of the variation files and of the regulatory tree that hold no
value, and a `vep` run reads the list to name what its command asks for.
`fasta::parse_protein_fasta` reads the FASTA once and maps
each header's `transcript:` token, whole, to its
peptide and to the protein identifier the header opens with and its version.
`fasta::record_of` finds a transcript's record: the one whose token is the transcript's
identifier, a dot and its version when the model states a version and the file holds that
token, else the one whose token is the identifier alone. No dot of the token is cut at,
so an identifier with dots of its own (`AT1G01020.2`, `Y74C9A.3.1`) takes its own record,
and a token that names another version of a transcript is not that transcript's.
A record is its transcript's peptide only when its length is one the transcript's coding
sequence translates to: the whole codons of the translateable sequence, the cDNA coding
span led by the start exon's phase padding, or one fewer, the stop codon's residue left
out. A stop character that closes a record is no residue: it is neither counted nor
stored, as Ensembl VEP's cache peptide ends at its last residue.
A record of another length is set aside and the transcript is built as one with no
record, its peptide translated from the genome when one is given and absent otherwise;
the build's one warning counts the translated transcripts without a record and the
records set aside. A transcript with no coding sequence keeps the record that names it.
`gff3::parse_gff3` reads the gene, transcript, exon and CDS (coding sequence)
lines into a `GeneRecord` holding `TranscriptRecord`s, each with its sorted `ExonRecord`s
and `CdsRecord`s. On an Ensembl GFF3 read for the cache's content (a build with a protein
FASTA or the core dumps, `gff3::Reading::EnsemblCache`) a gene is any line whose `ID` is
`gene:<id>` and a transcript any line whose `ID` is `transcript:<id>`, whatever the type
column says (`gene`, `ncRNA_gene`, `pseudogene`; `mRNA`, `lnc_RNA`, `miRNA`,
`unconfirmed_transcript` and the rest), so every gene and transcript type reaches the cache;
from the genome alone (`gff3::Reading::GffSource`) the same file is read by the type list
below, as Ensembl VEP's `--gff` source reads it, and its TEC (`unconfirmed_transcript`) and
scaRNA (`scRNA`) lines are skipped. The reader names each sequence as the genome FASTA
holds it (`gff3::SequenceNames`): the UCSC `chr` before a digit or a capital letter is
dropped (`chr1`, `chrUn_KI270302v1`), as is a `chr` whose shortened name the FASTA holds
where it does not hold the written one; any other name the three letters begin
(`chromosome_1`) is kept as written. `chrM` or `M` is `MT`, except under a FASTA whose
reader reaches no record as `MT` and which names an organelle genome of its own (`Mito`),
where the name is kept as written. On GRCh37 `chrM` and `M` also name the 16,571-base
mitochondrial sequence of the UCSC hg19 naming, so there the record read as `MT` counts
only where it has the assembly's 16,569 bases: over a record of another length the name
is kept as written, a build says so once a region, and a run that reads the file beside a
cache leaves the region's transcripts out with a warning. Neither the name nor the length
says which of the two sequences a file's models count on, so on GRCh37 a mitochondrial
region's coding transcripts are also tested on the FASTA record they are cut from
(`builder::grch37_mitochondrial_regions_off_their_record`, which names the transcripts
left untested): where at least one is tested and none of those opens on a start codon of
the region's codon table and runs to its last codon without a stop codon, a build ends
with an error that names the region, the record and its length and writes no cache, and a
run leaves the region's transcripts out with a warning. The reader
drops the `gene:` and `transcript:` identifier prefixes, and keeps what the
Ensembl core API writes on each line: the gene and transcript `version`, `ccdsid`,
`transcript_support_level`, the `tag` codes (`gencode_basic`, or `basic` on the GRCh37
release-87 file, `gencode_primary`, `Ensembl_canonical`, `MANE_Select`,
`MANE_Plus_Clinical`), the exon `ensembl_phase` and `ensembl_end_phase`, and the CDS
line's `protein_id` and translation `version`, which name the ENSP (Ensembl protein)
identifier. A transcript of biotype `artifact` or tagged readthrough (`readthrough_transcript`
on the release 116 GRCh38 file, `readthrough_tra` by the core API's attribute code) is a
transcript of the tree; `builder::build_cache` leaves both out of a build that reads a
protein FASTA or the core dumps, the cache's content, as Ensembl VEP's own dumper does
(`Database/Transcript.pm`), and keeps them on a build from the model and the genome
alone, as Ensembl VEP's `--gff` and `--gtf` sources do. A model holding no other
transcript leaves a build of the cache's content nothing to write: `build` refuses it
before anything is written, as it refuses a model from which no transcript is read, with
`No transcript was kept from <file>` and the transcripts read and set aside under each
reason.

Any other line, and every line of a GTF read by `gtf::parse_gtf_model`, is read as
Ensembl VEP reads a GFF3 or GTF given to it as an annotation source (`BaseGXF.pm`,
`GFF.pm`, `GTF.pm`): only the feature types of Ensembl VEP's list are taken (a GTF's
`gene`, `transcript`, `exon`, `CDS` and `stop_codon`; a GFF3's `gene`, `mRNA`, `ncRNA`,
`lnc_RNA`, `miRNA`, `primary_transcript`, `pseudogene`, the gene segments and the rest),
the skipped types logged once; a record is named by `ID`, `Name`, `id` or `name` (a
GTF's `<type>_id`) and placed under each record its `Parent` list names (a GTF's
`transcript_id` or `gene_id`); a record with no parent is a gene when its type names one
and otherwise a transcript with no gene line; a record whose parent list opens with a
parent in no line, or in a line of a type not read, is placed as its other parents place
it, else dropped with its children, the parents named in a warning, the two cases apart; a
stop codon is one more CDS segment. An attribute whose value is empty or `0` is
absent, as it is to Ensembl VEP, so each of these rules takes its next alternative, and a
GTF source column of `0` names no biotype. A record whose parent list opens with `0` and
names nothing set (`Parent=0,0`; a GTF transcript's `gene_id "0"`) is dropped without a
warning. One whose list opens with an empty field is a record of its own, or the
transcript of a gene a later field names. A transcript's biotype is
its `biotype`, `transcript_type` or `transcript_biotype`, else, for a GFF3, the RefSeq
rules (`mRNA` is `protein_coding`, `ncRNA` takes `ncrna_class`, a `V_gene_segment` is
`IG_V_gene`, a transcript of a gene described `microRNA ...` is `miRNA`, else `gbkey`)
and, for a GTF, the source column; a transcript without one is dropped with a warning, as
is a `protein_coding` transcript without an exon line. Its identifier is its
`transcript_id`, else the record's, either without a `gene:` or `transcript:` prefix. Its
gene identifier is its `gene_id`, else the `GeneID` of its `Dbxref`, else the gene line's;
its symbol its `gene_name` or `gene`, else the gene line's `name` or `Name`; its version
`version`, else `transcript_version`; its translation's version the CDS line's `version`,
else 1 (a build of the cache's content takes, when the line carries none, a GTF CDS line's
`protein_version`, else the core dumps' version when they hold the translation, else the
protein FASTA header's, else none); its CCDS identifier `ccdsid`, else a GTF's `ccds_id`,
as the line carries it. A `version` that is not a whole number, or one above 4294967295,
the largest the cache's field holds, on any gene, transcript or CDS line refuses the build
with one `ERROR: [cache-builder]` line naming the file, the line and the text or the limit
(`01` is the whole number 1), where Ensembl VEP prints the text as written; a protein FASTA
header's version above the limit is refused the same way, and an identifier ending in `.`
followed by anything but digits (`.cds1`, `.p1`) states no version, as does a header
identifier that is the CDS line's `protein_id` whole, whose last number is the
identifier's own (`AT1G01020.2`). The
tags read are the ones Ensembl VEP keeps: `Ensembl_canonical`, `MANE_Select`,
`MANE_Plus_Clinical`, `gencode_basic` (a GTF's `basic`, the GRCh37 release-87 and GENCODE
spelling, is read as `gencode_basic` through `gtf::KEPT_TAGS`, the table the tags beside a
GFF3 use) and `gencode_primary`, with `transcript_support_level` or `tsl`; the reader drops
no readthrough or `artifact` transcript. A `cds_start_NF` tag of a transcript's own line is
read beside them for the GRCh37 mitochondrial test alone
(`gff3::TranscriptRecord::cds_start_nf`) and gives the record no attribute.

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
minute and 3 GB of memory. A release the site publishes no variation database for leaves
the phenotype tables unread, with a warning, and every gene counts 0, the value Ensembl
VEP's dumper writes for such a release, so `gene_phenotype` is present on every transcript
and the `GENE_PHENO` column stays empty. `predictions::fetch_predictions`
hashes every peptide with MD5 (Message-Digest algorithm 5) and fetches the matching SIFT
and PolyPhen matrices from Ensembl's public MySQL server in batches of 500 hashes,
`gtf::parse_gtf_tags` reads the `#!` header keys and the per-transcript `tag` values of the
GTF (Gene Transfer Format), and `vep_fasta::IndexedFasta` opens `--genome-fasta`.
An `--auto` installation without `--include-predictions` takes the matrices of the same
hashes from `predictions::read_prediction_dumps`, which reads the release's dumps of the
variation database's `attrib_type`, `attrib`, `translation_md5` and
`protein_function_predictions` tables.
A `transcript` line's tags are held under its `transcript_id` whole, with its version
(`version`, else `transcript_version`) behind a dot when the line carries one, and a
transcript takes the tags of the line named by its identifier, a dot and its version when
the model states a version and the file holds that name, else of the line named by its
identifier alone, as it takes its protein FASTA record: transcripts that share a first
component (`T1.1`, `T1.2`) take their own lines' tags, and the line of another version of
a transcript tags nothing.
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
`codon_table` is the sequence region's `codon_table` attribute when the core dumps carry
one, else the table `--codon-table <region>=<table>` names for the region, else 1, the
table Ensembl VEP's dumper writes for a region without the attribute. The two are settled
before any shard is built: a region the dumps and the flag give different tables stops
the build, and so does a region named like an organelle genome (`MT`, `M`, `chrM`,
`chrMT`, `Mito`, `MtDNA`, `mitochondrion_genome`, `Mt`, `Pt`, `ChrC`, `ChrM`,
`chloroplast`, `plastid`, or any name that holds `mito`, `chloro` or `plastid`, such as
`OVOC_MITOCHONDRIAL`, in any letter case) whose table neither names, with
`ERROR: [cache-builder] codon table of sequence region <name> unknown: pass --core-dumps or --codon-table <name>=<n>`,
since the default of 1 would translate its transcripts under the nuclear code. A
translation's start, the 5' end of the first CDS record in transcript order, and its end,
the 3' end of the last, each lie inside an exon of the transcript: one that lies in none
refuses the build before any shard is written, naming the file, the line and span of the
CDS record that sets it, the transcript and the base. A CDS record that leaves its exon
between the two ends is read as written. A coding sequence that still opens outside its
start exon, as when the exons' ranks run against the strand, is refused the same way with
its whole span. Neither end is written at an offset outside the exon it is counted in:
where the last exon in transcript order that the coding sequence overlaps does not hold
its last base, as when the ranks run against position there or an exon nests in the last
coding one and comes after it, the build is refused, naming the offset, that exon and its
length. The coding sequence's last base is placed at or after its first: on the reverse
strand an exon nested in the one the translation starts in comes before it by position,
so a translation that ends in the nested exon is refused, naming the CDS record that sets
the end and both exons. Exon ranks that run against position anywhere else, as 1, 3, 2 over three exons with
the translation ending inside the exon ranked last, are also refused, naming the first exon in
rank order that is out of place, the exon ranked before it and both ranks: Ensembl VEP's
source places each exon at its `rank` and reads an exon's bases as an intron's there. A CDS
line's phase column is text that source takes for an exon's phase (`0`, `-1`, a number that
reads 1 or 2, or text holding a `.`); any other text refuses the build, naming it, where
Ensembl VEP ends its run. With the core dumps a one-exon transcript whose `translation`
row starts or ends the coding sequence at other positions than its CDS records is built
as the row places it (`translation_row::coding_records`), as Ensembl VEP's cache holds it,
and one line counts such transcripts; a transcript of several exons whose row differs keeps
its CDS records and is counted in a warning, since the dumps name the row's exons by
database identifier alone. The
peptide comes from the FASTA map or, for a transcript the
map lacks, from `translateable_seq` translated under that codon table as Ensembl VEP
translates a transcript read from a GFF3 or GTF (whole codons only, the terminal stop
dropped, a `CTG` or `TTG` initiation codon read as `M`, an ambiguous codon resolved when
every expansion agrees), the protein identifier from the CDS line and the translation
version from the CDS line (`version`; on a build of the cache's content a GTF's
`protein_version` too), else the core dumps when they hold the translation, else the
FASTA header (`builder::translation_version`), `ccds` from `ccdsid` or `ccds_id` (or the dumper's `-`), and
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
written (the chromosomes of the core dumps' karyotype first, by their `karyotype_rank`
attribute, then autosomes by number, X, Y, MT and the rest by name) and whose `variation_cols`
and `var_type tabix` name the variation files' columns and layout only when the files were
written. With `--core-dumps`,
`info::write_chr_synonyms` writes `chr_synonyms.txt` beside it from the `seq_region` and
`seq_region_synonym` tables, the pairs Ensembl VEP's dump pipeline writes (each
sequence region's name against its synonyms, linked one step further through the
synonyms already recorded), so `NC_000021.9` or `CM000683.2` in an input resolves to the
cache's `21`. Last, `vep_core::cache_manifest::CacheManifest::write` lists every file of the
output directory in `CHECKSUMS.sha256`, the file `vep-cli` checks each cache file it reads
against. `builder::build_cache` writes that file first, before its first shard, as the one
line of `CacheManifest::mark_unfinished`. A build that fails or is killed between the two
writes leaves a directory `CacheManifest::load` refuses as `ManifestError::Unfinished`,
which `vep` prints under the `cache` component; building into the same directory completes
it. The build writes into the output directory as it stands: a file the directory holds
under another name is kept and listed.

Without `--core-dumps` the cache states what the GFF3 states and no more: no symbol
source or HGNC identifier, no gene phenotype count, no MANE accession, no APPRIS, no
mature miRNA arms, no RefSeq or UniProt lists, no protein features and no sequence edits, so the `SYMBOL_SOURCE`,
`HGNC_ID`, `GENE_PHENO`, `MANE_SELECT`, `APPRIS`, `miRNA`, `RefSeq`, `SWISSPROT`, `TREMBL`,
`UNIPARC`, `UNIPROT_ISOFORM` and `DOMAINS` columns stay empty on such a cache, and an
organelle region's codon table is named with `--codon-table` instead. The
Ensembl cache lists a translation's protein features in the database's row order, which
no dump carries; a native cache lists each row of the `protein_feature` table once, by
feature identifier, so `DOMAINS` may name the features in another order. An Ensembl cache
can also hold other features than the release's table: one row for each `interpro` row a
hit joins, and the rows of another release's table (`docs/cache-setup.md` gives both for
`saccharomyces_cerevisiae`).
An Ensembl cache can also name an analysis by another label than the release's
`analysis_description` table, whose label a native cache writes. The Ensembl Genomes
site's release 63 cache of `arabidopsis_thaliana` names three analyses so, on 303 of the
2,396 rows that `--everything` writes for the input of
`tests/golden/116/TAIR10-arabidopsis` (`docs/cache-setup.md` tabulates them).

With `--include-variations`, `download::resolve_variation_inputs` scans
`--variation-vcf-dir` for the per-chromosome `*-chr*.vcf.gz` files and, under their Ensembl
names, the genome-wide `1000GENOMES-phase_3.vcf.gz`, `<species>_somatic.vcf.gz`,
`<species>_phenotype_associated.vcf.gz` and `<species>_clinically_associated.vcf.gz`, or
downloads them into the work directory (Variant Call Format, each with its `.csi` index; a
genome-wide file the site does not publish, the 1000 Genomes file on GRCh37, is left out
with a warning).
That layout is the one Ensembl publishes for `homo_sapiens` alone. Every other species
with variation data has one genome-wide germline file, `<species>.vcf.gz`: a directory
holding no per-chromosome file is read from it, and it is what the download and `--auto`
fetch for a species other than human, for every sequence the file holds records on.
`--chromosomes` narrows that to the chromosomes it names, in a build from named inputs and
with `--auto` alike.
A directory with neither form, or a site without the file a download fetches, ends the
build naming the files looked for. `--auto` reads the same files without the flag, with
the ClinVar tables of the release's variation database; where the site publishes none of
the files it leaves the variation data out and its closing lines name the listing, and it
refuses before any download a `homo_sapiens` listing that names some per-chromosome files
and lacks another, or lacks one `--chromosomes` names.
An installation also fills three cells no VCF carries (`variant_cells.rs`): `pubmed` from the
`variation`, `variation_citation` and `publication` table dumps of the variation database,
`var_synonyms` from `variation`, `variation_synonym` and `source`, without the sources whose
name holds `dbSNP`, as in Ensembl VEP's dump, and `failed` from `variation` and
`failed_variation`. The release's failed-variants GVF (`<species>_failed.gvf.gz`, read by
`gvf.rs`) holds every feature of every failed variant: one whose name a VCF record of its
chromosome carries flags that record, and any other becomes a record of its own, since the
VCFs leave most failed variants out. `variation` and `variation_synonym` hold a row for
every variant and for every name of one, so each is read once as a stream (`dump_rows.rs`)
and joined through sorted runs on disk (`sorted_runs.rs`), never held in memory. The
published `variation_synonym` table does not record the order in which the dump's query
returned the names of one source, so they stand in the reverse of the table's row order,
which can differ from the order of Ensembl VEP's cache. A list longer than the 1,024
characters Ensembl VEP's cache ends one at is stored at its whole items within that length,
the text `vep` prints from the cut value. For a species whose variation VCF is one
genome-wide file, a sequence of the cache that file holds no record on takes its known
variants from the displayed rows of the `variation_feature` table
(`variation/feature_table.rs`), fetched only when such a sequence stands; a sequence the
file holds records on is not read from the table. `--chromosomes` narrows none of the table
dumps. A build from named inputs reads none of these files: its variation files carry no
`pubmed`, `var_synonyms` or `failed`.
`variation::chromosome_records` reads one chromosome at a time: the per-chromosome file
whole, or the chromosome's records of the species' genome-wide germline file, and each
further genome-wide file through its index (or read through when it has none),
keeping each record whose `FILTER` is `PASS`, `.` or empty, adding the genome-wide
records the germline file lacks (a record is the same one when its name, position
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
  regulatory/
    binding_matrices.json
    1/
      1-1000000.json
```

The `variations/` tree exists only with `--include-variations` or in an `--auto`
installation of a genome the release publishes variation VCFs for, and the `regulatory/`
tree only with the regulatory flags or in an installation of a genome it publishes
regulatory files for. A transcript shard
carries no `chr` field: the loader takes the chromosome from the directory name.

## Modules

| Module               | Role                                                                                                                                                                                         |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src/lib.rs`         | The library `vep-cli` links: it documents `gene_model` alone, and the other modules it declares are the builder's own, hidden from the API documentation or private. |
| `src/gene_model.rs`  | `read`: a GFF3 or GTF that a run names beside a cache (`vep --gff`, `--gtf`, or a `--custom` source of format `gff` or `gtf` in transcript mode), read into the transcript records a cache build from the file and a genome alone writes, built in memory and never written; `Format`, `GeneModel` with the `OffRecord` values of its `off_record` field, and `GRCH37_MITOCHONDRIAL_NAMING`, the words of the warning for a region kept apart from GRCh37's `MT`. |
| `src/main.rs`        | The `clap` argument struct, the `BuildInputs` either path fills (the flags one input at a time, or `install::install` under `--auto`) and `build`, the run order: choose the transcript model (the GFF3, else the GTF, else the download) and whether a protein FASTA is read, resolve and read the core dumps, parse the FASTA and the model, refuse a model that leaves the build no transcript, settle each sequence region's codon table, fetch predictions, load the genome FASTA and the GTF tags, build the transcript cache, then the regulatory tree, then the variation files, then `chr_synonyms.txt` and `info.json` |
| `src/download.rs`    | The published tree as data (`Layout`: the Ensembl release, GRCh37 under `pub/grch37/` with its gene set capped at release 87, an Ensembl Genomes division, a collection, the rapid release; in each the directory and name of every input, `dna_index/` then `dna/`), the species-list reader, the database lookups, a typed "not published" answer, the answers of a server read as absent (404, and 403 for an index page or for a file a tree may lack), the work directory of a build (`tmp/<species>_<release>_<assembly>` under the cache root for `--auto` and `tmp/vep_cache_builder/<species>_<release>_<assembly>` under the output directory otherwise, its removal, there with every other build's, and the error line of a failed build, which names the directory when it holds a file and removes it when it holds none), the GFF3 and protein FASTA downloads into it with reuse, the `--core-dumps` directory fill, the local variation directory scan, the parallel per-chromosome VCF download of `homo_sapiens`, the one genome-wide germline VCF of any other species and the further genome-wide variation VCFs with their indexes, all into the work directory; the names an installation looks up for the cells of its parts: the failed-variants GVF and its directory (`failed_gvf_name`, `variation_gvf_dir`) and the funcgen database beside the core database (`funcgen_database`) |
| `src/install.rs`     | `--auto`: the `Site` (a URL or a local directory), `locate` (finds the genome's site, tree and core database: by the release's species list, else in a release tree that lists the genome's GFF3, a site of several that cannot be read passed over with a warning), the two refusals (a genome of a collection, an assembly of the rapid release), each directory's `CHECKSUMS` listing, the GFF3, the GTF, the protein FASTA, the core tables, the genome FASTA, the variation VCFs, the table dumps of the matrices, the regulatory GFF3 files and the files of the cells (`VariantCellSources`, `SequenceFeatureTables`, `RegulatoryCells`: the variation tables, the failed-variants GVF, the activity file and the funcgen tables) each looked up in its listing before the first fetch, the BSD `sum` check of every fetched file, the parts a site publishes no file of and the cells the built cache names as left out (`left_out`), the sequences of the built cache that hold no known variant (`sequences_without_known_variants`) and the closing lines that name them (`close`), the work directory, the genome FASTA choice (primary assembly, else toplevel, never a masked file) and the move of the built cache into `<dir_cache>/<species>/<release>_<assembly>` |
| `src/fasta_index.rs` | The cache directory's genome FASTA: BGZF detection by the `BC` extra subfield, gzip-to-BGZF recompression, and the call of the two index writers, `vep_fasta::index::write_fai` and `write_gzi`, which live in `vep-fasta` |
| `src/fasta.rs`       | Protein FASTA reader: each header's whole `transcript:` token to peptide, header protein identifier and its version, plain, gzipped or bgzipped; the lookup of a transcript's record by its identifier and version |
| `src/gff3.rs`        | GFF3 reader into `GeneRecord`, `TranscriptRecord`, `ExonRecord` and `CdsRecord`, and the gene tree a GFF3 or GTF builds by parent identifiers (Ensembl VEP's feature-type list, biotype rules, identifier and symbol fallbacks, and the truth test under which an empty or `0` value is absent); chromosome normalisation, attribute URL-decoding, the `tag` codes, `ccdsid` or `ccds_id`, `transcript_support_level`, the exon phases, the CDS `protein_id` and `version`, and the readthrough tag |
| `src/core_dumps.rs`  | The core and variation table dumps: MySQL escape decoding, the joins of `transcript_attrib`, `translation`, `protein_feature`, `xref`, `object_xref`, `gene`, `seq_region_attrib`, `meta`, `phenotype_feature` and `phenotype` into per-transcript and per-gene facts; the table reader `info.rs` uses for `seq_region_synonym` |
| `src/info.rs`        | `info.json` (species, assembly, version, the `sift` and `polyphen` markers, the `regulatory` marker and `regulatory_schema`, `cell_types` for a tree baked with an activity file, source versions, `variation_cols` and `var_type`, `contig_order`, `left_out`) and `chr_synonyms.txt` from the `seq_region_synonym` table |
| `src/left_out.rs`    | What a build leaves out of its cache: the names `info.json` lists under `left_out` and the warnings the build ends on, one for each part |
| `src/regulatory.rs`  | The regulatory tree from the Ensembl regulatory GFF files: each regulatory, EMAR and motif record read as Ensembl VEP's GFF regulatory source reads it (the accepted types, the identifier fallbacks, the promoter widening of `--extended-promoters` with a bound that is not a whole number ignored and warned about, a record the file holds twice written once and two distinct records of one file sharing an identifier both written with the collision warned about once), the shards per 1 Mb region of each feature's start in the shape `vep-cli` reads, each feature carrying the file it came from (`source`: the regulatory file 0, the EMAR file 1, the motif file 2) and its place among that file's accepted records (`ordinal`), and `binding_matrices.json` of one entry per motif matrix: its stub (span and factors, no counts) from the GFF files alone, or its counts, name, source, threshold and factor names where the build is given the funcgen table dumps. With an activity file each regulatory feature carries one activity code per cell type and the tree returns the `cell_types` list; with the genome FASTA each motif carries its sequence on its own strand; with the `regulatory_build` table the tree returns the build's version (`regbuild`) |
| `src/line.rs`        | One line of a GFF3, a GTF or a regulatory GFF for the three readers: read under the length limit, its byte-order marks taken off, told from a blank line and a comment, and split into columns twice, as the index build of Ensembl VEP's release ends them (at a tab or a NUL byte), which decides whether that build places the line, and as an annotation reader ends them (at a tab), which decides what a record holds |
| `src/count.rs`       | A count beside its noun as the log lines word them: the singular beside one, the plural beside any other count, and the verb or pronoun that agrees with it                                  |
| `src/gtf.rs`         | A GTF as the transcript model (`gene`, `transcript`, `exon`, `CDS` and `stop_codon` lines into the gene tree, the `#!` header keys) and, beside an Ensembl GFF3, its `transcript` lines' five tags kept as transcript attributes; one tag table serves both readings (`basic` and `gencode_basic` both become `gencode_basic`) |
| `src/mapper.rs`      | cDNA-to-genome mapper pairs, `cdna_coding_start` and `cdna_coding_end` for one transcript, with exon 1's phase as the start phase                                                              |
| `src/translation_row.rs` | The coding sequence a build that reads the core dumps gives a one-exon transcript whose `translation` row places it otherwise than the model's CDS records, and what becomes of a row the build cannot place |
| `src/builder.rs`     | The serialised cache types, the dumper's `artifact` and readthrough drop on a build that reads a protein FASTA or the core dumps, exon ordering and phases, introns, the translation's bounds and the refusal of a coding sequence that does not sit in its exons, the spliced cDNA cut into UTRs and `translateable_seq`, the peptide translated from it when the protein FASTA lacks the transcript, the codon table, seq edits, attributes and xref lists from the GFF3 or the dumps, the per-region HGNC identifier copy, 1 Mb binning and the transcript shards |
| `src/predictions.rs` | SIFT and PolyPhen matrix fetch from Ensembl's MySQL server, or read from the table dumps a release publishes of the variation database, keyed by the MD5 of the peptide and stored base64-encoded under the Ensembl matrix object's keys (`analysis`, `sub_analysis`, `translation_md5`) |
| `src/variation.rs`   | The variation VCF readers (a chromosome's own file whole, a genome-wide file through the index), each record in the variation database's form, the 1000 Genomes frequency join, the ClinVar columns from the variation dumps, the `pubmed`, `var_synonyms` and `failed` cells an installation gives a record, the failed variants it adds from the failed-variants GVF, and the `all_vars.gz` writer with its `.csi` |
| `src/variation/feature_table.rs` | The known variants of a sequence no variation VCF holds a record on, read from the displayed rows of the variation database's `variation_feature` table and sorted on disk by sequence and start |
| `src/variant_cells.rs` | The `pubmed`, `var_synonyms` and `failed` cells of the known variants, joined from the variation database's table dumps by the variation's name: each table read once as a stream, the rows that take part sorted in runs on disk and merged, a list of more than 1,024 bytes stored as `vep` prints its first 1,024 |
| `src/gvf.rs`         | The failed-variants GVF of a release, each feature read back as the record the variation database holds of it: its name, span and allele string |
| `src/dump_rows.rs`   | A table dump read row by row as bytes, keeping the values of the columns asked for, with the escapes `LOAD DATA` reads |
| `src/sorted_runs.rs` | Records sorted on disk: a writer that spills its buffer as a sorted run each time it fills, and a merge that reads the runs back as one stream in key order |

## Command line

The binary has no subcommands. `--release` is required, and so are `--assembly` and
`--output-dir` unless `--auto` names the species and assembly; every other flag is optional. [docs/cli-reference.md](../../docs/cli-reference.md) documents
the `vep` binary, not this one, so the flags are listed here.

| Flag                          | Role                                                                                                                                                      |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--release <N>`               | Ensembl release number; names the download paths and becomes `cache_version` in `info.json`. For a species of Ensembl Genomes it is the Ensembl Genomes release, the number that site's paths and file names carry (63 beside Ensembl release 116) |
| `--auto <SPECIES> <ASSEMBLY>` | Install the release in one command: fetch every input from the site, check each against its directory's `CHECKSUMS`, build the cache as `<dir_cache>/<species>/<release>_<assembly>` and place the indexed genome FASTA inside it; refused beside the local-file flags, `--species`, `--assembly` and `--output-dir`. The cache carries the known variants with the 1000 Genomes frequencies, the ClinVar columns, the PubMed identifiers, the synonyms and the failed variants, the SIFT and PolyPhen matrices, read from the variation database's table dumps, and the regulatory and motif features with their cell-type activity and matrix counts, wherever the release publishes a file of them for the genome; a part it publishes no file of is left out, and the installation ends on one line saying the cache holds each part or on one warning for each part left out, with the reason, then on a warning for the gnomAD frequencies, which no file of the release carries, and for each cell whose file the genome's tree lacks. The release's species list says which tree the species is published in and names its core database: a vertebrate installs from the Ensembl release tree, a species of the plants, metazoa, fungi or protists division from that division's tree on the Ensembl Genomes site. A genome of a collection (every bacterium and archaeon, and the fungi and protists published in collections), whose core database is shared, and an assembly of the rapid release (the human pangenome assemblies), which has no table dumps, are refused before any download with one line naming what the tree lacks, cache conversion and a build from local copies |
| `--dir_cache <DIR>`           | The cache root `--auto` installs under, default `$HOME/.vep`; the downloads land under its `tmp/` and are removed once the cache is in place (also `--dir-cache`) |
| `--site <URL or DIR>`         | Where `--auto` fetches from: a URL or a local directory holding the `pub/release-<N>/` tree of the Ensembl site or of the Ensembl Genomes site. Without it a species is looked up on `https://ftp.ensembl.org`, then on `https://ftp.ensemblgenomes.ebi.ac.uk` when the first does not hold it; one of the two that cannot be read is passed over with a warning. A server that answers 403, not 404, for a file it does not hold installs as the other does |
| `--assembly <NAME>`           | The assembly name as the release's file names carry it (`GRCh38`, `GRCh37`, `TAIR10`); names the downloads and becomes `assembly` in `info.json`. `GRCh37` selects the `grch37/` FTP tree, and `GRCh37` or `GRCh38` the prediction database |
| `--species <NAME>`            | Ensembl species directory name, default `homo_sapiens`                                                                                                    |
| `--output-dir <DIR>`          | The cache directory to write. A file the build downloads lands under its `tmp/vep_cache_builder/<species>_<release>_<assembly>/`; `tmp/vep_cache_builder/` is removed once the cache is built, and a build that fails with a file in its work directory keeps it, named in the error line |
| `--gff3 <FILE>`               | Local GFF3, plain or gzipped, as the transcript model instead of the download: Ensembl's own, or any GFF3, read as Ensembl VEP reads a `--gff` source     |
| `--fasta <FILE>`              | Local protein FASTA instead of the download; with `--genome-fasta`, no `--fasta` and no `--include-predictions`, no protein FASTA is read and each peptide is translated from the genome |
| `--gtf <FILE>`                | Local GTF. Beside `--gff3`, the Ensembl GTF of the same release: the `cds_start_NF` and `cds_end_NF` attributes and the `source_versions` keys that the GFF3 lacks, when no core dumps are read; it is read for those alone, so no coordinate of it is read and none refuses the build. Without `--gff3`, the transcript model itself, read as Ensembl VEP reads a `--gtf` source |
| `--core-dumps <DIR>`          | Directory of the release's core and variation MySQL table dumps (`<table>.txt.gz`, or `<table>.txt`); a table missing from it is downloaded there. Supplies every transcript attribute, xref, protein feature, symbol source, HGNC identifier, codon table, source version, chromosome synonym and gene phenotype count the GFF3 lacks |
| `--codon-table <REGION=TABLE>` | The codon table of a sequence region the core dumps state none for, one of the genetic codes `Bio::Tools::CodonTable` ships or 26; repeatable or comma-joined (`MT=2` for a vertebrate's mitochondrion built without `--core-dumps`, `Mito=3` yeast's, `mitochondrion_genome=5` the fly's, `Pt=11,Mt=1` the Arabidopsis plastid and mitochondrion). A region named like an organelle genome whose table neither the dumps nor this flag states stops the build |
| `--genome-fasta <FILE>`       | Indexed genome FASTA (`.fa` beside its `.fai`): `translateable_seq` for every coding transcript, the sequence behind codon-level consequences, and, with no `--fasta` and no `--include-predictions`, the peptides, each translated from it with no protein FASTA read |
| `--include-predictions`       | Fetch the SIFT and PolyPhen matrices from Ensembl's MySQL server; `homo_sapiens` on GRCh37 or GRCh38 only, any other species or assembly being refused by name before anything is downloaded, beside `--auto` too (a converted Ensembl cache carries the matrices Ensembl publishes for a species). An `--auto` installation without the flag reads the matrices of any species from the release's variation table dumps, and with it from the server instead: the route to the GRCh37 matrices (`--mysql-port 3337`), since the GRCh37 tree publishes the matrices table without the `translation_md5` table that keys it. The matrices are keyed by the peptide, so the protein FASTA is read (downloaded when `--fasta` is absent) whatever `--genome-fasta` says |
| `--mysql-host <HOST>`         | Prediction server, default `ensembldb.ensembl.org`                                                                                                        |
| `--mysql-port <PORT>`         | Prediction server port, default `3306`                                                                                                                    |
| `--include-variations`        | Also write `variations/<chr>/all_vars.gz` and its `.csi` from Ensembl's variation VCFs: for `homo_sapiens` the per-chromosome files with the genome-wide somatic, phenotype-associated, clinically-associated and 1000 Genomes frequency files; for any other species its one genome-wide file, `<species>.vcf.gz`. An `--auto` installation reads them without the flag, which changes nothing beside it |
| `--variation-vcf-dir <DIR>`   | Local directory of the `*-chr*.vcf.gz` files, or of the one genome-wide `<species>.vcf.gz` where it holds no per-chromosome file, and, under their Ensembl names, the further genome-wide files with their `.csi` indexes, instead of the download |
| `--variation-dumps <DIR>`     | Directory of the release's variation MySQL table dumps (`source`, `seq_region`, `attrib_type`, `phenotype`, `phenotype_feature`, `phenotype_feature_attrib`; nothing is downloaded): the ClinVar allele and somatic classification columns of the variation records, read with `--include-variations` for the chromosomes being built |
| `--chromosomes <A,B,...>`     | Chromosomes whose variation records to read; the default is 1 to 22, X, Y and MT under `homo_sapiens` and, for another species, every sequence its genome-wide variation VCF holds records on. The flag names the per-chromosome files of `homo_sapiens` and narrows the read of another species' one genome-wide file, with `--auto` as in a build from named inputs. It narrows no table dump an installation fetches; an installation of another species reads the `variation_feature` table for the named sequences its variation VCF holds no record on, and without the flag for every such sequence of the cache. The transcript build always covers every chromosome the model carries |
| `--regulatory-gff <FILE>`     | The Ensembl regulatory feature GFF3 (`regulatory_features`), written as the cache's regulatory shards: promoters, enhancers, CTCF binding sites and open chromatin regions, each with its type as `BIOTYPE`, read by a run with `--regulatory` |
| `--emars-gff <FILE>`          | The Ensembl EMAR GFF, regulatory features of type `EMAR` written beside those of `--regulatory-gff`                                                      |
| `--motif-gff <FILE>`          | The Ensembl motif feature GFF3: each motif names its binding matrix by identifier and the cache carries a stub of the matrix (its span and transcription factors, no counts), so motif rows carry `MOTIF_NAME`, `MOTIF_POS` and `TRANSCRIPTION_FACTORS` and no `HIGH_INF_POS` or `MOTIF_SCORE_CHANGE`; with `--genome-fasta` each motif carries its sequence on its own strand. An `--auto` installation reads the matrices' counts from the funcgen table dumps where the release publishes them, and its motif rows then carry both fields |
| `--extended-promoters`        | Write each promoter of `--regulatory-gff` at its `extended_start` and `extended_end`                                                                     |

Against Ensembl VEP's `INSTALL.pl`: `--auto <species> <assembly>` is its `--AUTO cf --SPECIES
--ASSEMBLY`, the cache and FASTA steps alone, since a `vep` binary needs no API, htslib or
BioPerl install and plugins are built in or loaded from `--dir_plugins`; `--dir_cache` is
`--CACHEDIR`, `--release` is `--CACHE_VERSION` and `--site` stands for both `--CACHEURL` and
`--FASTAURL`. The cache is built from the release files where `INSTALL.pl` downloads Ensembl
VEP's prebuilt cache archive. Of what that archive carries beside the transcripts, an `--auto`
install brings each part that a file of the release carries for the genome: the known variants
with their 1000 Genomes frequencies, PubMed identifiers, synonyms and failed flag, the SIFT
and PolyPhen matrices and the regulatory and motif features with their cell-type activity and
matrix counts. It ends on one warning for each part left out, with the reason, and on
warnings for the cells it did not bring: the gnomAD frequencies, which no file of the release
carries, and each cell whose file or table the genome's tree lacks, as the GRCh37 tree lacks
the tables of the PubMed identifiers and the variant synonyms; a conversion of that archive
carries each.
Between the two it names, in one warning, the sequences of the cache it read no variation
record on, where a conversion of that archive carries any known variant the archive holds.
`--include-variations` changes nothing beside `--auto`, and `--include-predictions` fetches the
matrices from Ensembl's MySQL server instead of the release's table dumps, the route to the
GRCh37 matrices, whose key table the GRCh37 tree does not publish; that tree publishes no 1000
Genomes file either. `--regulatory-gff`, `--emars-gff` and `--motif-gff` are refused beside
`--auto`, which reads no EMAR file, and a build with `--output-dir` takes them. An install
holds the Ensembl transcript set alone, since vep-rs builds no RefSeq or merged cache.
`INSTALL.pl` reaches a site other than Ensembl's through
`--CACHEURL` alone; `--auto` needs no `--site` for a species of an Ensembl Genomes division,
which the release's species list places on the Ensembl Genomes site, and `--release` is then
that site's release, the version in the name of the cache `INSTALL.pl` downloads from there.
A species' core database is the one the species list names, or the one the tree's `mysql/`
directory holds for it where the tree publishes no list; a name ending in the assembly's
number is tried only on a server that lists neither. Where `INSTALL.pl` downloads the cache of a bacterial genome or of a human
pangenome assembly, `--auto` refuses before any download, since neither tree has table dumps
of the genome's own, and names [cache conversion](../../docs/cache-conversion.md).
`--DESTDIR`, `--PLUGINS`, `--PLUGINSDIR`, `--PLUGINURL`,
`--NO_UPDATE`, `--NO_HTSLIB`, `--NO_BIOPERL`, `--NO_TEST`, `--CONVERT`, `--PREFER_BIN`,
`--QUIET`, `--TEST`, `--USE_HTTPS_PROTO`, `--GITHUB_TOKEN` and the `VEP_*` environment
variables have no counterpart. `--HELP` is `--help`; `--NO_PLUGINS` and `--SKIP_PLUGINS` have
no counterpart, because the supported plugins are built into `vep`, and `--BIOPERL` has none,
because nothing depends on BioPerl. The FASTA step takes `dna_index/` before `dna/` whatever the
protocol, keeps a `dna/` file that is already BGZF, reads GRCh37 from `pub/grch37/release-<N>/`
rather than release 75, and writes a bgzip-compressed file with `.fai` and `.gzi` only, never
the plain `.fa` with a `Bio::DB::Fasta` index `INSTALL.pl` writes without htslib.

The one-command install is on the [Quick start](../../docs/getting-started/quick-start.md) page;
[docs/cache-setup.md](../../docs/cache-setup.md) explains what each input decides and gives
the GRCh37 command from local copies. Progress is logged to standard error; a line the stream
refuses, on a device that refuses every write or a pipe whose reader has gone, is dropped, and
the build goes on and keeps its exit status.

## Run it

```bash
cargo build --release -p vep-cache-builder
target/release/vep-cache-builder --help
cargo test -p vep-cache-builder
```

The unit tests run in-process with no network and no MySQL server: URL construction, the
download destinations and a build's work directory on a local copy of the site's tree, the
lookup of a genome over two sites of which one cannot be read (a directory whose listing is
no text: a GRCh37 tree and a release tree with no species list install from the first site,
and a species on neither is refused naming the site that could not be read), GFF3
attribute decoding, the gene tree of a RefSeq-shaped GFF3 and of a GTF (biotype rules,
identifier fallbacks, parent lists, stop codons), the translated peptide, the table dump
joins on a hand-written slice, mapper arithmetic, GTF
tag selection, VCF parsing on temporary files and region binning. The integration test
`tests/native_cache_matches_ensembl_cache.rs` runs the binary on
`tests/fixtures/grch38_chr21/`, a release 116 slice of 27 chromosome 21 genes (the GFF3
lines, three GTF transcript lines, peptides and dump rows), and compares the identifier,
symbol, xref, attribute, protein feature and exon coordinate fields of the 63 transcripts
it writes with the same transcripts of the Ensembl VEP 116 GRCh38 cache. The integration test
`tests/variation_cache_from_dumps.rs` runs the binary with `--include-variations` and
`--variation-dumps` on indexed VCFs and dump rows written to a temporary directory, checks
the `all_vars.gz` lines and `info.json` it writes, and reads the records back through `vep`'s
variation source. `tests/auto_install.rs` lays the fixture out as the Ensembl site in a
temporary directory, with a `CHECKSUMS` per directory, and runs `--auto` against it: the
installed layout, the FASTA indexes byte for byte against the committed `samtools faidx`
ones, the shards against a build named one flag at a time, `vep` resolving the result
through `--dir_cache`, the variation files, a checksum mismatch and an existing release
directory. A `dna_index/` file published beside the `.gzi` of another compression of its
sequence installs with one warning, both indexes written from the file.
The same test installs a release of another species, laid out with its core
tables, its one genome-wide variation VCF and no variation database: its variation files
cover the sequences that file holds records on, and every transcript carries
`gene_phenotype` 0; with `--chromosomes` they cover the named sequences alone, a site
lacking that variation VCF installs without known variants and says so in its last lines,
a `homo_sapiens` listing lacking a per-chromosome file `--chromosomes` names is refused
before any download, and a failure that leaves no file in the work directory names none. It also serves the tree from
a loopback HTTP server that answers 403 for what it does not hold: the release installs
byte for byte as from a server answering 404, a directory a tree may lack reads as absent,
a tree without the core database is refused as under a 404, and a 403 for a listing the
installation needs is refused naming the file and the answer.
`tests/auto_install_further_sets.rs` lays out an Ensembl Genomes
division, a collection and a rapid-release assembly and runs `--auto` against each: the
division installs, with the one variation VCF its tree publishes or without known variants
where it publishes none, and the other two are refused by name.
`tests/auto_install_default_parts.rs` lays out a tree that publishes every part an
installation brings (the variation VCFs with the 1000 Genomes file, the four table dumps
of the matrices, the regulatory-feature and motif-feature GFF3) and the same tree without
one part at a time: the cache holds each part, a `vep --everything` run on it prints no
`disabling SIFT`, no `disabling PolyPhen` and no regulatory warning, and the last lines of
each installation are the lines `close` writes for its tree, whole.
Further trees publish the files of the cells (the activity file and the funcgen tables, the
variation tables and the failed-variants GVF, the `variation_feature` table of a species with
one variation VCF) or lack one at a time: the cache holds each cell its tree publishes, its
`info.json` names the others under `left_out`, the last lines name them with the listing that
lacks the file, and a `vep` run prints the cells and keeps a failed variant under `--failed 1`.
Two more trees hold their one variation record on a sequence the model names no transcript
on, and the last lines of each installation name the sequence the cache holds no known
variant on.
`tests/predictions_species_limit.rs` runs
`--include-predictions` for a species other than human, as a build from named inputs and
as an installation: each is refused before any input is resolved.
`tests/gtf_transcript_model.rs` builds a cache from a GTF and a genome
FASTA alone. `tests/custom_transcript_source.rs` builds one gene model as an Ensembl GFF3,
as a GTF and as a generic GFF3 and checks that `--transcript_filter` and
`--exclude_predicted` keep the same transcripts on each. `tests/build_that_ends_early.rs`
holds what a build leaves in `--output-dir`: a build that fails after its shards leaves a
directory `vep` refuses by name, a finished one the list of its files, and a file the
directory held is kept and listed. `tests/malformed_inputs.rs` runs
the binary on the files under `tests/fixtures/malformed/` and on written models: a line
whose columns NUL bytes end builds the cache of the file without it, a NUL byte inside a
column of a record is read where the index build ends it, and a CDS phase Ensembl VEP
ends its run on is refused. Exons mixing ranked and unranked
lines build in start order; a seqid naming a directory outside the cache, or a FASTA index
naming more bases than its file holds, is refused with one `ERROR: [cache-builder]` line
and leaves the output directory empty. `tests/refused_models.rs` runs the binary on a
transcript whose translation starts or ends in no exon (the CDS line that sets the end
starts before the exon it overlaps, ends past it or lies in none), as a GFF3 and as a
GTF on either strand, on exon ranks that run against the strand, on a translation offset
outside its exon (ranks that run against position at the translation's end, an exon nested
in the last coding one), on a reverse-strand translation that ends in an exon nested in
the exon it starts in, on exon ranks that run against position with both ends of the
translation inside their exons, and on a model whose every transcript a build of the
cache's content sets aside: each is refused with one such line and no output directory. It also
builds the models whose CDS lines leave the
exons between the translation's two ends, one more CDS line in an intron and one CDS
line across it, to the cache the sound model builds.
`tests/build_log.rs` holds the per-chromosome log lines to name order.
`tests/log_lines.rs` holds the whole log of a one-gene build, sent to a file and to a pipe,
to plain text: a log line carries colour codes only where standard error is a terminal and
`NO_COLOR` is unset or empty.
It also holds each line that counts to its nouns, the singular beside a count of one and
the plural beside any other: the whole log of the fixture's build with its core dumps, and
the regulatory, skipped-transcript and GTF tag lines of written models;
`tests/build_log.rs` holds the lines of an eight-chromosome model so, and
`tests/variation_cache_from_dumps.rs` the variation lines.
`tests/sequence_names.rs` builds the synthetic one-gene model of
`tests/fixtures/mito_codon_table/` with its sequence renamed: one named `chromosome_1`
keeps its name, a mitochondrion written `chrM` is kept as `chrM` only where the genome
FASTA names its own otherwise (`Mito`) or, in a build for GRCh37, gives `MT` a record of
another length than the assembly's, the region kept being the one a transcript's own line
names, and a region the FASTA alone shortens takes its `--codon-table` value under either
name. A build for GRCh37 over a record that holds the model's bases two positions on is
refused with the one error, under either mitochondrial name and either length, and the
same inputs build on another assembly and, on GRCh37, where the model's own line or a
GTF beside it tags the transcript `cds_start_NF`.
`tests/protein_fasta_identifiers.rs` builds `tests/fixtures/tair10_slice/`, the GFF3 lines
and protein FASTA records of two Arabidopsis thaliana genes of Ensembl Plants release 63
whose identifiers end in a dot and a number of their own, and checks that each of the seven
transcripts carries its own record's peptide and no translation version; on written models
it checks that transcripts sharing a first component take their own records and that a
versioned token goes to the model of that version alone.
`tests/gtf_tags_identifiers.rs` builds written GFF3 models beside a GTF read for its tags
and checks that sibling transcripts keep the tags of their own lines and that a line tags
the model of its version alone. `tests/protein_fasta_record_length.rs` builds one
transcript with records of four lengths, with a genome and without one, and checks which
are its peptide, what the cache carries for the others and the warning's count; a record
that ends in the stop character is the peptide without it, and the warning counts one
record in the singular. On Linux
the binary
links OpenSSL (`reqwest` uses
`native-tls` without ALPN; the MySQL client is pure Rust), so the build needs the OpenSSL
development package; macOS links Security.framework. The workspace `rust-version` in
`Cargo.toml` is set by the `noodles` crates; this crate's `mysql` dependency chain is on the
next tier down.
`scripts/validation/validate_builder_cache.sh` builds the
binary, builds a one-region cache and compares it with captured fixtures
(`--fixture-mode`) or an Ensembl VEP-derived JSON cache (`--baseline-cache`). The
[contributing guide](../../CONTRIBUTING.md) has the repository's standard build, test,
lint and format commands.

## Related

- [vep-fasta](../vep-fasta/README.md): `IndexedFasta`, the `.fa` plus `.fai` reader behind `--genome-fasta`.
- [vep-core](../vep-core/README.md): `codon::reverse_complement` for minus-strand CDS; `prediction.rs` decodes the matrices this crate embeds; `transcript.rs` holds the runtime `Transcript` the cache fills.
- [vep-cli](../vep-cli/README.md): `crates/vep-cli/src/json_cache.rs` loads the shards named by `--json_cache`: a bare array of transcripts, which this crate and `tools/vep-storable-to-json` each write, or a map of such arrays.
- [vep-cache-converter](../vep-cache-converter/README.md): the other route to a cache, validating an Ensembl VEP-converted JSON cache and copying it into the runtime layout.
- [docs/cache-setup.md](../../docs/cache-setup.md): what each input decides and the GRCh37 command; [docs/cache-variation-data.md](../../docs/cache-variation-data.md): the variation files a build reads and the records it writes; [docs/cache-from-gtf-or-gff3.md](../../docs/cache-from-gtf-or-gff3.md): a GTF or any GFF3 as the transcript model, built into a cache when it is a run's only model or read beside a cache in the same run (`vep --gff`, `vep --gtf`); [docs/cache-predictions.md](../../docs/cache-predictions.md): the SIFT and PolyPhen predictions.
- [docs/cli-reference.md](../../docs/cli-reference.md): the `vep` flags that read what this crate writes: `--json_cache`, `--sift`, `--polyphen`, `--check_existing`, `--af`, `--af_1kg` and `--max_af`.
- [Root README](../../README.md#crates): the index of every crate.
