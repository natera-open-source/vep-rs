# Plugin fixtures

The data files a plugin corpus under `tests/golden/116/` commits beside its
records, so vep-rs and Ensembl VEP's plugin read one copy of the same bytes.
A corpus names them in its manifest under `plugin_files`, per plugin, and its
flags reach them through two placeholders both engines expand: `{corpus}`, the
corpus directory, and `{plugins}`, the plugin directory each engine reads (the
corpus's `plugin_data/` for vep-rs, the directory of Ensembl VEP's plugin
modules in the reference run). `generate_reference_output.sh --corpus` checks
that every named file exists, mounts the corpus and records each file's sha256
in `provenance.json`, so a reference is tied to the fixture it was generated
over; a regenerated fixture means a regenerated reference.

Two kinds of fixture, decided by the data's terms:

- **A cut of the real file**, for a plugin whose data may be redistributed
  (LoFTEE's ancestral-allele FASTA and GERP bigWig, CC0 and Apache-licensed
  tables). `cut_plugin_regions.py` writes it from the corpus itself: the
  contigs keep their full length, so every coordinate stays genomic, and the
  bases a run never reads are masked to `N` in a FASTA or left out of a bigWig.
  The regions are the cached transcripts' spans (the run's reference), the
  records (the ancestral alleles the plugin reads at the variant) and the
  transcripts' exons (the bases LoFTEE's END_TRUNC walk sums). The bigWig is
  written without zoom levels, so Ensembl VEP's summary read of an interval
  totals the per-base values exactly, as vep-rs's interval read does; a real
  bigWig's zoom levels are a lossy summary, and the two engines' totals
  differ on it, a difference `docs/plugins/loftee.md` documents
  (`loftee_gerp_zoom_level_summary`) and the fixture does not carry.
- **A synthetic file a generator writes**, for a plugin whose data file is not
  redistributed here: PrimateAI, ClinPred, BayesDel, SpliceVault, the three CADD
  files, dbNSFP and dbscSNV, with SpliceAI taking the same route. The generator
  writes the file's preamble lines, header and column order, with fabricated
  scores derived from a seed at the corpus's own sites. The row a plugin must
  pick for a record carries the record's allele.
  `plugin_score_fixtures.py` writes the PrimateAI,
  ClinPred, BayesDel, SpliceVault, CADD and SpliceAI files;
  `make_plugin_fixture.py` the layouts a plugin prepares per assembly (dbNSFP,
  dbscSNV), reading the corpus's own VEP-format output for the transcripts and
  amino acid changes a row must carry; and `maxentscan_fixture_matrices.py`
  writes a seeded matrix directory in the layout the MaxEntScan plugin reads,
  which a corpus names in its `generated` block so both engines read one copy.
  A layout's header order is the one its generator lists, and these generators
  read no data file of the plugins they write for.

Both kinds stay inside the corpus budget `golden.rs` asserts (15 MiB on disk
with the run's reference counted separately). A fixture is committed with the
tool versions and digests under `plugin_fixtures` in `provenance.json`, and the
corpus's `README`-level description of its records lives in `provenance.json`'s
`input` field, as every corpus's does.

| Script                   | Purpose                                                                                                                                                                                                                                      |
| ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `cut_plugin_regions.py`  | Cut a reference FASTA, an ancestral-allele FASTA or a GERP bigWig, any of the three, to the regions a corpus reads: `--corpus <dir> [--reference <fa>] [--ancestor <fa>] [--gerp <bw>] [--gerp-regions exons\|records] [--gerp-name <file>]`; writes `reference.fa.gz` and `plugin_data/{human_ancestor.fa.gz,gerp.bw}` with their indexes for the files given and records the digests; `--gerp-name` writes the bigWig under the name a plugin derives from its directory (LoFTEE's `gerp_conservation_scores.homo_sapiens.GRCh38.bw`). |
| `test_cut_plugin_regions.py` | Tests for the cut on a hand-built contig, cache and bigWig.
| `plugin_score_fixtures.py` | Write a corpus's synthetic data file for PrimateAI, ClinPred, BayesDel, SpliceVault, CADD or SpliceAI: the shape's preamble, header and column layout with seeded scores at the corpus's sites (one row per alternate base; the BayesDel shape adds insertions, deletions and longer substitutions; the `cadd_snv` and `cadd_indel` shapes write the CADD v1.7 SNV and indel files whose first line names `--assembly` and says the scores are fabricated, the `cadd_sv` shape the CADD-SV fixture with one row per `<DEL>`, `<INS>` or `<DUP>` record and two decoys; the `spliceai_snv` and `spliceai_indel` shapes write the SpliceAI score VCFs with one record per gene symbol `--symbols` lists at a site; `--miss` leaves a site unscored; `--mitochondrial-name` writes the mitochondrial contig's rows under `M`, `MT`, `chrM` or `chrMT` whatever the input spells it, for a table whose files name that contig otherwise than the input does, and moves no score; `--second-stream-at` writes the file as two bgzip streams joined, the second opening at the named site, so the first stream's end marker, a block that holds no bytes, stands inside the file; `--join-entries`, `--no-information`, `--join-alternates` and `--any-base-reference` write a site's SpliceAI records in another form of the layout (several entries in one record, a missing-value INFO, the alternate bases listed in one record, `N` as the reference allele); `--repeat-lines` writes every record line of a SpliceAI site N times, next to each other, in the plain form or beside one of those forms; `--rescore-lines` writes them N times with the scores of every entry of each further copy drawn again; `--header-names` writes the given column names in the BayesDel or ClinPred header line, for a table whose header spells or names its columns otherwise than the shape's own header; `--dir` names the corpus directory written under), bgzip-compressed, indexed with `tabix -s 1 -b 2 -e 2` (`tabix -0 -s 1 -b 2 -e 3` for the CADD-SV file and `tabix -p vcf` for the SpliceAI files, the command SpliceAI.pm:114-122 gives), and recorded with its seed, row counts, digests and tool versions, and the mitochondrial name and header names when given, under `plugin_fixtures`.
| `test_plugin_score_fixtures.py` | Tests on hand-built inputs: the allele representation, the rows per shape, determinism, the PrimateAI columns, the mitochondrial contig's name, the header names, the CADD rows and headers, the SpliceAI records, their symbols file, their other forms, the lines written more than once and the copies whose scores are drawn again, the two joined streams, the written files and their records, the refusals.                                                                                                                                                                              |
| `make_plugin_fixture.py` | Write a corpus's synthetic data file for dbNSFP or dbscSNV, prepared for one assembly as Ensembl VEP's plugin documents (dbNSFP.pm:62-74, dbscSNV.pm:49-59): `--layout dbnsfp\|dbscsnv --assembly GRCh37\|GRCh38 --seed N --vcf <variants.vcf> --reference <expected/default.txt.gz> --out <plugin_data/file.gz> [--provenance <provenance.json>] [--miss chr:pos]`; the generator's dbNSFP header (458 names; the GRCh38 file indexed on columns 1 and 2, the GRCh37 file on 8 and 9; the rows of the input's own allele carry the amino acid changes and transcripts the output shows, one row per change; a row planted in `1:1008170-1082927` with two transcripts in every per-transcript column, the region the plugin reads to find the `;`-joined columns, dbNSFP.pm:162) and the generator's dbscSNV header (18 names, with no leading `#`, every line the generator writes ending in CRLF, which dbscSNV.pm:191 strips; the GRCh37 file indexed on columns 1 and 2, the GRCh38 file on 5 and 6); the other assembly's coordinate columns carry the position plus 100,000; `<file>.spec.json` records the layout, preparation, seed, input digests and row counts, and `--provenance` records the fixture under `plugin_fixtures`. CADD and SpliceAI are registered names without a writer. |
| `test_make_plugin_fixture.py` | Tests on hand-built inputs: the header order and preparation per layout, the sites and output rows read, the dbNSFP rows per amino acid change and the planted row, the dbscSNV rows per assembly, determinism, the provenance record, the refusals. |
| `maxentscan_fixture_matrices.py` | Write a MaxEntScan matrix directory in the layout the MaxEntScan plugin reads (MaxEntScan.pm:104-105, 824-826) from a seed (`me2x5`, `splicemodels/splice5sequences`, `splicemodels/me2x3acc1..9`, the line counts the plugin's hashing implies (MaxEntScan.pm:803-820, 851-859), fabricated scores, each number in a spelling drawn from the seed: its decimal or an exponent of either sign, the lines of `me2x5` between blanks, and files the plugin never opens beside the tables; `--plain-lines` writes the bare decimals and the eleven tables alone); `--score` prints the plugin's `score5` or `score3` of a k-mer over a written directory, and `--check-kmers` refuses a directory whose k-mers or paired differences sit at a `%.3f` rounding boundary.
| `test_maxentscan_fixture_matrices.py` | Tests for the generator: line counts, value range, the spellings of a line and the files beside the tables, the shuffled 7-mer permutation, determinism, the digests, the scoring arithmetic on hand-built tables and the boundary guard.
