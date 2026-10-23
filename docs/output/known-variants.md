# Co-located variants, frequencies and known structural variants

The keys on this page are the co-located variant keys, the frequency keys, `FREQS` and `SV`. They appear in the default format's Extra column, as `--tab` columns and as CSQ subfields of the VCF output, under the flags listed on the [allele frequencies and co-located variants](../cli/allele-frequencies.md) page. Their place among the Extra keys is given on the [VEP default and tab](vep-and-tab.md#extra-field) page.

## Co-located variant keys

The co-located variant keys follow Ensembl VEP's `add_colocated_variant_info`. A row lists the known variants matched on its own allele or on their span alone, in Ensembl VEP's order: germline before somatic, then dbSNP, HGMD, COSMIC, other. `SOMATIC` and `PHENO` carry one `0` or `1` per listed variant and are left out when every value is `0`. `PUBMED` lists every id of every variant, duplicates kept. `VAR_SYNONYMS` (under `--var_synonyms`, which switches `--check_existing` on) joins each variant's `source::id,id` string with `--`.

A multi-allelic record is matched whole, so a known variant matched on any of its alleles reaches every allele's rows and prints on the rows of the alleles it matched. Under `--minimal`, a record whose alleles trim to different spans is matched as Ensembl VEP matches it after splitting the record. The known variants of the first alternate allele, trimmed and on its trimmed span, are the record's, and every row keeps those matched on its own allele or on their span alone. A later allele's own known variants are not listed, as Ensembl VEP builds every row of the record from the first allele's feature, unless a `--check_frequency` drop [writes each allele as a record of its own](#filtering-by-frequency).

`CLINVAR_SOMATIC_CLASSIFICATION` (under `--clinvar_somatic_classification`) is the cache's classification string as written, `phenotype(classification)` entries joined by `;`, the last listed variant's. Ensembl VEP builds that string from ClinVar's somatic records when it dumps a cache or reads its database, and vep-rs reads a cache only, so a classification the cache lacks is never recovered.

`CLIN_SIG` reports, by default (`--clin_sig_allele 1`), the ClinVar entries recorded for the row's own allele when any listed variant carries one, `;`-joined across variants, each variant's terms `,`-joined. Where a variant's ClinVar reference allele differs from the reference `REF_ALLELE` prints (the input's, or under `--minimal` a split record's first allele's trimmed reference), its entry alleles are read reverse-complemented first. When no listed variant carries an entry for the row's allele, it lists the terms of the variants read before the first one with allele-specific entries, duplicates kept. Under `--clin_sig_allele 0` it lists every listed variant's terms, duplicates kept.

## Frequency keys

The frequency keys follow Ensembl VEP's `add_colocated_frequency_data`. Each flag in the table adds its columns and switches `--check_existing` on.

| Flag | Columns |
| --- | --- |
| `--af` | `AF` |
| `--af_1kg` | the 1000 Genomes columns: `AF`, `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, `SAS_AF` |
| `--af_gnomade`, or its older spelling `--af_gnomad` (one or the other, never both) | the gnomAD exome columns: `gnomADe_AF` and nine populations |
| `--af_gnomadg` | the gnomAD genome columns: `gnomADg_AF` and ten populations, `gnomADg_AMI_AF` included |

A row reports, for each listed variant, the frequency of the record allele matched to the row's own allele, as the cache wrote it (`0.005244`, `6.875e-05`), except `AF`, which prints with four decimals (`0.0068`). The values of every listed variant merge into one `,`-joined list, each value once. A biallelic record whose `AF` stores one allele yields the remainder `1 - total` for the other (`0.9932`). A variant matched on its span alone is read only under `--af` or `--max_af` and contributes only in that case, the `AF` entry of the row's own allele, stored or remainder.

`--max_af` adds `MAX_AF`, the largest frequency over every population column of every group whether its flag is set or not, leaving out the combined populations (`AF`, `gnomADe_AF`, `gnomADg_AF`). It also adds `MAX_AF_POPS`, the populations reaching it, in the order the groups sort (`af`, `af_1kg`, `af_gnomad`, `af_gnomade`, `af_gnomadg`), and each group lists its populations. A variant whose populations are all `0` writes `MAX_AF=0` with every zero population listed, and across variants an equal maximum appends its populations as they are.

`--af_gnomad` on a cache whose declared variation columns lack `gnomADe` ends the run with `ERROR: [cache] gnomAD data is not available in this cache`, Ensembl VEP's message (`CacheDir.pm` line 195) under the `cache` component. The check runs before the FASTA and the input open and under `--show_cache_info` too, reads that spelling alone and is skipped under `--everything`, which switches `--af_gnomade` on and is refused beside `--af_gnomad` unless `--safe`; so `--af_gnomade` never checks. Under `--old_maf` the `AFR_AF`, `AMR_AF` and `EUR_AF` columns print empty on a release-116 cache, as they do in Ensembl VEP, whose rewrite of those columns leaves no allele to match.

## Filtering by frequency

`--check_frequency` filters the input before any transcript is read, as Ensembl VEP's variation source filters its buffer. `--filter_common` switches it on with `--freq_pop 1KG_ALL --freq_freq 0.01 --freq_gt_lt gt --freq_filter exclude`. For each record, the population `--freq_pop` names is read from every listed variant: the column of that name, else the column named by the second of its `_`-separated parts, so `gnomADe_AFR` falls back to the 1000 Genomes `AFR`.

A variant on the other strand is read through its reverse complement, and the one allele a variant names no frequency for takes the remainder `1 - total`. Each alternate allele takes its own value, else the value of the record allele it was matched to, else `NA`. It passes when the value is at or above `--freq_freq` under `gt` or at or below it under `lt`; a value of `0` is no value.

Under `exclude` a record whose every allele passes is dropped, one with some passing alleles keeps the failing ones alone, and one with no data or no passing allele is kept whole. Under `include` a record is kept only when every allele passes, so one without data is dropped.

So is a non-variant record written through under `--allow_non_variant`, unless a known variant lying exactly at its span gives it data. That variant is a null-allele variant (`COSMIC_MUTATION`, `HGMD_MUTATION`), whose one allele takes the remainder `1`, except under `--exclude_null_alleles`, which matches none, or any variant with data under `--no_check_alleles`. The record is then kept, as Ensembl VEP keeps it, since it has no allele to fail. A non-variant record with a true INFO/SVTYPE (not empty, not `0`) is dropped under `include` whatever lies at its span, as it is in Ensembl VEP, which reads it as a structural variant and co-locates no known variant with one.

A record that loses alleles is renumbered (`ALLELE_NUM`) and its allele string rewritten to the kept alleles. The JSON `allele_string` and `AMBIGUITY` read the rewritten one, and so does `UPLOADED_ALLELE` unless the run is `--minimal` or the alleles differ in length with no `-` among them, when Ensembl VEP prints the uploaded spelling. `FREQS` lists, on every row of a record that had data, `POPULATION:ALLELE:frequency` for every alternate allele as uploaded, dropped alleles included, the population as `--freq_pop` was given in upper case (`AF:C:0.006789`, `GNOMADE_AFR:-:0.0001206`, `AF:A:NA`). The values come in the alleles' order, where Ensembl VEP writes them in hash order; the JSON format has no `FREQS`.

Under `--minimal` a multi-allelic record whose alleles trim is judged allele by allele, each on its own span with its own co-located variants and one frequency. Ensembl VEP's input buffer splits such a record into one feature per allele before the variation source runs, and each feature is judged alone. A batch in which the filter dropped no record and no split feature writes the record as one, every row carrying the first allele's `FREQS`, since the rejoin builds the rows from the first feature.

A batch in which a record or a split feature was dropped writes each surviving allele of a split record as a record of its own. Each such record has its own span, co-located variants and `FREQS`, `MINIMISED` on the first alone, its own `VARIANT_CLASS`, one VCF line and one JSON object. A split feature that writes no row (a `*`, a `<DEL:*>`, an allele equal to the reference once trimmed) counts as a drop when the filter drops it; a record `exclude` keeps with its failing alleles does not.

The first feature's `UPLOADED_ALLELE` is the record's anchor-chopped allele string, as every row's is when the record is written as one (`TG/-/AG` for `GTG/G/GAG`), and a later feature's is the uploaded spelling. Ensembl VEP's drop resets its buffer and with it the rejoin, so a reset record takes the same shape in both engines. A split record none of whose alleles writes a row (`GG *,GG`) stays one VCF line and one JSON object under a reset, where Ensembl VEP writes each of its features as its own.

The reach of the reset differs. Here it is the dropped record's run of one chromosome within a batch of `--buffer_size` input lines counted from the file's first data line, skipped lines included. In Ensembl VEP it is the buffer, one chromosome's records up to `--buffer_size` (`InputBuffer.pm` lines 191-199), or under `--fork` a chunk of the buffer holding at most `--buffer_size` over twice the fork count (`Runner.pm` line 486). Where the units differ, one engine rejoins a split record the other resets.

Under `--buffer_size 10`, take a split record on the tenth line after a line on a chromosome the cache lacks, with a dropped feature on a line between them. Ensembl VEP without `--fork` resets it, since its buffer never held the skipped line and so reaches the tenth record; here the batch ends one line before it, and the record is rejoined.

`--freq_pop` accepts `1KG_ALL`, a declared variation column in any case (`af`, `AFR`, `gnomade_nfe`) and `1KG_` before one. It refuses every other name, the gnomAD genome populations included, with `ERROR: [cache] Invalid population <name> specified`, Ensembl VEP's message (`BaseCacheVariation.pm` line 285), before the input opens. On a cache without variation data the filter is inactive, as Ensembl VEP builds no variation source there. `1KG_ALL` reads the cache's `AF` column where Ensembl VEP 116.2 finds no data and filters nothing ([intended divergences](../divergences/filter-common-1kg-all.md)). The records dropped whole are counted in the run's closing line as `filtered out`.

## SV

**SV** (under `--check_svs <file>`) lists the names of the known structural variants overlapping the record, sorted and unique, on every row of the record; the JSON key is `colocated_structural_variants`, a list at the record's top level. Ensembl VEP reads them from its variation database and refuses the flag offline. vep-rs takes the flag with a file, Ensembl's `homo_sapiens_structural_variations.vcf.gz` of the matching release with its `.csi` or `.tbi` beside it, and reads it offline (the flag is accepted under `--offline`).

That file is Ensembl's dump of the same database query, and vep-rs reads from it what the query returns. The supporting-evidence records (`Parent`) are left out, as the database query leaves them out. Records of the `probe` class are left out, as the query leaves the class out (the dump carries none). The failed and somatic structural variants are absent from the dump, as they are from the query. A record's span is read as the dump writes it, `POS + 1` to `END`, which is the database's span on every chromosome 21 record of release 116.

The dump lacks three kinds of mapping the database holds. First, the dump covers the 25 chromosomes (1 to 22, X, Y, MT) and no scaffold or patch, so a variant on an unplaced or unlocalised scaffold gets no name where the database holds a mapping. Release 116 holds 27,305 top-level features on 360 other regions, 4,903 of them on 75 `GL` and `KI` scaffolds, several of which the release-116 cache carries.

Second, a structural variant mapped to several locations is in the dump at one location per dump piece, so a further mapping the database holds may be missing. On chromosome 21 of release 116 that is 713 of the 97,616 database mappings, which lie under the same name's dump mapping except at 590 of the chromosome's 46,709,983 positions, where the database names a structural variant the dump lacks.

Third, a structural variant starting at the first base of its chromosome is absent from the dump. Release 116 has six, on the chromosomes 4, 9, 12, X and Y; two of them span the whole of X and of Y, so every variant on those two chromosomes lacks one name the database reports.

A record overlaps the input when their spans meet, an insertion between two bases overlapping the records that span both, on the input's own chromosome; a breakend allele's mate position is read as well, on the mate's chromosome. Ensembl VEP matches the records of a buffer to its variants without a chromosome check, so a bracket mate on another chromosome collects, beside the mate chromosome's variants, the buffer chromosome's variants that span the mate's number ([the divergence page](../divergences/check-svs-mate-chromosome.md)).

A record's alleles share the names, as Ensembl VEP keeps them on the record: a breakend record's alleles print the names either position reached. Under `--minimal`, a record whose alleles trim to different spans prints on every row the names overlapping the first alternate allele's trimmed span, as Ensembl VEP builds every row of the record from that allele's feature; a later allele's own span is not read.
