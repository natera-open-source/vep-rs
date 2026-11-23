# Output control

These flags choose which annotation fields appear in the output, report the known variants co-located with each input variant, and drop the rows outside a coding sequence or outside every transcript. `--everything` switches on exactly Ensembl VEP's `--everything` set, listed under [option sets](../cli-reference.md#option-sets-configpm351-593).

## Annotation fields

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--everything`, `-e` | bool | false | Shorthand for Ensembl VEP's `--everything` set: `--sift` b, `--polyphen` b, `--ccds`, `--hgvs`, `--symbol`, `--numbers`, `--domains`, `--regulatory`, `--canonical`, `--protein`, `--biotype`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--uniprot`, `--mane`, `--tsl`, `--appris`, `--variant_class`, `--gene_phenotype` and `--mirna`, each overriding an explicit value as in Ensembl VEP. `--regulatory` in the set writes the regulatory and motif rows of a cache that carries regulatory features, every motif field included where the cache carries the binding matrices' counts; CELL\_TYPE is added by `--cell_type` beside the set where it carries cell-type activity |
| `--symbol` | bool | false | Add gene symbol to output |
| `--biotype` | bool | false | Add transcript biotype to output |
| `--canonical` | bool | false | Flag canonical transcripts |
| `--mane` | bool | false | Add the MANE, MANE\_SELECT and MANE\_PLUS\_CLINICAL columns: the MANE sets the transcript belongs to and the RefSeq accession each set names (switches on `--mane_select`) |
| `--mane_select` | bool | false | Add the MANE and MANE\_SELECT columns alone: MANE lists MANE\_Select on a MANE Select transcript and MANE\_SELECT carries its RefSeq accession |
| `--tsl` | bool | false | Add the TSL column: the digits of the transcript's support level as its TSL attribute spells them, empty for tslNA (the GRCh37 gene set carries no level) |
| `--appris` | bool | false | Add the APPRIS column with the grade abbreviated as Ensembl VEP prints it: principal1 is P1, alternative2 is A2, principalM is PM |
| `--ccds` | bool | false | Add the CCDS column: the transcript's CCDS identifier (a native cache carries it from the GFF3; no RefSeq or merged cache is built natively) |
| `--protein` | bool | false | Add the ENSP column: the Ensembl protein identifier |
| `--protein_version` | bool | false | Add the translation version to the Ensembl protein ID (ENSP00000286808.3); needs `--protein` |
| `--uniprot` | bool | false | Add the SWISSPROT, TREMBL, UNIPARC and UNIPROT\_ISOFORM columns: the translation's UniProt cross-references (a native cache carries them only when built with `--core-dumps`) |
| `--xref_refseq` | bool | false | Add the RefSeq column: the RefSeq mRNA identifiers cross-referenced to the transcript (a native cache carries them only when built with `--core-dumps`) |
| `--transcript_version` | bool | false | Append the transcript's version to the Feature column (ENST00000366667.4) unless the identifier already ends in one |
| `--gene_version` | bool | false | Append the gene's version to the Gene column (ENSG00000142168.14) unless the identifier already ends in one |
| `--hgvs` | bool | false | Add HGVS coding and protein notations; needs a FASTA named with `--fasta` or found in the cache directory and stops the run without one |
| `--hgvsg` | bool | false | Add the HGVS genomic notation (HGVSg); needs a FASTA named with `--fasta` or found in the cache directory and stops the run without one |
| `--hgvsg_use_accession` | bool | false | Name the chromosome in HGVSg by its RefSeq accession: the first NC\_ name among the input name's own synonyms in the chromosome synonyms table, or else among the synonyms of the chromosome the name resolves to, the accession Ensembl VEP prints with a database for a human chromosome; on an input name with no NC\_ synonym of its own Ensembl VEP stops with an error offline. A scaffold or patch with no NC\_ synonym takes its RefSeq name of another prefix (NT\_, NW\_) the same way, and a sequence the table gives no RefSeq name takes the cache's name for it. An accession the table lists in more than one form is printed in its latest: the highest version, a versioned form over one without a version |
| `--spdi` | bool | false | Add the SPDI column: each allele in Sequence Position Deletion Insertion notation, the position 0-based, the deleted bases (their count when more than 20) and the inserted bases, as the input wrote the record; the chromosome is named as the cache names it, and on a row overlapping no transcript as the FASTA names it. Ensembl VEP needs a FASTA for those rows; this reads none. A record whose reference allele ends in a character outside `ACGTN-` (`AR/CA`) is annotated with no SPDI, where Ensembl VEP stops the run on it |
| `--ga4gh_vrs` | bool | false | Add the GA4GH\_VRS column: the allele's SPDI as a GA4GH Variation Representation Specification Allele object on the sequence refseq:\<NC\_ accession\>, the accession being the first NC\_ synonym of the cache's chromosome name in the chromosome synonyms table, in its latest form where the table lists it in more than one (Ensembl VEP writes the object only where its SPDI is named by an NC\_ accession: with a database, or on a row overlapping no transcript whose FASTA record carries that name; a chromosome with no NC\_ synonym gets none); JSON-encoded in the text outputs, nested in the JSON output. Cannot be combined with `--vcf` |
| `--hgvsp_use_prediction` | bool | false | Write HGVSp in the predicted form, the description after p. in parentheses (ENSP00000286808.3:p.(Arg97Ter)); the = of a synonymous description is still escaped as %3D in the text formats |
| `--remove_hgvsp_version` | bool | false | Leave the translation version off the HGVSp reference name (ENSP00000286808:p.Arg97Ter); HGVSc keeps its version |
| `--sift` | `s`, `p`, `b` | none | Add the SIFT column on every single-residue substitution inside the coding sequence (start\_lost included, a stop change never): p the prediction, s the score, b both as prediction(score). Needs a cache whose metadata names a SIFT version; a cache without one stops the run, or under `--everything` or `--safe` drops the flag with an INFO line (a native cache carries the matrices when its build read them: an `--auto` installation where the release publishes them, or `--include-predictions`) |
| `--polyphen` | `s`, `p`, `b` | none | Add the PolyPhen column on every single-residue substitution inside the coding sequence (start\_lost included, a stop change never): p the prediction, s the score, b both as prediction(score). Needs a cache whose metadata names a PolyPhen version; a cache without one stops the run, or under `--everything` or `--safe` drops the flag with an INFO line (a native cache carries the matrices when its build read them: an `--auto` installation where the release publishes them, or `--include-predictions`) |
| `--humdiv` | bool | false | Read the PolyPhen column from the HumDiv-trained matrix instead of the default HumVar one; a cache without the HumDiv matrix prints an empty column with no error, as in Ensembl VEP |
| `--condel` | String | none | Accepted without effect: Ensembl VEP declares it and reads nothing from it |
| `--nearest` | `transcript`, `gene`, `symbol` | none | Add the NEAREST column: the protein-coding transcription start site(s) nearest each variant as a transcript stable id, gene stable id or gene symbol; the site above the variant and the site below are measured to the nearer end of the variant, both are reported when they tie, and every transcript starting at a reported site is listed in stable-id order where Ensembl VEP lists one |
| `--numbers` | bool | false | Add exon/intron numbers |
| `--domains` | bool | false | Add the DOMAINS column: on a coding row, the protein features whose residues the variant overlaps, each as label:identifier (Pfam:PF00080, Gene3D:1.10.510.10; a native cache carries the features only when built with `--core-dumps`) |
| `--regulatory` | bool | false | Annotate overlaps with cached regulatory and motif features; warns when the cache has no regulatory data |
| `--cell_type` | String (list of the cache's epigenome names) | none | Report the activity of each regulatory feature in the listed epigenomes, the CELL\_TYPE field (`<epigenome>:<activity>` per listed epigenome the feature has a record for; empty on a motif row); switches on `--regulatory`. The list is read against the cache's epigenome names, which may contain commas (`A549_(m,_58_y),BJ_(m)` names two); a name the cache lacks is an error listing the names |
| `--regulatory_gff` | String | none | Use regulatory annotation from Ensembl regulatory GFF3 files, instead of data from a cache; overrides the cache's regulatory data when both are present. A comma-separated list of key=value pairs, in the same style as `--custom`, providing one or more of: file (the regulatory-feature GFF3), motifs (a transcription factor binding motif GFF3) and emars (an epigenetically modified accessible region GFF3); a bare file name is the file. Each file is bgzip-compressed and tabix-indexed. Features from file and emars have a Feature type of RegulatoryFeature with the feature's type under BIOTYPE; features from motifs have a Feature type of MotifFeature. HIGH\_INF\_POS and MOTIF\_SCORE\_CHANGE require the position weight matrix, which the GFF does not carry, and are omitted. `--cell_type` is ignored, with a warning, as a GFF carries no epigenome activity data |
| `--extended_promoters` | bool | false | Use the extended boundaries of promoter features when annotating with `--regulatory_gff`: the extended\_start and extended\_end attributes, which are never inside the core region, so the flag only ever adds annotation. Affects promoter features only, and changes nothing without `--regulatory_gff` |
| `--variant_class` | bool | false | Add variant class from SO |
| `--overlaps` | bool | false | Add the OverlapBP and OverlapPC columns to the tab, VCF and Parquet outputs: the bases of a structural variant inside the transcript and the percentage of the transcript they cover (also reported for cached regulatory and motif features) |
| `--hgnc` | bool | false | Accepted without effect: Ensembl VEP declares it and reads nothing from it |
| `--mirna` | bool | false | Add the miRNA column: the secondary-structure elements of a miRNA transcript the variant overlaps, as the SO terms miRNA\_stem and miRNA\_loop (a native cache carries the structure only when built with `--core-dumps`) |

## Co-located variants

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--check_existing` | bool | false | Check for co-located known variants |
| `--no_check_alleles` | bool | false | Report a known variant as co-located when its start and end equal the input variant's, whatever its alleles; such a variant prints on every allele's row |
| `--exclude_null_alleles` | bool | false | Never report a known variant whose alleles the cache does not hold (COSMIC and HGMD records, dbSNP records without variation); without the flag such a variant is reported when its start and end equal the input variant's |
| `--failed` | i64 | 0 | Report known variants whose failed value is at most this: 0 excludes variants that failed Ensembl's quality checks, 1 includes them |
| `--clin_sig_allele` | i64 | 1 | Report the ClinVar significance of the row's own allele where a co-located variant carries allele-specific entries (1, the default); 0 reports every variant's significance whatever the allele |
| `--clinvar_somatic_classification` | bool | false | Add the ClinVar somatic classification of co-located variants (CLINVAR\_SOMATIC\_CLASSIFICATION), the cache's string as written |
| `--var_synonyms` | bool | false | Add the known synonyms of co-located variants and their sources (VAR\_SYNONYMS); switches on `--check_existing` |
| `--check_svs` | String | none | Report the known structural variants overlapping each variant (SV) from Ensembl's indexed structural-variant VCF of the species (\<species\>\_structural\_variations.vcf.gz, for human homo\_sapiens\_structural\_variations.vcf.gz, with its .csi or .tbi beside it), Ensembl's dump of the database Ensembl VEP reads them from; the file is read under whatever name it is given; Ensembl VEP takes the flag without a value and refuses it offline |

## Dropped rows

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--coding_only` | bool | false | Keep only the rows whose allele lies in a transcript's coding sequence (a CDS position or a coding consequence term, a variant covering a translated transcript whole included); intergenic and non-coding rows go, and a record left with no row is written as an unannotated one |
| `--no_intergenic` | bool | false | Write no row for a record that overlaps no transcript: no intergenic row in the default and tab formats, the VCF line without a CSQ entry, and in JSON a most\_severe\_consequence of "?" with no intergenic\_consequences |

