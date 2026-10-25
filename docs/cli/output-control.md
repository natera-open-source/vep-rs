# Output control

These flags choose which annotation fields appear in the output. `--everything` switches on exactly Ensembl VEP's `--everything` set, listed under [option sets](../cli-reference.md#option-sets-configpm351-593).

## Flags

| Flag | Type | Default | Description |
| --- | --- | --- | --- |
| `--everything`, `-e` | bool | false | Shorthand for Ensembl VEP's `--everything` set: `--sift` b, `--polyphen` b, `--ccds`, `--hgvs`, `--symbol`, `--numbers`, `--domains`, `--regulatory`, `--canonical`, `--protein`, `--biotype`, `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--uniprot`, `--mane`, `--tsl`, `--appris`, `--variant_class`, `--gene_phenotype` and `--mirna`, each overriding an explicit value as in Ensembl VEP. Every column of Ensembl VEP's `--everything` except the regulatory columns: `--regulatory` is switched on and annotates nothing, so no RegulatoryFeature or MotifFeature row is written |
| `--symbol` | bool | false | Add gene symbol to output |
| `--biotype` | bool | false | Add transcript biotype to output |
| `--canonical` | bool | false | Flag canonical transcripts |
| `--mane` | bool | false | Add the MANE, MANE\_SELECT and MANE\_PLUS\_CLINICAL columns: the MANE sets the transcript belongs to and the RefSeq accession each set names (switches on `--mane_select`) |
| `--mane_select` | bool | false | Add the MANE and MANE\_SELECT columns alone: MANE lists MANE\_Select on a MANE Select transcript and MANE\_SELECT carries its RefSeq accession |
| `--tsl` | bool | false | Add the TSL column: the digit of the transcript's support level, empty for tslNA (the GRCh37 gene set carries no level) |
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
| `--hgvsg_use_accession` | bool | false | Name the chromosome in HGVSg by its RefSeq accession: the first NC\_ name among the input name's own synonyms in the chromosome synonyms table (one hop: only a line pairing the name as written with the accession renames it); an input name with no NC\_ synonym keeps its own name, where Ensembl VEP stops with an error |
| `--spdi` | bool | false | Add the SPDI column: each allele in Sequence Position Deletion Insertion notation, the position 0-based, the deleted bases (their count when more than 20) and the inserted bases, as the input wrote the record; the chromosome is named as the cache names it, and on a row overlapping no transcript as the FASTA names it. Ensembl VEP needs a FASTA for those rows; this reads none. A record whose reference allele ends in a character outside `ACGTN-` (`AR/CA`) is annotated with no SPDI, where Ensembl VEP stops the run on it |
| `--ga4gh_vrs` | bool | false | Add the GA4GH\_VRS column: the allele's SPDI as a GA4GH Variation Representation Specification Allele object on the sequence refseq:\<NC\_ accession\>, the accession being the first NC\_ synonym of the cache's chromosome name in the chromosome synonyms table (Ensembl VEP writes the object only where its SPDI is named by an NC\_ accession: with a database, or on a row overlapping no transcript whose FASTA record carries that name; a chromosome with no NC\_ synonym gets none); JSON-encoded in the text outputs, nested in the JSON output. Cannot be combined with `--vcf` |
| `--hgvsp_use_prediction` | bool | false | Write HGVSp in the predicted form, the description after p. in parentheses (ENSP00000286808.3:p.(Arg97Ter)); the = of a synonymous description is still escaped as %3D in the text formats |
| `--remove_hgvsp_version` | bool | false | Leave the translation version off the HGVSp reference name (ENSP00000286808:p.Arg97Ter); HGVSc keeps its version |
| `--sift` | `s`, `p`, `b` | none | Add the SIFT column on every single-residue substitution inside the coding sequence (start\_lost included, a stop change never): p the prediction, s the score, b both as prediction(score). Needs a cache whose metadata names a SIFT version; a cache without one stops the run, or under `--everything` or `--safe` drops the flag with an INFO line (a native cache carries the matrices only when built with `--include-predictions`) |
| `--polyphen` | `s`, `p`, `b` | none | Add the PolyPhen column on every single-residue substitution inside the coding sequence (start\_lost included, a stop change never): p the prediction, s the score, b both as prediction(score). Needs a cache whose metadata names a PolyPhen version; a cache without one stops the run, or under `--everything` or `--safe` drops the flag with an INFO line (a native cache carries the matrices only when built with `--include-predictions`) |
| `--humdiv` | bool | false | Read the PolyPhen column from the HumDiv-trained matrix instead of the default HumVar one; a cache without the HumDiv matrix prints an empty column with no error, as in Ensembl VEP |
| `--condel` | String | none | Accepted without effect: Ensembl VEP declares it and reads nothing from it |
| `--nearest` | `transcript`, `gene`, `symbol` | none | Add the NEAREST column: the protein-coding transcription start site(s) nearest each variant as a transcript stable id, gene stable id or gene symbol; the site above the variant and the site below are measured to the nearer end of the variant, both are reported when they tie, and every transcript starting at a reported site is listed in stable-id order where Ensembl VEP lists one |
| `--numbers` | bool | false | Add exon/intron numbers |
| `--domains` | bool | false | Add the DOMAINS column: on a coding row, the protein features whose residues the variant overlaps, each as label:identifier (Pfam:PF00080, Gene3D:1.10.510.10; a native cache carries the features only when built with `--core-dumps`) |
| `--regulatory` | bool | false | Accepted for VEP compatibility; regulatory-feature annotation is not implemented and the flag warns |
| `--regulatory_gff` | String | none | Refused with an error: `--regulatory_gff` is not supported: it would annotate regulatory features and motifs from a GFF3 file in place of the cache, and regulatory annotation is not supported |
| `--extended_promoters` | bool | false | Refused with an error: `--extended_promoters` is not supported: it would widen the promoters read with `--regulatory_gff` to their extended bounds, and regulatory annotation is not supported |
| `--variant_class` | bool | false | Add variant class from SO |
| `--overlaps` | bool | false | Add the OverlapBP and OverlapPC columns to the tab, VCF and Parquet outputs: the bases of a structural variant inside the transcript and the percentage of the transcript they cover (transcript rows only; regulatory features are not annotated) |
| `--check_existing` | bool | false | Check for co-located known variants |
| `--no_check_alleles` | bool | false | Report a known variant as co-located when its start and end equal the input variant's, whatever its alleles; such a variant prints on every allele's row |
| `--exclude_null_alleles` | bool | false | Never report a known variant whose alleles the cache does not hold (COSMIC and HGMD records, dbSNP records without variation); without the flag such a variant is reported when its start and end equal the input variant's |
| `--failed` | i64 | 0 | Report known variants whose failed value is at most this: 0 excludes variants that failed Ensembl's quality checks, 1 includes them |
| `--clin_sig_allele` | i64 | 1 | Report the ClinVar significance of the row's own allele where a co-located variant carries allele-specific entries (1, the default); 0 reports every variant's significance whatever the allele |
| `--clinvar_somatic_classification` | bool | false | Add the ClinVar somatic classification of co-located variants (CLINVAR\_SOMATIC\_CLASSIFICATION), the cache's string as written |
| `--var_synonyms` | bool | false | Add the known synonyms of co-located variants and their sources (VAR\_SYNONYMS); switches on `--check_existing` |
| `--check_svs` | String | none | Report the known structural variants overlapping each variant (SV) from Ensembl's indexed structural-variant VCF (homo\_sapiens\_structural\_variations.vcf.gz with its .csi or .tbi beside it), Ensembl's dump of the database Ensembl VEP reads them from; Ensembl VEP takes the flag without a value and refuses it offline |
| `--coding_only` | bool | false | Keep only the rows whose allele lies in a transcript's coding sequence (a CDS position or a coding consequence term, a variant covering a translated transcript whole included); intergenic and non-coding rows go, and a record left with no row is written as an unannotated one |
| `--no_intergenic` | bool | false | Write no row for a record that overlaps no transcript: no intergenic row in the default and tab formats, the VCF line without a CSQ entry, and in JSON a most\_severe\_consequence of "?" with no intergenic\_consequences |
| `--hgnc` | bool | false | Accepted without effect: Ensembl VEP declares it and reads nothing from it |

