# Other species

Any species with an Ensembl or Ensembl Genomes cache converts and annotates the same way as human ([Cache setup](cache-setup.md), [Cache conversion](cache-conversion.md)).

In each such cache:

- the species directory is the Latin name, with a strain or an assembly accession where Ensembl's name carries one;
- each transcript's genetic code is its record's `codon_table`;
- the chromosome names are the cache's.

An Ensembl Genomes cache carries that division's own release number, so `--cache_version` names it rather than the Ensembl release.

vep-rs is verified by concordance against Ensembl VEP 116.3 on human GRCh37 and GRCh38, the whole-genome datasets of the current release, and on 35 other species, among them mouse GRCm39, fruit fly BDGP6.54, Arabidopsis thaliana TAIR10 and yeast R64-1-1.

Each of the 35 has a golden corpus of its own under `tests/golden/116/`, cut from a conversion of the cache Ensembl ships and compared under several flag sets. A larger comparison dataset stands behind the corpora of mouse, the fruit fly, Arabidopsis and yeast alone. That dataset is one chromosome with the organelle genomes for mouse, the fruit fly and Arabidopsis, and every chromosome for yeast. On each corpus and each of those datasets vep-rs matches Ensembl VEP exactly except for the documented divergences ([Intended divergences](intended-divergences.md)).

Every other genome with a ready-made Ensembl VEP cache is compared once, under default flags, on a small generated input. Those genomes are the species of the Ensembl and Ensembl Genomes release directories, the genomes of the bacterial, fungal and protist collections and the human pangenome assemblies. [Supported genomes](supported-genomes.md) lists every cache with its outcome. One pangenome assembly, T2T-CHM13v2.0, has a golden corpus as well. No genome other than human is timed.

`--sift` and `--polyphen` refuse on a cache without their data, and `--af_gnomad` on a cache that declares variation columns without the gnomAD exome column; every other flag that reads data Ensembl ships for human alone runs there as [the flag table](cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone) states.

## Species with a corpus

| Species and assembly | Cache release and `--cache_version` | Genetic codes | Corpus |
| --- | --- | --- | --- |
| `bos_taurus` ARS-UCD2.0 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/ARS-UCD2.0-cattle` |
| `bos_taurus` ARS-UCD2.0, the merged cache | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/ARS-UCD2.0-cattle-merged` |
| `caenorhabditis_elegans` WBcel235 | Ensembl 116; `--cache_version 116`, the default | `MtDNA` table 5; every other region table 1 | `tests/golden/116/WBcel235-worm` |
| `canis_lupus_familiaris` ROS_Cfam_1.0 | Ensembl 116; `--cache_version 116`, the default | every region table 1; the cache holds no mitochondrion | `tests/golden/116/ROS_Cfam_1.0-dog` |
| `capra_hircus` ARS1 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/ARS1-goat` |
| `ciona_intestinalis` KH | Ensembl 116; `--cache_version 116`, the default | `MT` table 13; every other region table 1 | `tests/golden/116/KH-ciona` |
| `danio_rerio` GRCz11 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/GRCz11-zebrafish` |
| `drosophila_melanogaster` BDGP6.54 | Ensembl 116; `--cache_version 116`, the default | `mitochondrion_genome` table 5; every other region table 1 | `tests/golden/116/BDGP6.54-fly` |
| `equus_caballus` EquCab3.0 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/EquCab3.0-horse` |
| `felis_catus` F.catus_Fca126_mat1.0 | Ensembl 116; `--cache_version 116`, the default | every region table 1; the cache holds no organelle sequence | `tests/golden/116/F.catus_Fca126_mat1.0-cat` |
| `gallus_gallus` bGalGal1.mat.broiler.GRCg7b | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/bGalGal1.mat.broiler.GRCg7b-chicken` |
| `macaca_mulatta` Mmul_10 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/Mmul_10-macaque` |
| `mus_musculus` GRCm39 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/GRCm39-mouse` |
| `mus_musculus` GRCm39, the RefSeq cache | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/GRCm39-mouse-refseq` |
| `mus_musculus` GRCm39, the merged cache | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/GRCm39-mouse-merged` |
| `ovis_aries` ARS-UI_Ramb_v3.0 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/ARS-UI_Ramb_v3.0-sheep` |
| `pan_troglodytes` Pan_tro_3.0 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/Pan_tro_3.0-chimpanzee` |
| `rattus_norvegicus` GRCr8 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/GRCr8-rat` |
| `saccharomyces_cerevisiae` R64-1-1 | Ensembl 116; `--cache_version 116`, the default | `Mito` table 3; every other region table 1 | `tests/golden/116/R64-1-1-yeast` |
| `salmo_salar` Ssal_v3.1 | Ensembl 116; `--cache_version 116`, the default | every region table 1; the cache names no mitochondrion | `tests/golden/116/Ssal_v3.1-salmon` |
| `sus_scrofa` Sscrofa11.1 | Ensembl 116; `--cache_version 116`, the default | `MT` table 2; every other region table 1 | `tests/golden/116/Sscrofa11.1-pig` |
| `acanthaster_planci_gca001949145v1` OKI_Apl_1.0 | Ensembl Metazoa 63; `--cache_version 63` | `NC_007788.1` table 9; every scaffold table 1 | `tests/golden/116/OKI_Apl_1.0-acanthaster` |
| `aedes_aegypti_lvpagwg` AaegL5 | Ensembl Metazoa 63; `--cache_version 63` | `MT` table 5; every other region table 1 | `tests/golden/116/AaegL5-aedes` |
| `anopheles_gambiae` AgamP4 | Ensembl Metazoa 63; `--cache_version 63` | `Mt` table 5; every other region table 1 | `tests/golden/116/AgamP4-anopheles` |
| `apis_mellifera` Amel_HAv3.1 | Ensembl Metazoa 63; `--cache_version 63` | `CM009947.2` table 5; every other region table 1 | `tests/golden/116/Amel_HAv3.1-honeybee` |
| `phlebotomus_papatasi` PpapI1 | Ensembl Metazoa 63; `--cache_version 63` | every region table 1 | `tests/golden/116/PpapI1-phlebotomus` |
| `arabidopsis_thaliana` TAIR10 | Ensembl Plants 63; `--cache_version 63` | `Pt` table 11; `Mt` table 1, which the cache states; the nuclear chromosomes table 1 | `tests/golden/116/TAIR10-arabidopsis` |
| `hordeum_vulgare` MorexV3_pseudomolecules_assembly | Ensembl Plants 63; `--cache_version 63` | every region table 1 | `tests/golden/116/MorexV3_pseudomolecules_assembly-barley` |
| `oryza_sativa` IRGSP-1.0 | Ensembl Plants 63; `--cache_version 63` | `Pt` table 11; `Mt` and the twelve chromosomes table 1 | `tests/golden/116/IRGSP-1.0-rice` |
| `sorghum_bicolor` Sorghum_bicolor_NCBIv3 | Ensembl Plants 63; `--cache_version 63` | every region table 1 | `tests/golden/116/Sorghum_bicolor_NCBIv3-sorghum` |
| `triticum_aestivum` IWGSC | Ensembl Plants 63; `--cache_version 63` | every region table 1 | `tests/golden/116/IWGSC-wheat` |
| `zea_mays` Zm-B73-REFERENCE-NAM-5.0 | Ensembl Plants 63; `--cache_version 63` | every region table 1 | `tests/golden/116/Zm-B73-REFERENCE-NAM-5.0-maize` |
| `candida_albicans` GCA000182965v3 | Ensembl Fungi 63; `--cache_version 63` | the eight chromosomes table 12; `Ca22chrM_C_albicans_SC5314` table 4 | `tests/golden/116/GCA000182965v3-candida` |
| `schizosaccharomyces_pombe` ASM294v2 | Ensembl Fungi 63; `--cache_version 63` | `MT` table 4; every other region table 1 | `tests/golden/116/ASM294v2-fission-yeast` |
| `dictyostelium_discoideum` dicty_2.7 | Ensembl Protists 63; `--cache_version 63` | every region table 1, which the cache states | `tests/golden/116/dicty_2.7-dictyostelium` |
| `phytophthora_infestans` ASM14294v1 | Ensembl Protists 63; `--cache_version 63` | every region table 1; the cache holds no organelle | `tests/golden/116/ASM14294v1-phytophthora` |
| `plasmodium_falciparum` GCA000002765v3 | Ensembl Protists 63; `--cache_version 63` | `Pf3D7_API_v3` table 11; `Pf3D7_MIT_v3` table 4; the fourteen chromosomes table 1, which the cache states | `tests/golden/116/GCA000002765v3-plasmodium` |
| `spironucleus_salmonicida` GCA000497125v2 | Ensembl Protists 63; `--cache_version 63` | every region that holds a coding transcript table 6 | `tests/golden/116/GCA000497125v2-spironucleus` |
| `homo_sapiens_gca009914755v4` T2T-CHM13v2.0 | Ensembl rapid release; `--cache_version 107` | every region table 1 | `tests/golden/116/T2T-CHM13v2.0-human-pangenome` |

The genetic codes are the tables the cache names per sequence region. They are read from the core database's `codon_table` attribute by Ensembl's dumper and carried by the conversion. The translation tables themselves are NCBI's genetic codes ([The codon table of each sequence region](cache-setup.md#the-codon-table-of-each-sequence-region)). Genetic code 26 has a page of its own: [A transcript on genetic code 26, which Ensembl VEP reads under the standard code](divergences/pachysolen-nuclear-code-read-as-standard.md). The 35 genomes exercise tables 1, 2, 3, 4, 5, 6, 9, 11, 12 and 13.

Any species in Ensembl release 116 or Ensembl Genomes release 63 with a GFF3 and a peptide FASTA also builds a cache from those files ([Build a cache](cache-setup.md#build-a-cache)). A cache built that way annotates under the codon tables it names, without a concordance measurement. A GFF3 and a FASTA state no genetic code. A build from them names an organelle genome's table with `--codon-table`, or reads it from the release's core dumps.

## Long sequences and CSI-indexed files

A record past position 536,870,912 (2^29), the reach of a `.tbi` index, is annotated from a cache like any other record. 150 ready-made caches hold a sequence longer than that, the longest of 1,805,244,829 bases. On 31,314 records past that position, spread over the 150 caches, vep-rs matches Ensembl VEP exactly except for the documented divergences. One shows there, on five rows of the Tasmanian devil's merged cache: [Two transcripts of a merged cache that carry one `dbID`](divergences/merged-same-dbid-collapse.md).

A cache's own known-variant file and a `--custom` file are read through a `.csi` index, which reaches past that position. Run beside a `--custom` BED file that carries a `.csi` index, the same 31,314 records carry the file's value alike in both engines, on every row both write. On 90 ready-made caches of genomes other than human that carry known variants in such a file, ten records a cache, 89 name the same known variants in both engines. The Anopheles funestus AfunF3 cache differs on nine rows: [Known variants behind a cache variation line with no start](divergences/cache-variation-line-without-a-start.md).

## What each cache carries

The species of this table have no section below. A row gives how the cache names its sequences, what it carries beside its transcripts, and the sequences and the `--hgvs` run of its corpus. Every corpus is also compared under the default flags and `--everything`. A flag whose data a cache lacks runs as in Ensembl VEP: `--everything` leaves SIFT or PolyPhen disabled, and `--check_existing` prints empty columns. A cache whose variation files hold no record carries no known variant.

| Species | Sequence names | Data the cache carries | The corpus |
| --- | --- | --- | --- |
| `caenorhabditis_elegans`, the nematode worm | Chromosomes `I` to `V` and `X`, mitochondrion `MtDNA`. The synonyms give `III` and `MtDNA` their GenBank and RefSeq accessions. `chrIII` resolves through the `chr` rule; `MT`, `M` or `chrM` needs a `--synonyms` line. | No SIFT or PolyPhen matrix, no known variant, no regulatory build. | `III` and `MtDNA`; `--hgvs` |
| `canis_lupus_familiaris`, the dog | Chromosomes `1` to `38`, `X` and `Y`, scaffolds by GenBank accession, no mitochondrion. Synonyms: RefSeq and GenBank accessions, a scaffold's assembly name. | SIFT; known variants without a frequency column. No PolyPhen matrix or regulatory build. | Chromosome 38 and one scaffold; `--hgvs --sift b` |
| `capra_hircus`, the goat | Chromosomes `1` to `29` and `MT`, scaffolds by accession. Synonyms: three accessions a chromosome, a scaffold's RefSeq accession and assembly name, none for `MT`. | SIFT; known variants by dbSNP identifier. No PolyPhen matrix, regulatory build or allele frequency. | Chromosome 25, `MT` and one scaffold; `--hgvs --sift b` |
| `ciona_intestinalis`, the sea squirt | Chromosomes `1` to `14`, `MT`, scaffolds by accession. Synonyms: a sequence's other accessions. | UniProt accessions and protein domains. The variation files hold no record; no SIFT or PolyPhen matrix, no regulatory feature. | `MT`, chromosome 6 and one scaffold; `--hgvs`. The mitochondrial records reach each of the codons table 13 reads unlike the standard code. |
| `danio_rerio`, the zebrafish | Chromosomes `1` to `25`, `MT`, scaffolds by accession, alternate sequences (`ALT_CTG9_3_1`). Synonyms: GenBank and RefSeq accessions, none for `MT`. | SIFT; known variants. No PolyPhen matrix, regulatory build or allele frequency. | Chromosome 25, `MT` and one scaffold; `--hgvs --sift b` |
| `equus_caballus`, the horse | Chromosomes `1` to `31`, `X`, `MT`, scaffolds by accession. Synonyms: `chr` names, RefSeq and GenBank accessions. | SIFT; known variants with `rs` identifiers. No PolyPhen matrix or regulatory build. | Chromosome 31, `MT` and one scaffold; `--hgvs --sift b` |
| `felis_catus`, the cat | Chromosomes `A1` to `F2` and `X`, scaffolds by accession, no mitochondrion. Synonyms: each chromosome's GenBank accession, a scaffold's assembly name. | SIFT; known variants with dbSNP identifiers. No PolyPhen matrix or regulatory build. | Chromosome `F1` and one scaffold; `--hgvs --sift b` |
| `gallus_gallus`, the chicken | Chromosomes `1` to `39`, `W`, `Z`, `MT`, scaffolds by accession. Synonyms: RefSeq and GenBank accessions, a scaffold's assembly name. `--assembly` takes the dotted name whole. | SIFT; known variants; a regulatory build. No PolyPhen matrix. | Chromosome 25, `MT` and one scaffold; `--hgvs --sift b` |
| `macaca_mulatta`, the rhesus macaque | Chromosomes `1` to `20`, `X`, `Y`, `MT`, scaffolds by accession. Synonyms: GenBank and RefSeq accessions, a scaffold's `chrUn` name, none for `MT`. | The variation files hold no record; no SIFT or PolyPhen matrix, no regulatory build. | Chromosome 19, `MT` and one scaffold; `--hgvs` |
| `ovis_aries`, the sheep | Chromosomes `1` to `26`, `X`, `Y`, `MT`, scaffolds by accession. Synonyms: GenBank and RefSeq accessions, a scaffold's `scaffold_` name. | SIFT; known variants. No PolyPhen matrix or regulatory build. | Chromosome 24, `MT` and one scaffold; `--hgvs --sift b` |
| `pan_troglodytes`, the chimpanzee | Chromosomes `1`, `2A`, `2B`, `3` to `22`, `X`, `Y`, `MT`, scaffolds by accession. Synonyms: RefSeq and GenBank accessions, a scaffold's contig names. | The variation files hold no record; no SIFT or PolyPhen matrix, no regulatory build. | Chromosome 21, `MT` and one scaffold; `--hgvs` |
| `rattus_norvegicus`, the rat | Chromosomes `1` to `20`, `X`, `Y`, `MT`, scaffolds by accession. Synonyms: RefSeq and GenBank accessions, a scaffold's assembly name. | SIFT; known variants. No PolyPhen matrix or regulatory build. | Chromosome 12, `MT` and one scaffold; `--hgvs --sift b` |
| `salmo_salar`, the Atlantic salmon | Chromosomes `1` to `29`, scaffolds by accession, no mitochondrion. Synonyms: the `ssa` names and GenBank accessions of the chromosomes, the `contig` names of the scaffolds. | A regulatory build. The variation files hold no record; no SIFT or PolyPhen matrix. | Chromosome 8 and one scaffold; `--hgvs` |
| `sus_scrofa`, the pig | Chromosomes `1` to `18`, `X`, `Y`, `MT`, scaffolds by accession. Synonyms: `chr` names, GenBank and RefSeq accessions, a scaffold's UCSC and contig names, none for `MT`. | SIFT; known variants; a regulatory build. No PolyPhen matrix. | Chromosome 18, `MT` and one scaffold; `--hgvs --sift b` |
| `acanthaster_planci_gca001949145v1`, the crown-of-thorns starfish | Scaffolds by accession (`DF978489.1`) and the mitochondrion `NC_007788.1`, which the cache states circular. Synonyms: four names a scaffold; `MT` and `NC_007788` for the mitochondrion. | UniProt accessions and protein domains. No known variant, no SIFT or PolyPhen matrix, no regulatory build. | Two scaffolds and the mitochondrion; `--hgvs`. Records reach each of the codons table 9 reads unlike table 1. |
| `aedes_aegypti_lvpagwg`, the yellow fever mosquito | Chromosomes `1`, `2`, `3`, `MT`, scaffolds by accession without its version (`NIGP01000988`). Synonyms: `Mt` and the mitochondrion's accessions; a scaffold's versioned accession, RefSeq accession and assembly name. | SIFT; known variants. No PolyPhen matrix or regulatory build. | `MT` and 24 scaffolds; `--hgvs --sift b` |
| `anopheles_gambiae`, the malaria mosquito | Chromosome arms `2L`, `2R`, `3L`, `3R` and `X`, mitochondrion `Mt`, `UNKN` and `Y_unplaced`. Synonyms: INSDC accessions. | SIFT; known variants on all eight sequences. No PolyPhen matrix or regulatory build. | `X`, `Mt` and `Y_unplaced`; `--hgvs --sift b` |
| `apis_mellifera`, the honey bee | Every sequence by GenBank accession: chromosomes `CM009931.2` to `CM009946.2`, mitochondrion `CM009947.2`, scaffolds. Synonyms: linkage group names (`LG16`, `Group16`), `MT`, RefSeq accessions. | UniProt accessions and protein domains. No known variant, no SIFT or PolyPhen matrix, no regulatory data. | `CM009946.2`, the mitochondrion and one scaffold; `--hgvs` |
| `phlebotomus_papatasi`, the sand fly | Scaffolds alone, by GenBank accession with its version (`JH660948.1`); no chromosome or organelle. Synonyms: the assembly's scaffold names (`Scaffold1`). | SIFT; known variants. No PolyPhen matrix or regulatory build. | Seven scaffolds; `--hgvs --sift b` |
| `hordeum_vulgare`, barley | Chromosomes `1H` to `7H`, scaffolds by accession, no organelle. Synonyms: a chromosome's INSDC accession. `--assembly` takes `MorexV3_pseudomolecules_assembly` whole. | SIFT; known variants on the seven chromosomes and four scaffolds. No PolyPhen matrix or regulatory build. | Chromosome `1H`; `--hgvs --sift b` |
| `oryza_sativa`, rice | Chromosomes `1` to `12`, `Mt`, `Pt`, no scaffold. Synonyms: `Chr1` to `Chr12`, GenBank and RefSeq accessions. | SIFT; known variants on the twelve chromosomes, none on `Mt` or `Pt`. No PolyPhen matrix or regulatory build. | Chromosome 9, `Mt` and `Pt`; `--hgvs --sift b` |
| `sorghum_bicolor`, sorghum | Chromosomes `1` to `10`, scaffolds as `super_` and a number, no organelle. The synonyms link chromosome `9` to `chr9` and two accessions. | SIFT; known variants. No PolyPhen matrix or regulatory build. | Chromosome 9 and one scaffold; `--hgvs --sift b` |
| `triticum_aestivum`, bread wheat | Chromosomes `1A` to `7D` and `Un`, no organelle. Synonyms: `chr` names and INSDC accessions. | SIFT; known variants, indexed with `.csi`: eighteen of the 22 sequences are longer than 536,870,912 bases, the reach of a `.tbi` index. No PolyPhen matrix or regulatory build. | Chromosome `5D`, with 39 of its 272 records past position 536,870,912; `--hgvs --sift b` |
| `zea_mays`, maize | Chromosomes `1` to `10`, scaffolds (`scaf_100`), no organelle. The synonyms file is empty; a `chr` prefix resolves. `--assembly` takes `Zm-B73-REFERENCE-NAM-5.0` whole. | SIFT; known variants on the ten chromosomes and two scaffolds. No PolyPhen matrix or regulatory build. | Chromosome `10` and one scaffold; `--hgvs --sift b` |
| `schizosaccharomyces_pombe`, fission yeast | Chromosomes `I`, `II`, `III`, the sequences `MTR` and `AB325691`, mitochondrion `MT`. Synonyms: `chr3` and an accession for `III`, accessions for `MT` and `AB325691`. | Known variants on all six sequences, with their citations. No SIFT or PolyPhen matrix, no regulatory build. | `III`, `MT` and `AB325691`; `--hgvs`. Records change the codon table 4 reads as tryptophan and table 1 as a stop. |
| `dictyostelium_discoideum`, the slime mould | Chromosomes `1` to `6`, unplaced sequences by accession (`CH709148`), no organelle. Synonyms: versioned accessions. | UniProt accessions and protein domains. No known variant, no SIFT or PolyPhen matrix, no regulatory feature. | Chromosome 6 and two unplaced sequences; `--hgvs` |
| `phytophthora_infestans`, the potato blight pathogen | Supercontigs alone (`supercont1.1`), no chromosome or organelle. The synonyms give a sequence a second name, the INSDC accession for the three of the corpus. | Known variants. No SIFT or PolyPhen matrix, no regulatory build, no frequency column. | Three supercontigs; `--hgvs` |

The sections below give a worked example for each set of ready-made caches, and the species whose cache needs more than a row:

- the Ensembl release directory: [Mus musculus](#mus-musculus), with [Bos taurus](#bos-taurus) for a merged cache;
- a plant, an insect and a fungus: [Arabidopsis thaliana](#arabidopsis-thaliana) from Ensembl Plants and, from the Ensembl release directory, [Drosophila melanogaster](#drosophila-melanogaster) and [Saccharomyces cerevisiae](#saccharomyces-cerevisiae);
- Ensembl Protists: [Plasmodium falciparum](#plasmodium-falciparum), with the installer's command beside the conversion;
- a nuclear genetic code other than the standard one: [Candida albicans](#candida-albicans) and [Spironucleus salmonicida](#spironucleus-salmonicida);
- the collections and the human pangenome: [a bacterial genome](#a-genome-of-a-bacterial-collection), [a fungal genome](#a-genome-of-a-fungal-collection) and [a pangenome assembly](#a-human-pangenome-assembly).

A cache larger than those of the worked examples converts and runs by the same two commands under its own species directory, assembly and `--cache_version`.

## Mus musculus

The mouse at release 116, whose cache carries SIFT matrices, a regulatory build and the mitochondrion as `MT` under genetic code table 2:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/mus_musculus/116_GRCm39" /path/to/mouse_cache

# Annotate with SIFT predictions and known variants
vep -i mouse.vcf -o mouse.txt --offline --json_cache /path/to/mouse_cache \
  --species mus_musculus --assembly GRCm39 \
  --fasta Mus_musculus.GRCm39.dna.primary_assembly.fa \
  --sift b --check_existing
```

The cache names the chromosomes `1` to `19`, `X`, `Y` and `MT` and 28 unplaced scaffolds by accession (`GL456210.1`, `JH584296.1`). Its synonyms resolve the `chr`-prefixed names, `chrM` and the RefSeq and GenBank accessions (`NC_000085.7` and `CM001012.3` for chromosome 19, `NC_005089.1` for `MT`). Under `--dir_cache` the cache sits at `mus_musculus/116_GRCm39`, which the default `--cache_version` selects.

The mouse cache has SIFT matrices and no PolyPhen matrix: `--everything` runs with PolyPhen disabled, as in Ensembl VEP, and SIFT populated.

Among the mouse, fruit fly, Arabidopsis and yeast caches, some data is in the mouse cache alone. The flags that read it print values here and empty cells on the other three ([the flag table](cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone)):

| Data | Which of those caches carry it |
| --- | --- |
| The TSL, CCDS and GENCODE attributes, gene phenotypes and variant citations | the mouse cache alone |
| APPRIS | the mouse and fruit fly caches |
| UniProt accessions and protein domains | each of them |
| Variant synonyms | the mouse and Arabidopsis caches |

The corpus `tests/golden/116/GRCm39-mouse` compares both engines on chromosome 19, the mitochondrion and the unplaced scaffolds under the default flags, `--everything` and `--hgvs --sift b`.

Ensembl also ships a RefSeq and a merged cache for the mouse at release 116, `mus_musculus_refseq_vep_116_GRCm39.tar.gz` and `mus_musculus_merged_vep_116_GRCm39.tar.gz`. Each converts with the same converter, given the directory inside Ensembl VEP's layout, into `mus_musculus_refseq/116_GRCm39` or `mus_musculus_merged/116_GRCm39` under the cache root, where `--species mus_musculus` with `--refseq` or `--merged` finds it. The corpora `tests/golden/116/GRCm39-mouse-refseq` and `tests/golden/116/GRCm39-mouse-merged` compare both engines on chromosome 19, the mitochondrion and unplaced scaffolds under `--refseq` and `--merged`, each with and without `--hgvs`.

## Bos taurus

Cattle at release 116, whose cache carries SIFT matrices, known variants, a regulatory build and the mitochondrion as `MT` under genetic code table 2. Ensembl also ships a merged cache for cattle, `bos_taurus_merged_vep_116_ARS-UCD2.0.tar.gz`: the RefSeq annotation's transcripts beside Ensembl's own, with known variants, SIFT matrices and a regulatory build.

The cache names the chromosomes `1` to `29`, `X`, `Y` and `MT` and its unplaced scaffolds by GenBank accession (`NKLS02000103.1`). Its synonyms resolve the `chr`-prefixed names and the RefSeq and GenBank accessions (`NC_037352.1` and `CM008192.2` for chromosome 25). Under `--dir_cache` the cache sits at `bos_taurus/116_ARS-UCD2.0`, which the default `--cache_version` selects. The merged cache converts with the same converter, given the directory inside Ensembl VEP's layout, so that `info.json` names `source_type` `merged`. Under `--dir_cache` it sits at `bos_taurus_merged/116_ARS-UCD2.0`, which `--species bos_taurus` with `--merged` selects.

The cattle cache has SIFT matrices and no PolyPhen matrix: `--everything` runs with PolyPhen disabled, as in Ensembl VEP, and SIFT populated. Its known variants carry no frequency, so `--check_existing` fills `Existing_variation` and the frequency columns stay empty. A run on the merged cache prints `SOURCE`, `Ensembl` or `RefSeq`, on its transcript rows. That cache's `info.txt` names no alignment file, so no model is BAM-edited and a run prints no `GIVEN_REF`, `USED_REF`, `BAM_EDIT` or `REFSEQ_MATCH` column.

The corpus `tests/golden/116/ARS-UCD2.0-cattle` compares both engines on chromosome 25, the mitochondrion and one unplaced scaffold under the default flags, `--everything` and `--hgvs --sift b`. The corpus `tests/golden/116/ARS-UCD2.0-cattle-merged` compares them on chromosome 25 and one unplaced scaffold under `--merged`, `--merged --everything` and `--merged --hgvs --sift b`. It sits on no organelle.

## Drosophila melanogaster

The fruit fly at release 116, whose cache names the mitochondrion `mitochondrion_genome` under genetic code table 5 while its variant sets say `MT`:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/drosophila_melanogaster/116_BDGP6.54" /path/to/fly_cache

# Annotate; the synonyms file links the input's MT to the cache's name
printf 'mitochondrion_genome\tMT\n' > synonyms.txt
vep -i fly.vcf -o fly.txt --offline --json_cache /path/to/fly_cache \
  --species drosophila_melanogaster --assembly BDGP6.54 \
  --synonyms synonyms.txt
```

The cache names the chromosome arms `2L`, `2R`, `3L` and `3R`, the chromosomes `4`, `X` and `Y`, `mitochondrion_genome`, `rDNA` and 20 scaffolds (`211000022278760`, `Unmapped_Scaffold_8_D1580_D1567`). Its synonyms resolve the `chr`-prefixed names, `chrM` and the GenBank and RefSeq accessions (`AE014134.6` and `NT_033779.5` for `2L`, `KJ947872.2` and `NC_024511.2` for the mitochondrion). `MT` is the name Ensembl's variant files give the fly's mitochondrion, and it is neither a chromosome nor a synonym of the cache. So an input spelling it that way is skipped in both engines unless a `--synonyms` file links the two names. The flag replaces the cache's own table, so the file also lists any other synonym the input needs.

The fly cache has no SIFT or PolyPhen matrices and no variant. So `--everything` runs with both disabled, as in Ensembl VEP, and `--check_existing` prints empty columns. The corpus `tests/golden/116/BDGP6.54-fly` compares both engines on 2L, the mitochondrion, 4, Y, rDNA and seven of the scaffolds under the default flags, `--everything` and `--hgvs`.

## Saccharomyces cerevisiae

Yeast at release 116, whose cache names the chromosomes `I` to `XVI` and `Mito`, the mitochondrion under genetic code table 3, and carries SIFT matrices and the SGRP variants:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/saccharomyces_cerevisiae/116_R64-1-1" /path/to/yeast_cache

# Annotate; the SGRP identifiers come back through --check_existing
vep -i yeast.vcf -o yeast.txt --offline --json_cache /path/to/yeast_cache \
  --species saccharomyces_cerevisiae --assembly R64-1-1 \
  --check_existing --sift b
```

The cache's synonyms are the GenBank accessions alone (`BK006935.2` for `I`, `AJ011856.1` for `Mito`). A `chrI` spelling is reached through the `chr` rule. A mitochondrion spelled `MT`, `M` or `chrM` is skipped unless a `--synonyms` file links it to `Mito`. Under `--dir_cache` the cache sits at `saccharomyces_cerevisiae/116_R64-1-1`, which the default `--cache_version` selects.

Under table 3 a `CTN` codon reads as threonine, `TGA` as tryptophan and `ATA` as a start codon. So a mitochondrial change the nuclear code calls missense can be synonymous, and a stop codon changed to `TGA` is a stop lost. The cache has no PolyPhen matrix: `--everything` runs with PolyPhen disabled, as in Ensembl VEP, and the SIFT column populated. The corpus `tests/golden/116/R64-1-1-yeast` compares both engines on the mitochondrion and twelve nuclear chromosomes under the default flags, `--everything` and `--hgvs`.

## Arabidopsis thaliana

Arabidopsis thaliana at Ensembl Plants release 63, whose plastid `Pt` translates under genetic code table 11 and whose mitochondrion `Mt` under table 1:

```bash
# Convert the Ensembl VEP cache, extracted under $VEP_CACHE
# The converter is built once
cargo build --release --locked \
  --manifest-path tools/vep-storable-to-json/Cargo.toml
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/arabidopsis_thaliana/63_TAIR10" /path/to/arabidopsis_cache

# Annotate; the cache version is the Plants release
vep -i arabidopsis.vcf -o arabidopsis.txt --offline \
  --json_cache /path/to/arabidopsis_cache \
  --species arabidopsis_thaliana --assembly TAIR10 --cache_version 63
```

The cache names the chromosomes `1` to `5`, `Mt` and `Pt`. Its synonyms resolve the `Chr1` to `Chr5` spellings and the GenBank and RefSeq accessions, so an input spelled any of those ways annotates. Among the accessions are `CP002687.1` and `NC_003075.7` for chromosome 4, `Y08501.2` and `NC_001284.2` for `Mt`, and `AP000423.1` and `NC_000932.1` for `Pt`. Under `--dir_cache` the cache sits at `arabidopsis_thaliana/63_TAIR10`, which `--cache_version 63` selects. `--json_cache` names the directory itself and needs no version.

The Arabidopsis cache carries SIFT matrices, no PolyPhen matrices and variants on the five nuclear chromosomes alone. So `--everything` runs with PolyPhen disabled, as in Ensembl VEP, and `--check_existing` prints empty columns on `Mt` and `Pt`. The corpus `tests/golden/116/TAIR10-arabidopsis` compares both engines on chromosome 4, `Mt` and `Pt` under the default flags, `--everything` and `--hgvs`. The plastid's start codon shapes are among its records.

[Build a cache](cache-setup.md#build-a-cache) says what the installer fetches and what each input supplies. The installer builds a cache of the same genome from the release's own files and reads both organelle tables from the release's core dumps:

```bash
# Install the cache in one command
vep-cache-builder --auto arabidopsis_thaliana TAIR10 \
  --release 63 --dir_cache /path/to/vep_dir

# Annotate from the installed cache
vep -i arabidopsis.vcf -o arabidopsis_installed.txt --offline \
  --dir_cache /path/to/vep_dir --species arabidopsis_thaliana \
  --cache_version 63 --assembly TAIR10
```

A `dna_index/` file that its published `.fai` and `.gzi` do not open is indexed by the builder, with a warning. At release 63 the installation prints a warning that opens `The .fai and .gzi published beside Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz do not open it, so both are written here from the file itself`. The installer writes both indexes from the genome FASTA it fetched, and the installation ends with exit status 0. Another line counts the four one-exon transcripts of the next paragraph, each built as the release's `translation` table places its coding sequence, and no protein FASTA record is set aside.

The installed cache holds the known variants, with synonyms on 181,418 of them as in Ensembl VEP's own cache, and the SIFT matrices the release publishes for the genome. The installation's last lines name what the cache lacks: PolyPhen predictions, regulatory features and known variants on `Mt` and `Pt`, on which it read no variation record, in the variation VCF or in the `variation_feature` table.

The release's GFF3 gives four one-exon mitochondrial transcripts a `CDS` line over the whole exon: `ATMG00410.1`, `ATMG00730.1`, `ATMG01080.1` and `ATMG01170.1`. The release's `translation` table starts each coding sequence after the bases the file's own `five_prime_UTR` line covers, as Ensembl VEP's cache does. A build that reads the core dumps, the installer's among them, starts each coding sequence where the table does.

A build without the core dumps takes each `CDS` line as written and sets the four protein FASTA records aside: each record's length is neither its transcript's whole codons nor one fewer. With `--genome-fasta` that build translates the four peptides from the genome; without it the four transcripts carry no peptide.

## Plasmodium falciparum

The malaria parasite Plasmodium falciparum 3D7 at Ensembl Protists release 63, whose apicoplast `Pf3D7_API_v3` translates under genetic code table 11 and whose mitochondrion `Pf3D7_MIT_v3` under table 4:

```bash
# Download and extract the cache Ensembl ships
site=https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63
archive=plasmodium_falciparum_vep_63_GCA000002765v3.tar.gz
curl -fO "$site/variation/indexed_vep_cache/$archive"
tar -xzf "$archive" -C "$VEP_CACHE"

# Convert it
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/plasmodium_falciparum/63_GCA000002765v3" \
  /path/to/plasmodium_cache

# Annotate; the cache version is the Protists release
vep -i plasmodium.vcf -o plasmodium.txt --offline \
  --json_cache /path/to/plasmodium_cache \
  --species plasmodium_falciparum --assembly GCA000002765v3 \
  --cache_version 63
```

The cache names the fourteen chromosomes `Pf3D7_01_v3` to `Pf3D7_14_v3`, the apicoplast `Pf3D7_API_v3` and the mitochondrion `Pf3D7_MIT_v3`. Its synonyms give each chromosome two shorter names, `MAL5` and `Pf3D7_05` for `Pf3D7_05_v3`, and the two organelles none. Under `--dir_cache` the cache sits at `plasmodium_falciparum/63_GCA000002765v3`, which `--cache_version 63` selects.

The Plasmodium cache carries no known variant, no SIFT or PolyPhen matrix and no regulatory feature. So `--everything` runs with SIFT and PolyPhen disabled, as in Ensembl VEP, and prints empty known-variant and frequency columns.

The installer builds a cache of the same genome from the release's own files and reads both organelle tables from the release's core dumps:

```bash
# Install the cache in one command
vep-cache-builder --auto plasmodium_falciparum GCA000002765v3 \
  --release 63 --dir_cache /path/to/vep_dir

# Annotate from the installed cache
vep -i plasmodium.vcf -o plasmodium_installed.txt --offline \
  --dir_cache /path/to/vep_dir --species plasmodium_falciparum \
  --cache_version 63 --assembly GCA000002765v3
```

The corpus `tests/golden/116/GCA000002765v3-plasmodium` compares both engines on chromosome 5, the apicoplast and the mitochondrion under the default flags, `--everything` and `--hgvs`. Nine of its mitochondrial records make the codon that table 4 reads as tryptophan and the standard code as a stop.

## Candida albicans

The yeast Candida albicans SC5314 at Ensembl Fungi release 63, whose chromosomes translate under genetic code table 12, the alternative yeast nuclear code, and whose mitochondrion translates under table 4.

The cache names the eight chromosomes `Ca22chr1A_C_albicans_SC5314` to `Ca22chr7A_C_albicans_SC5314` and `Ca22chrRA_C_albicans_SC5314`, and the mitochondrion `Ca22chrM_C_albicans_SC5314`. Its synonyms resolve a chromosome's shorter name, `CalbSC5314_Chr7` for `Ca22chr7A_C_albicans_SC5314`. Under `--dir_cache` the cache sits at `candida_albicans/63_GCA000182965v3`, which `--cache_version 63` selects.

Every one of the eight chromosomes names table 12, which reads one codon as serine where the standard code reads leucine. The cache carries gene symbols, UniProt accessions and protein domains. It carries no known variant, no SIFT or PolyPhen matrix and no regulatory feature, so `--everything` runs with SIFT and PolyPhen disabled, as in Ensembl VEP.

The corpus `tests/golden/116/GCA000182965v3-candida` compares both engines on chromosome 7 and the mitochondrion under the default flags, `--everything` and `--hgvs`. Twenty-four records of chromosome 7, on both strands, make that codon, keep its residue, alter it or turn the closing stop codon into it.

## Spironucleus salmonicida

The diplomonad Spironucleus salmonicida at Ensembl Protists release 63, whose cache translates every coding transcript under genetic code table 6, the ciliate, dasycladacean and hexamita nuclear code.

The cache names its chromosomes and unplaced sequences by accession (`CM035318`, `CM035324`, `AUWU02000032`). Its synonyms resolve a chromosome's number with and without the `chr` prefix (`3` and `chr3` for `CM035318`), versioned accessions (`CM035324.1`) and the assembly's contig names (`tig00000055` for `AUWU02000032`). Under `--dir_cache` the cache sits at `spironucleus_salmonicida/63_GCA000497125v2`, which `--cache_version 63` selects.

Table 6 has one stop codon and translates as glutamine the two other codons that the standard code reads as stops. A substitution that makes one of those two inside a coding sequence is a synonymous or a missense change on this genome, where the standard code would name a stop gained. The closing stop codon turned into one of them is a stop lost.

The cache carries UniProt accessions and protein domains. It carries no known variant, no SIFT or PolyPhen matrix and no regulatory feature, it names no organelle, and none of its transcripts has an untranslated region. The corpus `tests/golden/116/GCA000497125v2-spironucleus` compares both engines on the chromosomes `CM035318` and `CM035324` and the unplaced sequence `AUWU02000032` under the default flags, `--everything` and `--hgvs`.

Under `--everything` and `--hgvs`, 20 frameshift rows of the corpus carry the documented HGVSp difference of [a transcript under another genetic code](divergences/organelle-peptide-translation-table.md). Ensembl VEP translates the altered sequence with the standard code, which reads the two glutamine codons as stops.

## A genome of a bacterial collection

Every bacterial and archaeal genome ships inside a collection, a folder of single-genome archives on the Ensembl Bacteria site. One archive converts to a cache like any other. Escherichia coli K-12 MG1655, of `bacteria_0_collection` at release 63:

```bash
# Download and extract the cache Ensembl ships
genome=escherichia_coli_str_k_12_substr_mg1655_gca_000005845
site=https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-63
dir=$site/variation/indexed_vep_cache/bacteria_0_collection
curl -fO "$dir/${genome}_vep_63_ASM584v2.tar.gz"
tar -xzf "${genome}_vep_63_ASM584v2.tar.gz" -C "$VEP_CACHE"

# Convert it
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/$genome/63_ASM584v2" /path/to/ecoli_cache

# Annotate; the cache version is the Ensembl Bacteria release
vep -i ecoli.vcf -o ecoli.txt --offline --json_cache /path/to/ecoli_cache \
  --species "$genome" --assembly ASM584v2 --cache_version 63
```

The species directory is the genome's whole name, with its strain and assembly accession, and `--species` takes it whole. The cache names one sequence, `Chromosome`, states it circular and translates every transcript under genetic code table 11. It carries no known variant, no SIFT or PolyPhen matrix and no regulatory build.

The synonyms file of such a cache is the collection's own: it also pairs the sequences of the collection's other genomes with their accessions. For this genome the INSDC accession `U00096.3` resolves to `Chromosome`, and the RefSeq accession `NC_000913.3` does not. Neither `--is_multispecies` nor `--genomes` is needed, and both are refused by name: they describe an Ensembl database, and vep-rs reads a local cache alone. The installer refuses a collection genome before any download and names the conversion.

## A genome of a fungal collection

A fungal or protist genome without a core database of its own ships the same way, in a collection folder of the Ensembl Fungi or Ensembl Protists site ([Genome checks](genome-checks.md#where-the-cache-comes-from) gives each directory). Candida auris, of `fungi_ascomycota4_collection` at release 63:

```bash
# Download and extract the cache Ensembl ships
genome=candida_auris_gca_003013715
site=https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-63
dir=$site/variation/indexed_vep_cache/fungi_ascomycota4_collection
curl -fO "$dir/${genome}_vep_63_ASM301371v2.tar.gz"
tar -xzf "${genome}_vep_63_ASM301371v2.tar.gz" -C "$VEP_CACHE"

# Convert it
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/$genome/63_ASM301371v2" /path/to/candida_cache

# Annotate; the cache version is the Ensembl Fungi release
vep -i candida.vcf -o candida.txt --offline \
  --json_cache /path/to/candida_cache \
  --species "$genome" --assembly ASM301371v2 --cache_version 63
```

The cache names seven sequences, `1` to `7`, none of them a mitochondrion, and translates every transcript under genetic code table 12. It carries no known variant, no SIFT or PolyPhen matrix and no regulatory build, and its synonyms file is the collection's.

## A human pangenome assembly

Ensembl ships a cache for each assembly of the human pangenome in its rapid-release directory. The species directory carries the assembly's accession, and the cache version is 107 rather than the Ensembl release. T2T-CHM13v2.0:

```bash
# Download and extract the cache Ensembl ships
site=https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens
dir=$site/GCA_009914755.4/ensembl/variation/2022_10/indexed_vep_cache
archive=Homo_sapiens-GCA_009914755.4-2022_10.tar.gz
curl -fO "$dir/$archive"
tar -xzf "$archive" -C "$VEP_CACHE"

# Convert it
tools/vep-storable-to-json/target/release/vep-storable-to-json \
  "$VEP_CACHE/homo_sapiens_gca009914755v4/107_T2T-CHM13v2.0" \
  /path/to/chm13_cache

# Annotate; the cache version is the one the cache's directory names
vep -i chm13.vcf -o chm13.txt --offline --json_cache /path/to/chm13_cache \
  --species homo_sapiens_gca009914755v4 --assembly T2T-CHM13v2.0 \
  --cache_version 107
```

The cache names the sequences `1` to `22`, `X` and `Y`, with no `chr` prefix and no mitochondrion. Its synonyms resolve each sequence's RefSeq and GenBank accessions (`NC_060925.1` and `CP068277.2` for `1`), and a `chr1` spelling is reached through the `chr` rule. Every transcript reads genetic code table 1. Both engines read the species directory as a species of its own and not as `homo_sapiens`. The directory has no RefSeq or merged cache, so `--refseq` and `--merged` are refused by name.

The cache carries no SIFT or PolyPhen matrix and no regulatory build, and its variation files hold no record, so `--check_existing` names no known variant. [Supported genomes](supported-genomes.md#human-pangenome-assemblies) lists every assembly by its species directory. The installer refuses a pangenome assembly before any download and names the conversion.

The corpus `tests/golden/116/T2T-CHM13v2.0-human-pangenome` compares both engines on chromosome 21 under the default flags, `--everything` and `--hgvs`.

## Where next

- [Cache setup](cache-setup.md): building a cache, naming it, the reference FASTA and the chromosome names.
- [Cache conversion](cache-conversion.md): what a converted Ensembl VEP cache carries.
- [The flag table](cli/cache-species-runtime.md#flags-that-read-data-ensembl-ships-for-human-alone): what each flag that reads data Ensembl ships for human alone does on another species.
- [Intended divergences](intended-divergences.md): where the two engines differ, and why.
