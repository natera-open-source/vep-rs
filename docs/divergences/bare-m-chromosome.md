# A bare `M` under a cache whose synonyms name `chrM`

A difference in which rows a record gets, on an input chromosome no measured dataset uses; nothing is set aside, and a record on `M` is charged to vep-rs in both figures.

Ensembl VEP validates a record's chromosome against the cache's chromosomes and the synonyms of each (`Parser.pm` lines 770-800, `_have_chr`): the name itself, then the name with `chr` removed or added. When every one of those misses and `MT` is a cache chromosome or a synonym of one, it renames `M` to `MT` on the record (lines 794-796). An Ensembl cache's `chr_synonyms.txt` lists `chrM` as a synonym of `MT`, so under it an input `M` is accepted as `chrM` and never renamed.

The annotation sources then resolve the name on their own (`BaseVEP.pm` lines 689-738, `get_source_chr_name`), through the cache's chromosomes and the synonyms of the input name alone, with no rename, so `M` reaches no cache chromosome. The record passes validation, overlaps nothing and is written as `intergenic_variant` with Location `M:<pos>`. Under a cache without that synonym the record is renamed and annotated as `MT`, with Location `MT:<pos>`.

vep-rs prints the name Ensembl VEP prints in both cases: `M` under the Ensembl cache, `MT` without the synonym. For its lookups it resolves a bare `M` the cache does not know to the cache's `MT`, under whatever name the cache holds it. So the record is annotated against the mitochondrial transcripts in both cases. `chrM` and `MT` are read alike by both engines.
