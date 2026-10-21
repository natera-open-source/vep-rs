# A bare `M` under a cache whose synonyms name `chrM`

A difference in which rows a record gets, on an input chromosome no measured dataset uses; nothing is set aside for a record on `M`, which is charged to vep-rs in both figures. The mate rows of a breakend mate written `M`, or `chrM` under a cache that links it to `MT`, are set aside under [class 24's mask](defect-p.md#247-how-the-comparator-treats-this-class).

Ensembl VEP validates a record's chromosome against the cache's chromosomes and the synonyms of each (`Parser.pm` lines 770-800, `_have_chr`): the name itself, then the name with `chr` removed or added. When every one of those misses and `MT` is a cache chromosome or a synonym of one, it renames `M` to `MT` on the record (lines 794-796). An Ensembl cache's `chr_synonyms.txt` lists `chrM` as a synonym of `MT`, so under it an input `M` is accepted as `chrM` and never renamed.

The annotation sources then resolve the name on their own (`BaseVEP.pm` lines 689-738, `get_source_chr_name`), through the cache's chromosomes and the synonyms of the input name alone, with no rename, so `M` reaches no cache chromosome. The record passes validation, overlaps nothing and is written as `intergenic_variant` with Location `M:<pos>`. Under a cache without that synonym the record is renamed and annotated as `MT`, with Location `MT:<pos>`.

vep-rs prints the name Ensembl VEP prints in both cases: `M` under the Ensembl cache, `MT` without the synonym. For its lookups it resolves a bare `M` the cache does not know to the cache's `MT`, under whatever name the cache holds it. So the record is annotated against the mitochondrial transcripts in both cases. `chrM` and `MT` are read alike by both engines.

vep-rs resolves a breakend mate's chromosome from its written name against the cache's chromosomes and their synonyms, with no rename of its own; a bare `M` mate reaches the cache's `MT` the same way the record does. Its mate rows are therefore annotated against the mitochondrial transcripts. Ensembl VEP writes none for it: its region fetch resolves the mate through `get_source_chr_name` alone (`BaseVEP.pm` lines 689-738), where `M` reaches no cache chromosome, and a transcript is kept only for a mate whose name, `chr` stripped, spells the transcript's chromosome (`AnnotationType/Transcript.pm` lines 133-134).

A `chrM` mate under a cache whose synonyms do not link it to `MT` stays unresolved and gets no mate rows in either engine. Under a cache that links them vep-rs writes the mate rows and Ensembl VEP still none: `chrM`, `chr` stripped, is `M`, no spelling of `MT`.
