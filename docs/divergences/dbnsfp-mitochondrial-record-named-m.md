# dbNSFP on GRCh37: a mitochondrial record named `M`

**The defect in one sentence.** Ensembl VEP's dbNSFP plugin names every record on `MT` or `chrMT` `M` before it reads its table and matches a row on position and alternate allele alone. On GRCh37, where `M` names another sequence than `MT`, it so writes the value of whatever row a table named `M` holds at the record's position; vep-rs reads such a table for such a record only where its rows state the assembly's bases.

**Why this is a defect.** On GRCh37 `MT` and `chrMT` name the assembly's mitochondrial sequence, of 16,569 bases, and `chrM` and `M` the sequence of the UCSC hg19 naming, of 16,571 bases, whose positions do not stand at the assembly's ([the two sequences](plugin-contig-mitochondrial-spelling.md)). A table that names its contig `M` and counts on that sequence holds, at a record's position on `MT`, the row of a neighbouring base. The module asks for that row by the name `M`, and its match never compares the row's reference allele with the record's, so nothing stops the neighbour's value. The row it prints carries the record's position and another base's value, with no warning.

**What vep-rs does instead.** vep-rs tries the module's name `M` first, as the module does. On GRCh37 it reads a contig named `M` or `chrM` for a record on `MT` or `chrMT` only after [the check of the table against the run's reference sequence](plugin-contig-mitochondrial-spelling.md): the table is read when its rows state the bases `--fasta` holds at their own positions. A table whose rows state another sequence's bases is not read for the record, and the field stays empty. A record written `chrM` reads the table as it stands, in both engines.

Without `--fasta` the comparison cannot be made. vep-rs then reads no contig of the other pair of names, writes no field for a record on `MT` and warns once, naming the table; Ensembl VEP writes the row's value. That difference is vep-rs's own bound and no part of the class: pass `--fasta`.

**Where it is in Ensembl VEP's source** (`VEP_plugins` at release/116, `ensembl-variation` at release/116):

- `dbNSFP.pm`, line 313: a record whose chromosome matches `MT` in any letter case is named `M`; line 314 reads the table under that name.
- `dbNSFP.pm`, lines 316-325: on GRCh37 the row's position is its `hg19_pos(1-based)` cell; lines 72-74 document the GRCh37 preparation, sorted and indexed on the `hg19_chr` and `hg19_pos(1-based)` columns.
- `dbNSFP.pm`, lines 335-338: a row matches when that position is the record's start and its `alt` cell is the record's allele. No line compares the row's `ref` cell.
- `dbNSFP.pm`, lines 340-344: under `pep_match`, the default, the row's `aaref` and `aaalt` must be the row's amino-acid change.
- `Bio/EnsEMBL/Variation/Utils/BaseVepTabixPlugin.pm`, lines 693-708: the name `M` is tried as written, then as `chrM`.
- `gnomADMt.pm`, lines 256-261, names a mitochondrial record `M` in the same way. Its header names a GRCh38 file (line 50), vep-rs supports the plugin on GRCh38 alone, and the class counts no gnomADMt field.

**What correcting it would change.** A record on `MT` would carry a dbNSFP value only from a table whose mitochondrial rows count on the assembly's sequence.

**The class it produces.** `dbnsfp_mitochondrial_record_named_m`: a dbNSFP field that the reference alone fills on a row of a record spelt `MT` or `chrMT`, in a run on GRCh37. The comparison reads the two outputs and the run's manifest or command line. A field is dbNSFP's when its name opens with the plugin's `dbNSFP_` prefix or, in a run without `output_prefix`, is a column the run's `--plugin dbNSFP` words name. The columns of an `ALL` run are not read.

The comparison can tell the assembly, the record's spelling and which output fills the field. It cannot tell the name the table gives its contig, whether the table's rows state the assembly's bases, or whether the run passed `--fasta`. Nor can it tell a selected column from a field of that name another source of the run writes. A comparison case therefore names the class only for a run whose table counts on another sequence and whose command passes `--fasta`.

**How often the shape occurs.** Tier 3 under the rule of [Intended divergences from Ensembl VEP](../intended-divergences.md#differences-outside-the-classes). The shape needs a GRCh37 run of dbNSFP, a record spelt `MT` or `chrMT` with a missense, stop-gained, stop-lost or start-lost row, and a table that names the contig `M` or `chrM` on another sequence than the assembly's. Five committed flag sets run dbNSFP on GRCh37, and no committed corpus carries the shape: the GRCh37 dbNSFP corpus holds no mitochondrial record.

## 1. A made-up table whose rows state the bases one position on

**In short.** `crates/vep-cli/tests/fixtures/mitochondrial-frame-check` holds a made-up reference of 16,569 bases named `MT` and a made-up dbNSFP-shaped table, `dbNSFP4.9a_one_position_on_synthetic.tsv.gz`, that names its contig `M`. Each of its rows states, as its reference allele, the reference's base one position on. The row at 3577 with the alternate allele `C` states `G`, the reference's base at 3578, and carries `0.14313`, the score the fixture derives for position 3578.

The record is `MT 3577 other_next_missense A C`, a missense variant on `ENST00000361390`. Under `--plugin dbNSFP,<table>,synthetic_score` on GRCh37, Ensembl VEP writes the transcript's row as `other_next_missense MT:3577 C ENST00000361390 missense_variant IMPACT=MODERATE;STRAND=1;synthetic_score=0.14313` (the identifier, location, allele, feature, consequence and extra columns): the value of a row whose reference base is not the record's. vep-rs, given `--fasta` on the fixture's reference, writes the same row ending `IMPACT=MODERATE;STRAND=1`.

**Try it yourself.** Run both engines on the record with the flag above and `--fasta <fixture>/reference_synthetic.fa`, under the caches of [Reproducing these examples](reproducing.md#reproducing-these-examples); the cache slice `crates/vep-cli/tests/fixtures/bayesdel-mitochondrial-contig/json_cache` serves vep-rs. The fixture's `README.md` says how its files were made.
