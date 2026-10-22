# A breakend mate on another chromosome under `--check_svs`

A difference in the text of the `SV` field, on a breakend whose bracket mate names another chromosome, under a flag no measured run names; nothing is set aside.

Under `--check_svs` Ensembl VEP marks each record with the structural variants of its database that overlap it. The variants fetched for a buffer are those overlapping its records' window on their chromosome and, for a bracket breakend, those overlapping the mate's position on the mate's chromosome (`AnnotationSource.pm` lines 249-257 add the mate's chromosome and position to the regions fetched).

For every variant fetched it asks the buffer for the records overlapping the variant's span, with no chromosome (`AnnotationSource/Database/StructuralVariation.pm` line 87). The buffer's test then accepts a record on any chromosome (`InputBuffer.pm` line 351). A bracket breakend sits in the buffer's tree at its mate's coordinate as well as its own (lines 387-389), and is matched there when a fetched variant's span holds the mate's number (lines 339-341), whatever chromosome the variant lies on.

So a breakend on chromosome 21 whose mate is `22:5035000` collects three sets of names. The first two are the variants spanning its own position and the chromosome 22 variants spanning the mate's position. The third is every fetched chromosome 21 variant that also spans position 5035000, a position the record never names on that chromosome.

vep-rs reads the same records from the file `--check_svs` names and matches a bracket mate on the mate's own chromosome at the mate's position, so the breakend collects the first two sets and not the third. The two engines agree on every record whose mate lies on its own chromosome, and on every record without a bracket mate.

The records below are a substitution and, one base after it, a breakend whose mate lies on chromosome 22. The structural variant `nsv6146554` spans 5029000 to 5047000 on chromosome 21, the first record's position and the mate's number alike:

```
21	5047000	.	T	A	.	.	.
21	5047001	.	T	T[22:5035000[	.	.	SVTYPE=BND
```

Both engines write the substitution with six names, the measured Ensembl VEP 116.2 row of the `GRCh38-checksvs` corpus:

```
SV=nsv2770691,nsv3109934,nsv3163927,nsv3167562,nsv533679,nsv6146554
```

For the breakend's two rows (`T.` and `T[22:5035000[`), read together with the substitution in one buffer, Ensembl VEP's lines above give the same six names, `nsv6146554` among them through the mate's number. That row is derived from the source: the corpus holds no breakend record, so no reference row was written for one. The file holds chromosome 21 alone, so no chromosome 22 variant enters either engine's list. vep-rs writes the five that span the breakend's own position:

```
SV=nsv2770691,nsv3109934,nsv3163927,nsv3167562,nsv533679
```

The flag names the structural variants a record overlaps, and a position on another chromosome overlaps no chromosome 21 variant, so the chromosome-blind match is a defect of the lookup and not a reading vep-rs reproduces.

To reproduce in Ensembl VEP 116.2: the flag reads the structural variation database, so run [the command of the closing section](reproducing.md#reproducing-these-examples) with `--cache --db_version 116 --check_svs` in place of `--offline`, the two records in one file. In vep-rs: the same records with `--check_svs tests/golden/116/GRCh38-checksvs/structural_variants.vcf.gz` over the GRCh38 cache.
