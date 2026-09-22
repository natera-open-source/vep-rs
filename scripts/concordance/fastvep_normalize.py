#!/usr/bin/env python3
"""Normalize fastVEP VEP-default-tab output so it matches Perl VEP semantics.

fastVEP's VEP-default-tab writer emits `Feature_type=Transcript` on rows
where `Feature=-` (i.e. intergenic). Perl VEP emits `Feature_type=-` on
those same rows. The concordance comparator keys on
(location, allele, feature, feature_type, consequence), so without
normalization every intergenic row appears as a discordant tuple.

The Feature_type rewrite is not the only tab-output difference between the
two engines. Differences in the column COUNT (`DISTANCE` column presence,
etc.) live outside the 5-tuple key and do not affect concordance; a
CONSEQUENCE-VOCABULARY difference does:

fastVEP emits VCF 4.4 copy-number terms that VEP release 115 never emits
(`copy_number_increase`, `copy_number_decrease`, `copy_number_change`,
`transcript_variant`; no VEP tuple carries any of them). It emits them
ALONGSIDE the term VEP uses, so `{copy_number_increase, feature_elongation}`
cannot match `{feature_elongation}` on a key that includes the whole
consequence set. Per dataset, its `feature_elongation` and
`copy_number_increase` counts are equal (174 on the GRCh37 1000 Genomes
Phase 3 chr21 dataset, 1,774 on the GRCh38 1000 Genomes high-coverage
dataset; `manuscript/data/f1_by_consequence_class_by_dataset.csv`), which
identifies the added term rather than a missing annotation as the cause.

That difference is NOT normalized here, and the distinction from the
Feature_type rewrite above is the reason. Rewriting `Transcript` to `-` on an
intergenic row corrects a FORMATTING difference in a field both engines agree
about. Rewriting a third party's consequence vocabulary would change what is
being measured: the reported quantity is concordance with VEP's output on that
key because the claim is drop-in replacement, and a pipeline whose filters are
written against VEP's consequence sets genuinely gets a different answer from
fastVEP's. The discordance is real and belongs in the result, so it is
reported rather than normalized away. A vocabulary map here would change the
quantity being measured.

Usage:
    python3 fastvep_normalize.py <input.fastvep.tab> <output.vep.tab>
"""

import sys


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: fastvep_normalize.py <in> <out>", file=sys.stderr)
        return 1
    infile, outfile = sys.argv[1], sys.argv[2]
    rewritten = 0
    total = 0
    with open(infile) as f_in, open(outfile, "w") as f_out:
        for line in f_in:
            if not line or line.startswith("#"):
                f_out.write(line)
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            # Perl emits "-" for Feature_type on intergenic; fastVEP emits "Transcript".
            # Fields are 1-based: 1=Uploaded, 2=Location, 3=Allele, 4=Gene, 5=Feature,
            # 6=Feature_type, 7=Consequence. In 0-based slicing: Gene=fields[3],
            # Feature=fields[4], Feature_type=fields[5].
            if (
                len(fields) >= 6
                and fields[3] == "-"
                and fields[4] == "-"
                and fields[5] == "Transcript"
            ):
                fields[5] = "-"
                rewritten += 1
            f_out.write("\t".join(fields) + "\n")
    print(
        f"normalize: rewrote {rewritten}/{total} intergenic Feature_type=Transcript → -",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
