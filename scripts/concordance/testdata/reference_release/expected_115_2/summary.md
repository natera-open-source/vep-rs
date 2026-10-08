# Concordance Summary

## Aggregate

- Perl tuples: 9
- Rust tuples: 14
- Intersection: 4
- Precision: 0.285714
- Recall: 0.444444
- F1: 0.347826
- Divergence buckets: sv_transcript_selection: 5 (3 excluded), intronic_overcall_swap: 1 (0 excluded), splice_lastwrite_swap: 1 (1 excluded), start_cooccurrence_swap: 1 (1 excluded)
- **Adjusted F1: 0.500000** (excluded 2P + 5R intended-divergence tuples)
- Context-collapsed F1: 0.347826 (collapsed terms: NMD_transcript_variant, coding_transcript_variant, non_coding_transcript_variant)

## Per File

| File | Perl Tuples | Rust Tuples | Intersection | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| output.txt | 9 | 14 | 4 | 0.285714 | 0.444444 | 0.347826 |

## Per File (Context-Collapsed)

| File | Precision | Recall | F1 | Missing | Extra |
|---|---:|---:|---:|---:|---:|
| output.txt | 0.285714 | 0.444444 | 0.347826 | 5 | 10 |
