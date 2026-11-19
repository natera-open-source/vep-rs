# Concordance Summary

## Aggregate

- Perl tuples: 15
- Rust tuples: 20
- Intersection: 8
- Precision: 0.400000
- Recall: 0.533333
- F1: 0.457143
- Divergence buckets: intronic_overcall_swap: 3 (1 excluded), frameshift_splice_cascade_swap: 1 (0 excluded), start_cooccurrence_swap: 1 (1 excluded), sv_transcript_selection: 1 (1 excluded)
- **Adjusted F1: 0.533333** (excluded 2P + 3R intended-divergence tuples)
- Context-collapsed F1: 0.514286 (collapsed terms: NMD_transcript_variant, coding_transcript_variant, non_coding_transcript_variant)

## Per File

| File | Perl Tuples | Rust Tuples | Intersection | Precision | Recall | F1 |
|---|---:|---:|---:|---:|---:|---:|
| sample_a.txt | 10 | 14 | 6 | 0.428571 | 0.600000 | 0.500000 |
| sample_b.txt | 5 | 6 | 2 | 0.333333 | 0.400000 | 0.363636 |

## Per File (Context-Collapsed)

| File | Precision | Recall | F1 | Missing | Extra |
|---|---:|---:|---:|---:|---:|
| sample_a.txt | 0.428571 | 0.600000 | 0.500000 | 4 | 8 |
| sample_b.txt | 0.500000 | 0.600000 | 0.545455 | 2 | 3 |
