# SV Concordance Report

## Overall

- **F1 (adjusted)**: 0.4074  *(excluded 5 Rust-extra + 5 Perl-extra transcript divergences)*
- **F1 (raw)**: 0.3438
- **Total input variants**: 25
- **Perl tuples**: 27
- **Rust tuples**: 37
- **Intersection**: 11
- **Only Perl**: 16
- **Only Rust**: 26

SV adjusted F1 additionally excludes the `<CNV:TR>` literal-expansion class (3 of 9 pairs GRCH37, 03_tandem_repeats only). Raw F1 unaffected.

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|
| 01_breakends | 6 | 9 | 14 | 8 | 0.6957 | 0.7619 |
| 02_one_sided | 12 | 11 | 11 | 1 | 0.0909 | 0.1000 |
| 03_tandem_repeats | 6 | 5 | 10 | 1 | 0.1333 | 0.2222 |
| 04_small_variants | 1 | 2 | 2 | 1 | 0.5000 |  |

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|
| BND | 7 | 11 | 16 | 3 | 8 | 0.5926 ** |
| Multi_allelic | 2 | 2 | 5 | 2 | 5 | 0.0000 |
| RefOnly | 1 | 1 | 0 | 1 | 0 | 0.0000 |
| SNV | 1 | 1 | 0 | 1 | 0 | 0.0000 |
| Small_INS | 1 | 2 | 2 | 1 | 1 | 0.5000 ** |
| Spanning | 1 | 2 | 0 | 2 | 0 | 0.0000 |
| Symbolic_CNV | 4 | 2 | 6 | 2 | 6 | 0.0000 |
| Symbolic_DEL | 8 | 6 | 7 | 5 | 6 | 0.1538 ** |

