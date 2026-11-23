# SV Concordance Report

## Overall

- **F1 (adjusted)**: 0.8000  *(excluded 4 Rust-extra + 2 Perl-extra transcript divergences)*
- **F1 (raw)**: 0.5714
- **Total input variants**: 8
- **Perl tuples**: 8
- **Rust tuples**: 13
- **Intersection**: 6
- **Only Perl**: 2
- **Only Rust**: 7

## Per-File Metrics

| File | Input | Perl Tuples | Rust Tuples | Intersect | Raw F1 | Adj F1 |
|------|-------|-------------|-------------|-----------|--------|--------|
| sv | 8 | 8 | 13 | 6 | 0.5714 | 0.8000 |

## Per-Variant-Type Metrics

| Variant Type | Input | Perl Tuples | Rust Tuples | Only Perl | Only Rust | F1 |
|--------------|-------|-------------|-------------|-----------|-----------|------|
| BND | 2 | 4 | 5 | 0 | 1 | 0.8889 ** |
| Complex_SV | 2 | 1 | 4 | 0 | 3 | 0.4000 ** |
| NON_REF | 1 | 1 | 1 | 1 | 1 | 0.0000 |
| Symbolic_DEL | 3 | 2 | 3 | 1 | 2 | 0.4000 ** |

