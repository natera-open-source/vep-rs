#!/usr/bin/env python3
"""Classify discordant tuples into actionable categories.

Supports both the SV concordance format and the `compare_vep_outputs.py`
tuple discordance format.

Classifies each discordant into one of:

1. coordinate_only: Same start, same transcript, same consequence, only the end
   position differs (VCF parser SVLEN/END calculation).
2. consequence_swap: Same location+transcript, different consequence set
   (consequence logic).
3. transcript_selection: Same variant location+allele, different transcripts
   (transcript lookup or filtering).
4. true_missing/true_extra: Variant not annotated at all by the other tool.

Output: discordant_classified.tsv with a mismatch_type column added.
Also prints a summary report to stdout.
"""

import argparse
import csv
import collections
import sys
from pathlib import Path


def parse_start(loc: str) -> str:
    """Extract start position from 'chr:start-end' or 'chr:pos'."""
    parts = loc.split(':')
    if len(parts) < 2:
        return loc
    coords = parts[1]
    return coords.split('-', 1)[0]


def _detect_schema(fieldnames: list[str]) -> tuple[str, str, str]:
    """Return (source_col, source_file_col, consequence_col) for known schemas."""
    source_col = 'direction' if 'direction' in fieldnames else 'source'
    source_file_col = 'source_file' if 'source_file' in fieldnames else 'file_name'
    consequence_col = 'consequence' if 'consequence' in fieldnames else 'consequence_set'

    required = {'location', 'allele', 'feature', source_col, source_file_col, consequence_col}
    missing = sorted(required - set(fieldnames))
    if missing:
        raise ValueError(
            f"Unsupported discordant TSV schema; missing required columns: {', '.join(missing)}"
        )

    return source_col, source_file_col, consequence_col


def _is_perl_only(row: dict[str, str], source_col: str) -> bool:
    """Return True when the row exists only on the Perl side."""
    return row[source_col] in {'only_perl', 'missing_in_rust'}


def classify_discordants(input_path: str, output_path: str):
    """Read discordant.tsv, classify each row, write classified output."""

    # Load all rows and separate by direction
    rows = []
    with open(input_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = list(reader.fieldnames)
        source_col, source_file_col, consequence_col = _detect_schema(fieldnames)
        for row in reader:
            row['_start'] = parse_start(row['location'])
            rows.append(row)

    perl_rows = [(i, r) for i, r in enumerate(rows) if _is_perl_only(r, source_col)]
    rust_rows = [(i, r) for i, r in enumerate(rows) if not _is_perl_only(r, source_col)]

    # === Pass 1: Coordinate-only matches ===
    # Same (start, allele, feature, consequence) but different location end
    # Build multimap for Rust side, consume from Perl side

    # Key: (start, allele, feature, consequence) -> list of (global_idx, location)
    rust_coord_index = collections.defaultdict(list)
    for gi, row in rust_rows:
        key = (row['_start'], row['allele'], row['feature'], row[consequence_col])
        rust_coord_index[key].append((gi, row['location']))

    # Track which global indices are matched
    matched = set()
    # For each Rust key, track consumption pointer
    rust_consumed = collections.defaultdict(int)

    for gi, row in perl_rows:
        key = (row['_start'], row['allele'], row['feature'], row[consequence_col])
        candidates = rust_coord_index.get(key, [])
        ptr = rust_consumed[key]
        while ptr < len(candidates):
            rgi, rloc = candidates[ptr]
            if rgi not in matched and rloc != row['location']:
                # Coordinate-only match found
                matched.add(gi)
                matched.add(rgi)
                rust_consumed[key] = ptr + 1
                break
            ptr += 1
        else:
            rust_consumed[key] = ptr

    # === Pass 2: Consequence swaps ===
    # Same (location, allele, feature), different consequence
    rust_cons_index = collections.defaultdict(list)
    for gi, row in rust_rows:
        if gi not in matched:
            key = (row['location'], row['allele'], row['feature'])
            rust_cons_index[key].append(gi)

    rust_cons_consumed = collections.defaultdict(int)

    for gi, row in perl_rows:
        if gi in matched:
            continue
        key = (row['location'], row['allele'], row['feature'])
        candidates = rust_cons_index.get(key, [])
        ptr = rust_cons_consumed[key]
        if ptr < len(candidates):
            rgi = candidates[ptr]
            matched.add(gi)
            matched.add(rgi)
            rust_cons_consumed[key] = ptr + 1
            # Mark both as consequence_swap
            rows[gi]['mismatch_type'] = 'consequence_swap'
            rows[rgi]['mismatch_type'] = 'consequence_swap'

    # === Pass 3: Transcript selection ===
    # Same start+allele exists in both directions but different features
    perl_features_by_start_allele = collections.defaultdict(set)
    rust_features_by_start_allele = collections.defaultdict(set)
    for gi, row in perl_rows:
        perl_features_by_start_allele[(row['_start'], row['allele'])].add(row['feature'])
    for gi, row in rust_rows:
        rust_features_by_start_allele[(row['_start'], row['allele'])].add(row['feature'])

    # === Assign classifications ===
    for i, row in enumerate(rows):
        if 'mismatch_type' in row:
            continue  # Already classified as consequence_swap

        if i in matched:
            row['mismatch_type'] = 'coordinate_only'
            continue

        start_allele = (row['_start'], row['allele'])
        if _is_perl_only(row, source_col):
            other_feats = rust_features_by_start_allele.get(start_allele, set())
            if other_feats and row['feature'] not in other_feats:
                row['mismatch_type'] = 'transcript_selection'
            elif other_feats:
                # The same feature exists in Rust but was not matched: a consequence plus
                # coordinate combination is the likely cause
                row['mismatch_type'] = 'true_missing'
            else:
                row['mismatch_type'] = 'true_missing'
        else:
            other_feats = perl_features_by_start_allele.get(start_allele, set())
            if other_feats and row['feature'] not in other_feats:
                row['mismatch_type'] = 'transcript_selection'
            elif other_feats:
                row['mismatch_type'] = 'true_extra'
            else:
                row['mismatch_type'] = 'true_extra'

    # Write output
    output_fields = fieldnames + ['mismatch_type']
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=output_fields, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    # Summary report
    type_counts = collections.Counter()
    type_dir_counts = collections.Counter()
    file_type_counts = collections.Counter()

    for row in rows:
        mt = row.get('mismatch_type', 'unclassified')
        d = row[source_col]
        sf = row[source_file_col]
        type_counts[mt] += 1
        type_dir_counts[(mt, d)] += 1
        file_type_counts[(sf, mt)] += 1

    total = len(rows)
    print("=" * 70)
    print("DISCORDANT CLASSIFICATION SUMMARY")
    print("=" * 70)
    print(f"\nTotal discordants: {total}")
    print(f"\nBy mismatch type:")
    for mt, count in type_counts.most_common():
        pct = 100 * count / total
        print(f"  {mt:25} {count:>8} ({pct:5.1f}%)")

    print(f"\nBy mismatch type + direction:")
    for (mt, d), count in sorted(type_dir_counts.items(), key=lambda x: -x[1]):
        print(f"  {mt:25} {d:12} {count:>8}")

    print(f"\nBy file + mismatch type (top 30):")
    items = sorted(file_type_counts.items(), key=lambda x: -x[1])
    for (sf, mt), count in items[:30]:
        pct = 100 * count / total
        print(f"  {sf:35} {mt:25} {count:>8} ({pct:5.1f}%)")


def main():
    parser = argparse.ArgumentParser(description='Classify discordant tuples')
    parser.add_argument('input', help='Path to discordant.tsv')
    parser.add_argument('-o', '--output', help='Output path (default: discordant_classified.tsv)',
                       default=None)
    args = parser.parse_args()

    if args.output is None:
        args.output = str(Path(args.input).parent / 'discordant_classified.tsv')

    classify_discordants(args.input, args.output)


if __name__ == '__main__':
    main()
