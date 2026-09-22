#!/usr/bin/env python3
"""Compare plugin output fields between Perl VEP and Rust VEP.

Supports batch mode (multiple plugins), per-directory processing (all matching
file pairs), and JSON output for integration with the concordance harness.

Usage:
    # Single plugin, single file pair
    python3 compare_plugin_outputs.py \
        --perl-output perl_out.txt --rust-output rust_out.txt \
        --plugins CADD

    # Multiple plugins, directory mode
    python3 compare_plugin_outputs.py \
        --perl-dir tmp/concordance/perl --rust-dir tmp/concordance/rust \
        --plugins CADD,REVEL,gnomADc \
        --json-output tmp/concordance/reports/plugin_concordance.json

    # All known plugins with custom tolerance
    python3 compare_plugin_outputs.py \
        --perl-dir perl/ --rust-dir rust/ \
        --plugins all --tolerance 1e-4
"""

import argparse
import gzip
import json
import sys
from pathlib import Path
from urllib.parse import unquote

# Plugin field registry

PLUGIN_FIELDS = {
    "HGVS": ["HGVSc", "HGVSp"],
    "CADD": ["CADD_PHRED", "CADD_RAW"],
    "REVEL": ["REVEL_score"],
    "SpliceAI": [
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
        "SpliceAI_pred_DP_AG",
        "SpliceAI_pred_DP_AL",
        "SpliceAI_pred_DP_DG",
        "SpliceAI_pred_DP_DL",
        "SpliceAI_pred_SYMBOL",
    ],
    "gnomADc": [
        "gnomADc_mean",
        "gnomADc_median",
        "gnomADc_over_1",
        "gnomADc_over_5",
        "gnomADc_over_10",
        "gnomADc_over_15",
        "gnomADc_over_20",
        "gnomADc_over_25",
        "gnomADc_over_30",
        "gnomADc_over_50",
        "gnomADc_over_100",
    ],
    "AlphaMissense": ["am_class", "am_pathogenicity"],
    "GWAS": [
        "GWAS_accessions",
        "GWAS_associated_gene",
        "GWAS_beta_coef",
        "GWAS_odds_ratio",
        "GWAS_p_value",
        "GWAS_pmid",
        "GWAS_risk_allele",
        "GWAS_study",
    ],
    "dbNSFP": [],  # dynamic: fields prefixed with dbNSFP_ are auto-discovered
    "dbscSNV": ["dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE"],
    "LoFtool": ["LoFtool_percentile"],
    "pLI": ["pLI_gene_value"],
    # LoFTEE emits the same field names as the Perl LoF plugin, so vep-rs
    # (--plugin LoFTEE) and Perl (--plugin LoF) compare directly. Compare the
    # LoF (HC/LC) verdict and LoF_filter per-field; LoF_flags/LoF_info are
    # informational and differ in formatting (use --fields LoF to isolate).
    "LoFTEE": ["LoF", "LoF_filter", "LoF_flags", "LoF_info"],
}

# `HGVS` is in PLUGIN_FIELDS because this script doubles as the HGVS comparator,
# but it is not a plugin: `--plugins all` leaves it out, because a pseudo-plugin
# evaluated against plugin-run outputs can only report zero.
NON_PLUGIN_FIELD_SETS = {"HGVS"}
ALL_PLUGIN_NAMES = [k for k in PLUGIN_FIELDS if k not in NON_PLUGIN_FIELD_SETS]


def parse_extra_field(extra_str):
    """Parse VEP Extra column (key=value;key=value) into a dict."""
    if extra_str == "-" or not extra_str:
        return {}
    result = {}
    for pair in extra_str.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            result[key] = value
    return result


def is_vcf_path(filepath):
    lower = str(filepath).lower()
    return lower.endswith(".vcf") or lower.endswith(".vcf.gz")


def open_text_auto(path):
    """Open plain or gzip text file in UTF-8 with replacement mode."""
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def decode_csq_value(value):
    """Decode VCF CSQ percent-encoding (VEP-compatible)."""
    return unquote(value)


def parse_csq_header_fields(filepath):
    """Parse CSQ field names from a VCF header."""
    with open_text_auto(filepath) as f:
        for line in f:
            if not line.startswith("##INFO=<ID=CSQ"):
                continue
            marker = "Format: "
            idx = line.find(marker)
            if idx == -1:
                return []
            tail = line[idx + len(marker) :]
            # Header line is quoted: ... Format: A|B|C">
            fields = tail.split('"')[0]
            return fields.split("|")
    return []


def extract_info_value(info_str, key):
    """Extract key=value from VCF INFO; returns None if missing."""
    prefix = f"{key}="
    for field in info_str.split(";"):
        if field.startswith(prefix):
            return field[len(prefix) :]
    return None


def discover_plugin_fields(filepath, plugin_name):
    """Auto-discover fields for dynamic plugins (e.g., dbNSFP_*)."""
    prefix = f"{plugin_name}_"
    found = set()
    if is_vcf_path(filepath):
        csq_fields = parse_csq_header_fields(filepath)
        for field in csq_fields:
            if field.startswith(prefix):
                found.add(field)
        return sorted(found)

    with open_text_auto(filepath) as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 14:
                continue
            extra = parse_extra_field(parts[13])
            for key in extra:
                if key.startswith(prefix):
                    found.add(key)
            # Sample first 1000 data lines for discovery
            if i > 1000:
                break
    return sorted(found)


def parse_vep_output(filepath, plugin_fields, all_prefix_fields=None):
    """Parse a VEP text output file and extract plugin fields.

    Returns a dict mapping (location, allele, feature) -> {field: value}.
    If all_prefix_fields is set, also collects any field starting with that prefix.
    """
    records = {}
    with open_text_auto(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 14:
                continue

            location = parts[1]
            allele = parts[2]
            feature = parts[4]
            extra = parse_extra_field(parts[13])

            key = (location, allele, feature)
            plugin_data = {}

            # Collect registered fields
            for field in plugin_fields:
                if field in extra:
                    plugin_data[field] = extra[field]

            # Collect prefix-matched fields (for dbNSFP etc.)
            if all_prefix_fields:
                for field_name, value in extra.items():
                    if any(field_name.startswith(p) for p in all_prefix_fields):
                        plugin_data[field_name] = value

            if plugin_data:
                records[key] = plugin_data

    return records


def parse_vep_vcf_output(filepath, plugin_fields, all_prefix_fields=None):
    """Parse VCF output and extract plugin fields from CSQ entries.

    Returns a dict mapping (location, allele, feature) -> {field: value}.
    """
    records = {}
    csq_fields = parse_csq_header_fields(filepath)
    if not csq_fields:
        return records

    with open_text_auto(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom = parts[0]
            pos = parts[1]
            info = parts[7]
            csq_value = extract_info_value(info, "CSQ")
            if not csq_value:
                continue

            entries = csq_value.split(",")
            for entry in entries:
                raw_values = entry.split("|")
                # Preserve empty trailing columns by padding to header width.
                if len(raw_values) < len(csq_fields):
                    raw_values += [""] * (len(csq_fields) - len(raw_values))
                elif len(raw_values) > len(csq_fields):
                    raw_values = raw_values[: len(csq_fields)]

                decoded_values = [decode_csq_value(v) for v in raw_values]
                csq_map = dict(zip(csq_fields, decoded_values))

                location = f"{chrom}:{pos}"
                allele = csq_map.get("Allele", "")
                feature = csq_map.get("Feature", "")
                key = (location, allele, feature)

                plugin_data = {}
                for field in plugin_fields:
                    if field in csq_map and csq_map[field] not in ("", "."):
                        plugin_data[field] = csq_map[field]

                if all_prefix_fields:
                    for field_name, value in csq_map.items():
                        if value in ("", "."):
                            continue
                        if any(field_name.startswith(p) for p in all_prefix_fields):
                            plugin_data[field_name] = value

                if plugin_data:
                    records[key] = plugin_data

    return records


def compare_records(perl_records, rust_records, plugin_fields, tolerance):
    """Compare plugin fields between Perl and Rust outputs for a single plugin."""
    all_keys = set(perl_records.keys()) | set(rust_records.keys())

    stats = {
        "total_keys": len(all_keys),
        "match": 0,
        "mismatch": 0,
        "perl_only": 0,
        "rust_only": 0,
        "both_empty": 0,
    }
    mismatches = []

    for key in sorted(all_keys):
        perl_data = perl_records.get(key, {})
        rust_data = rust_records.get(key, {})

        if not perl_data and not rust_data:
            stats["both_empty"] += 1
            continue

        if perl_data and not rust_data:
            stats["perl_only"] += 1
            mismatches.append(
                {
                    "key": list(key),
                    "type": "perl_only",
                    "perl": dict(perl_data),
                    "rust": {},
                }
            )
            continue

        if rust_data and not perl_data:
            stats["rust_only"] += 1
            mismatches.append(
                {
                    "key": list(key),
                    "type": "rust_only",
                    "perl": {},
                    "rust": dict(rust_data),
                }
            )
            continue

        # Both have data -- compare field values
        # Use union of fields present in either side
        check_fields = set(perl_data.keys()) | set(rust_data.keys())
        if plugin_fields:
            check_fields &= set(plugin_fields)

        field_match = True
        for field in check_fields:
            perl_val = perl_data.get(field, "")
            rust_val = rust_data.get(field, "")
            if perl_val == rust_val:
                continue
            # Numeric tolerance
            try:
                if abs(float(perl_val) - float(rust_val)) < tolerance:
                    continue
            except (ValueError, TypeError):
                pass
            field_match = False
            break

        if field_match:
            stats["match"] += 1
        else:
            stats["mismatch"] += 1
            mismatches.append(
                {
                    "key": list(key),
                    "type": "mismatch",
                    "perl": dict(perl_data),
                    "rust": dict(rust_data),
                }
            )

    annotated = (
        stats["match"] + stats["mismatch"] + stats["perl_only"] + stats["rust_only"]
    )
    # A zero-annotation run is NOT 100% concordant: it is unmeasured. Reporting
    # 100.0 here would let a completely dead plugin read as a pass: a contig-name
    # mismatch, an unparsed header, or a wrong position column all produce zero
    # records and no error. `None` forces callers to treat it as "no data".
    stats["concordance_pct"] = (
        (stats["match"] / annotated * 100) if annotated > 0 else None
    )
    stats["annotated_keys"] = annotated
    # Per-ENGINE annotated counts. The union count above cannot distinguish "both
    # engines annotated the same 12 records" from "only one engine annotated
    # anything", and a one-sided run is a failure even when every shared key
    # agrees.
    stats["perl_annotated"] = stats["match"] + stats["mismatch"] + stats["perl_only"]
    stats["rust_annotated"] = stats["match"] + stats["mismatch"] + stats["rust_only"]

    return stats, mismatches


def process_file_pair(perl_path, rust_path, plugins, tolerance):
    """Process a single Perl/Rust output file pair for all requested plugins."""
    results = {}

    for plugin_name in plugins:
        registered_fields = PLUGIN_FIELDS.get(plugin_name, [])
        prefix_fields = []

        # For dynamic plugins (empty registered fields), auto-discover
        if not registered_fields:
            prefix_fields = [f"{plugin_name}_"]
            # Discover from both files
            discovered_perl = discover_plugin_fields(perl_path, plugin_name)
            discovered_rust = discover_plugin_fields(rust_path, plugin_name)
            registered_fields = sorted(set(discovered_perl) | set(discovered_rust))

        # Also add gnomADc prefix matching for auto-discovered fields
        if plugin_name == "gnomADc":
            prefix_fields = ["gnomADc_"]

        if is_vcf_path(perl_path):
            perl_records = parse_vep_vcf_output(
                perl_path, registered_fields, prefix_fields
            )
        else:
            perl_records = parse_vep_output(perl_path, registered_fields, prefix_fields)

        if is_vcf_path(rust_path):
            rust_records = parse_vep_vcf_output(
                rust_path, registered_fields, prefix_fields
            )
        else:
            rust_records = parse_vep_output(rust_path, registered_fields, prefix_fields)

        stats, mismatches = compare_records(
            perl_records, rust_records, registered_fields, tolerance
        )

        results[plugin_name] = {
            "fields": registered_fields,
            "perl_annotated_records": len(perl_records),
            "rust_annotated_records": len(rust_records),
            "stats": stats,
            "mismatches_sample": mismatches[:10],  # Keep first 10 for reporting
            "total_mismatches": len(mismatches),
        }

    return results


def canonical_stem(path_obj):
    """Normalize benchmark file names for perl/rust pair matching."""
    name = path_obj.name
    for suffix in [".vep.vcf.gz", ".vep.vcf", ".vcf.gz", ".vcf", ".txt", ".gz"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.startswith("perl_"):
        name = name[len("perl_") :]
    elif name.startswith("rust_"):
        name = name[len("rust_") :]
    return name


def candidate_output_files(root):
    root_path = Path(root)
    files = []
    for pattern in ("*.txt", "*.vep.vcf", "*.vep.vcf.gz", "*.vcf", "*.vcf.gz"):
        files.extend(root_path.glob(pattern))
    return files


def find_matching_pairs(perl_dir, rust_dir):
    """Find matching file pairs in perl and rust directories."""
    perl_files = {canonical_stem(f): f for f in candidate_output_files(perl_dir)}
    rust_files = {canonical_stem(f): f for f in candidate_output_files(rust_dir)}
    common = sorted(set(perl_files.keys()) & set(rust_files.keys()))
    return [(str(perl_files[name]), str(rust_files[name]), name) for name in common]


def print_report(all_results, plugins, min_annotated=1):
    """Print human-readable summary to stdout.

    `min_annotated` is the per-engine floor of annotated records below which a
    plugin FAILS. It defaults to 1 because the interesting failure is zero: a
    contig-name mismatch, an unparsed header, or a wrong position column all make
    a plugin annotate nothing while exiting 0.
    """
    any_failure = False

    for plugin_name in plugins:
        print(f"\n{'=' * 60}")
        print(f"  Plugin: {plugin_name}")
        print(f"{'=' * 60}")

        # Aggregate across files
        total_match = 0
        total_mismatch = 0
        total_perl_only = 0
        total_rust_only = 0
        total_annotated = 0
        total_perl_annotated = 0
        total_rust_annotated = 0

        for file_name, file_results in all_results.items():
            if plugin_name not in file_results:
                continue
            pr = file_results[plugin_name]
            s = pr["stats"]
            total_match += s["match"]
            total_mismatch += s["mismatch"]
            total_perl_only += s["perl_only"]
            total_rust_only += s["rust_only"]
            total_annotated += s["annotated_keys"]
            total_perl_annotated += s.get("perl_annotated", 0)
            total_rust_annotated += s.get("rust_annotated", 0)

            print(f"\n  File: {file_name}")
            print(
                f"    Fields: {', '.join(pr['fields'][:5])}{'...' if len(pr['fields']) > 5 else ''}"
            )
            print(f"    Perl records: {pr['perl_annotated_records']}")
            print(f"    Rust records: {pr['rust_annotated_records']}")
            print(
                f"    Match: {s['match']}  Mismatch: {s['mismatch']}  "
                f"Perl-only: {s['perl_only']}  Rust-only: {s['rust_only']}"
            )
            if s["concordance_pct"] is None:
                print("    Concordance: n/a (no records annotated by either engine)")
            else:
                print(f"    Concordance: {s['concordance_pct']:.2f}%")

            if pr["mismatches_sample"]:
                print("    First mismatches:")
                for m in pr["mismatches_sample"][:3]:
                    loc, allele, feat = m["key"]
                    print(f"      {loc} {allele} {feat}: {m['type']}")
                    if m["perl"]:
                        print(f"        Perl: {m['perl']}")
                    if m["rust"]:
                        print(f"        Rust: {m['rust']}")

        if total_annotated > 0:
            agg = f"{total_match / total_annotated * 100:.2f}%"
        else:
            agg = "n/a"

        print(f"\n  AGGREGATE: {total_match}/{total_annotated} match ({agg})")
        print(
            f"  ANNOTATED: perl={total_perl_annotated}  rust={total_rust_annotated}"
            f"  (floor {min_annotated})"
        )

        # Reasons a plugin FAILS. Each is a shape that would otherwise pass silently:
        #   - either engine annotated fewer than `min_annotated` records: the run
        #     measured nothing, which is not agreement;
        #   - rust_only > 0: vep-rs over-annotated where Perl did not.
        reasons = []
        if total_perl_annotated < min_annotated:
            reasons.append(
                f"Perl annotated {total_perl_annotated} records (< {min_annotated})"
            )
        if total_rust_annotated < min_annotated:
            reasons.append(
                f"vep-rs annotated {total_rust_annotated} records (< {min_annotated})"
            )
        if total_mismatch > 0:
            reasons.append(f"{total_mismatch} value mismatches")
        if total_perl_only > 0:
            reasons.append(f"{total_perl_only} Perl-only records")
        if total_rust_only > 0:
            reasons.append(f"{total_rust_only} vep-rs-only records")

        if reasons:
            print(f"  STATUS: FAIL ({'; '.join(reasons)})")
            any_failure = True
        else:
            print("  STATUS: PASS")

    return any_failure


def main():
    parser = argparse.ArgumentParser(
        description="Compare VEP plugin outputs between Perl and Rust"
    )
    # Input modes: single file pair or directory
    parser.add_argument("--perl-output", help="Single Perl VEP output file")
    parser.add_argument("--rust-output", help="Single Rust VEP output file")
    parser.add_argument("--perl-dir", help="Directory of Perl VEP output files")
    parser.add_argument("--rust-dir", help="Directory of Rust VEP output files")

    # Plugin selection
    parser.add_argument(
        "--plugins",
        required=True,
        help="Comma-separated plugin names, or 'all' for all known plugins",
    )
    parser.add_argument(
        "--fields",
        default=None,
        help="Override: comma-separated field names (for single-plugin manual use)",
    )

    # Output
    parser.add_argument(
        "--json-output", default=None, help="Write JSON results to this file"
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=10,
        help="Max mismatches to include in JSON output per plugin per file (default: 10)",
    )

    # Tolerance
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="Numeric tolerance for floating-point comparison (default: 1e-6)",
    )
    parser.add_argument(
        "--min-annotated",
        type=int,
        default=1,
        help=(
            "Minimum annotated records EACH engine must produce, or the plugin "
            "FAILS (default: 1). A plugin that annotates nothing has not been "
            "measured -- a contig-name mismatch, an unparsed header, or a wrong "
            "position column all yield zero records without an error."
        ),
    )

    args = parser.parse_args()

    # Resolve plugins
    if args.plugins.lower() == "all":
        plugins = ALL_PLUGIN_NAMES
    else:
        plugins = [p.strip() for p in args.plugins.split(",")]
        unknown = [p for p in plugins if p not in PLUGIN_FIELDS]
        if unknown:
            # A typo would otherwise fall through to `PLUGIN_FIELDS.get(name, [])`,
            # which is the same signal as an intentionally dynamic plugin, so the
            # run would "pass" having compared nothing.
            print(
                f"ERROR: unknown plugin(s): {', '.join(unknown)}. "
                f"Known: {', '.join(sorted(PLUGIN_FIELDS))}",
                file=sys.stderr,
            )
            return 2

    # Override fields for single-plugin manual use
    if args.fields and len(plugins) == 1:
        PLUGIN_FIELDS[plugins[0]] = [f.strip() for f in args.fields.split(",")]

    # Resolve file pairs
    if args.perl_dir and args.rust_dir:
        pairs = find_matching_pairs(args.perl_dir, args.rust_dir)
        if not pairs:
            print(
                f"ERROR: no matching file pairs found in {args.perl_dir} and {args.rust_dir}",
                file=sys.stderr,
            )
            sys.exit(1)
    elif args.perl_output and args.rust_output:
        name = Path(args.perl_output).stem
        pairs = [(args.perl_output, args.rust_output, name)]
    else:
        print(
            "ERROR: provide either --perl-dir/--rust-dir or --perl-output/--rust-output",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Plugins: {', '.join(plugins)}")
    print(f"File pairs: {len(pairs)}")
    print(f"Tolerance: {args.tolerance}")

    # Process all file pairs
    all_results = {}
    for perl_path, rust_path, file_name in pairs:
        all_results[file_name] = process_file_pair(
            perl_path, rust_path, plugins, args.tolerance
        )

    # Print human-readable report
    any_failure = print_report(all_results, plugins, args.min_annotated)

    # Write JSON output
    if args.json_output:
        json_data = {
            "plugins": plugins,
            "tolerance": args.tolerance,
            "min_annotated": args.min_annotated,
            "file_count": len(pairs),
            "results": all_results,
        }

        # Add per-plugin aggregate stats
        aggregates = {}
        for plugin_name in plugins:
            agg = {
                "match": 0,
                "mismatch": 0,
                "perl_only": 0,
                "rust_only": 0,
                "annotated": 0,
                "perl_annotated": 0,
                "rust_annotated": 0,
            }
            for file_results in all_results.values():
                if plugin_name in file_results:
                    s = file_results[plugin_name]["stats"]
                    agg["match"] += s["match"]
                    agg["mismatch"] += s["mismatch"]
                    agg["perl_only"] += s["perl_only"]
                    agg["rust_only"] += s["rust_only"]
                    agg["annotated"] += s["annotated_keys"]
                    agg["perl_annotated"] += s.get("perl_annotated", 0)
                    agg["rust_annotated"] += s.get("rust_annotated", 0)
            # None, not 100.0, when nothing was annotated: unmeasured is not perfect.
            agg["concordance_pct"] = (
                agg["match"] / agg["annotated"] * 100 if agg["annotated"] > 0 else None
            )
            # Mirrors print_report's predicate: both engines must clear the floor,
            # and rust_only counts against the run.
            agg["pass"] = (
                agg["mismatch"] == 0
                and agg["perl_only"] == 0
                and agg["rust_only"] == 0
                and agg["perl_annotated"] >= args.min_annotated
                and agg["rust_annotated"] >= args.min_annotated
            )
            aggregates[plugin_name] = agg

        json_data["aggregates"] = aggregates

        Path(args.json_output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_output).write_text(
            json.dumps(json_data, indent=2), encoding="utf-8"
        )
        print(f"\nJSON results written to: {args.json_output}")

    # Exit code
    if any_failure:
        print("\nOVERALL: FAIL (one or more plugins have discordances)")
        sys.exit(1)
    else:
        print("\nOVERALL: PASS")
        sys.exit(0)


if __name__ == "__main__":
    # main() returns a status for its early-exit paths (e.g. an unknown plugin
    # name returns 2). Called bare, that status is discarded: a typo'd --plugins
    # prints ERROR yet exits 0, and a script gating on exit status reads a
    # misconfigured run as a pass. The success/failure paths call sys.exit()
    # directly, so this only ever propagates an early return.
    sys.exit(main() or 0)
