#!/usr/bin/env python3
"""Validate ALL test VCFs in the vep-rs SV validation test suite.

Checks:
  1. REF allele correctness (requires FASTA + samtools)
  2. VCF spec compliance via bcftools
  3. Assembly consistency (contig lengths vs declared assembly)
  4. Variant count report per file/assembly
  5. Real-world filename symmetry across assemblies
  6. VCF version header check

Usage:
    python3 scripts/validation/validate_test_vcfs.py \
        --fasta-grch37 /path/to/GRCh37.fa --fasta-grch38 /path/to/GRCh38.fa \
        [--test-dir tests/sv_validation]

    # Without FASTA (skips REF validation, runs all other checks):
    python3 scripts/validation/validate_test_vcfs.py

Exit code 0 if all checks pass, 1 if any fail.
"""

import argparse
import gzip
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Optional

# Constants

EXPECTED_CONTIG_LENGTHS = {
    "GRCh37": {"21": 48129895},
    "GRCh38": {"21": 46709983},
}

# Real-world file basenames expected in BOTH assemblies (minus version
# differences), when the real-world VCFs are present at all.
#
# They are NOT committed: they are large public downloads, seeded by
# `scripts/data/download_real_world_vcfs.sh` into `${VEP_VCF_DIR:-./.vep/vcf}/<asm>/sv/`,
# and `<test_dir>/<asm>/real_world/` is a symlink to (or copy of) that directory
# (tests/sv_validation/README.md). A tree without them is a normal
# state, so their absence is reported as SKIPPED rather than failed; only an
# ASYMMETRIC set (present for one assembly, missing for the other) is an error,
# because a per-assembly comparison over different sources is not a comparison.
EXPECTED_REAL_WORLD_STEMS = {"1kg_sv_chr21", "clinvar_sv_chr21", "gnomad_sv"}

# Bracket-notation breakend ALT.
BND_ALT_RE = re.compile(r"[\[\]]")


# RefFasta -- same pattern as generate_sv_test_vcfs.py


class RefFasta:
    """Fetch reference bases with full-chromosome pre-loading.

    Loads the entire primary contig (chr21, ~48MB) into memory once so all
    lookups are O(1) string slicing with zero subprocess overhead.
    """

    def __init__(self, fasta_path: str, contig: str = "21"):
        self.fasta_path = fasta_path
        self._seq: str = ""
        self._contig = contig
        self._load_contig(contig)

    def _load_contig(self, contig: str) -> None:
        result = subprocess.run(
            ["samtools", "faidx", self.fasta_path, contig],
            capture_output=True,
            text=True,
            check=True,
        )
        self._seq = "".join(result.stdout.strip().split("\n")[1:]).upper()
        print(f"  Loaded {contig}: {len(self._seq):,} bp", file=sys.stderr)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch bases from chrom:start-end (1-based inclusive)."""
        if chrom == self._contig and 1 <= start and end <= len(self._seq):
            return self._seq[start - 1 : end]
        # Fallback for other contigs (e.g. inter-chromosomal BND)
        region = f"{chrom}:{start}-{end}"
        result = subprocess.run(
            ["samtools", "faidx", self.fasta_path, region],
            capture_output=True,
            text=True,
            check=True,
        )
        return "".join(result.stdout.strip().split("\n")[1:]).upper()

    def fetch_base(self, chrom: str, pos: int) -> str:
        return self.fetch(chrom, pos, pos)

    @property
    def contig(self) -> str:
        """The contig loaded into memory."""
        return self._contig


# Data classes for results


@dataclass
class FileResult:
    """Validation results for a single VCF file."""

    path: str
    assembly: str
    category: str  # "synthetic" or "real_world"
    vcf_version: str = ""
    variant_count: int = 0
    bcftools_ok: bool = True
    bcftools_error: str = ""
    contig_ok: bool = True
    contig_errors: list = field(default_factory=list)
    ref_checked: bool = False
    ref_total: int = 0
    ref_mismatches: int = 0
    ref_errors: list = field(default_factory=list)  # first N mismatches

    @property
    def passed(self) -> bool:
        return (
            self.bcftools_ok
            and self.contig_ok
            and (not self.ref_checked or self.ref_mismatches == 0)
        )


# Check 1: REF allele correctness


def _should_check_ref(chrom: str, alt: str, primary_contig: str) -> bool:
    """Whether a record's REF is validated against the FASTA.

    Every record on the primary contig is checked: the full REF of an explicit
    sequence variant, the single anchor base of a symbolic allele, and the REF
    of a ref-only or spanning record. The one exemption is a breakend mate on
    another contig: `generate_sv_test_vcfs.py` loads only the primary contig
    and writes a placeholder base for those mates, so their REF is not a
    reference claim.
    """
    return chrom == primary_contig or not BND_ALT_RE.search(alt)


def check_ref_alleles(
    vcf_path: str, fasta: RefFasta, max_errors: int = 20
) -> tuple[int, int, list[str]]:
    """Check every variant's REF allele against the reference FASTA.

    Returns (total_checked, mismatches, error_details).
    """
    total = 0
    mismatches = 0
    errors: list[str] = []

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3].upper()
            alt = fields[4]

            if not _should_check_ref(chrom, alt, fasta.contig):
                continue

            total += 1

            # Fetch expected REF from FASTA
            try:
                expected = fasta.fetch(chrom, pos, pos + len(ref) - 1)
            except (subprocess.CalledProcessError, ValueError):
                # Can't fetch (e.g. contig not in FASTA) -- skip
                continue

            if ref != expected:
                mismatches += 1
                if len(errors) < max_errors:
                    errors.append(
                        f"  {chrom}:{pos} REF={ref} expected={expected} "
                        f"(len={len(ref)})"
                    )

    return total, mismatches, errors


# Check 2: bcftools spec compliance


def check_bcftools(vcf_path: str) -> tuple[bool, str]:
    """Run bcftools view -H and check exit code."""
    try:
        result = subprocess.run(
            ["bcftools", "view", "-H", vcf_path],
            capture_output=True,
            text=True,
            timeout=60,
        )
        if result.returncode != 0:
            return False, result.stderr.strip()[:200]
        return True, ""
    except FileNotFoundError:
        return True, "(bcftools not installed, skipped)"
    except subprocess.TimeoutExpired:
        return False, "bcftools timed out after 60s"


# Check 3: Assembly consistency (contig lengths)


def check_contig_lengths(vcf_path: str, assembly: str) -> tuple[bool, list[str]]:
    """Verify contig lengths in VCF header match expected assembly values."""
    expected = EXPECTED_CONTIG_LENGTHS.get(assembly, {})
    if not expected:
        return True, []

    errors: list[str] = []
    contig_re = re.compile(r"##contig=<ID=([^,>]+),length=(\d+)>")

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                break
            m = contig_re.search(line)
            if m:
                cid, clen = m.group(1), int(m.group(2))
                if cid in expected and clen != expected[cid]:
                    errors.append(
                        f"  contig {cid}: length={clen} "
                        f"expected={expected[cid]} for {assembly}"
                    )

    return len(errors) == 0, errors


# Check 4: Variant count + Check 6: VCF version


def count_variants_and_version(vcf_path: str) -> tuple[int, str]:
    """Count data lines and extract ##fileformat version."""
    count = 0
    version = "(unknown)"

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("##fileformat="):
                version = line.strip().split("=", 1)[1]
            elif not line.startswith("#"):
                count += 1

    return count, version


# Check 5: Real-world symmetry


def check_real_world_symmetry(
    test_dir: str,
) -> tuple[bool, list[str], dict[str, list[str]]]:
    """Verify grch37/real_world/ and grch38/real_world/ have matching sources.

    Returns (ok, errors, {assembly: [filenames]}).
    """
    errors: list[str] = []
    rw_files: dict[str, list[str]] = {}
    present: list[str] = []

    for assembly in ("grch37", "grch38"):
        rw_dir = os.path.join(test_dir, assembly, "real_world")
        if not os.path.isdir(rw_dir):
            rw_files[assembly] = []
            continue
        present.append(assembly)
        vcfs = sorted(
            f
            for f in os.listdir(rw_dir)
            if f.endswith(".vcf.gz") and not f.endswith(".tbi")
        )
        rw_files[assembly] = vcfs

    # Neither assembly has real-world VCFs staged. That is the state of a fresh
    # clone, so there is nothing to check and nothing to fail.
    if not present:
        return True, [], rw_files

    # Check that both have the same 3 sources (by stem)
    def stem_set(filenames: list[str]) -> set[str]:
        stems = set()
        for f in filenames:
            # Normalize: strip version numbers to get source stem
            # e.g. gnomad_sv_v2.1_chr21.vcf.gz -> gnomad_sv
            # e.g. clinvar_sv_chr21.vcf.gz -> clinvar_sv_chr21
            # e.g. 1kg_sv_chr21.vcf.gz -> 1kg_sv_chr21
            for expected_stem in EXPECTED_REAL_WORLD_STEMS:
                if f.startswith(expected_stem):
                    stems.add(expected_stem)
                    break
            else:
                stems.add(f.replace(".vcf.gz", ""))
        return stems

    grch37_stems = stem_set(rw_files.get("grch37", []))
    grch38_stems = stem_set(rw_files.get("grch38", []))

    if grch37_stems != grch38_stems:
        only_37 = grch37_stems - grch38_stems
        only_38 = grch38_stems - grch37_stems
        if only_37:
            errors.append(f"  Sources only in GRCh37: {sorted(only_37)}")
        if only_38:
            errors.append(f"  Sources only in GRCh38: {sorted(only_38)}")

    expected_count = len(EXPECTED_REAL_WORLD_STEMS)
    for assembly in ("grch37", "grch38"):
        actual_stems = stem_set(rw_files.get(assembly, []))
        missing = EXPECTED_REAL_WORLD_STEMS - actual_stems
        if missing:
            errors.append(
                f"  {assembly.upper()} missing expected sources: " f"{sorted(missing)}"
            )

    return len(errors) == 0, errors, rw_files


# Detect assembly from VCF header


def detect_assembly_from_header(vcf_path: str) -> Optional[str]:
    """Try to detect assembly from ##reference= line."""
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                break
            if line.startswith("##reference="):
                ref_val = line.strip().split("=", 1)[1]
                if "GRCh38" in ref_val or "grch38" in ref_val.lower():
                    return "GRCh38"
                if (
                    "GRCh37" in ref_val
                    or "grch37" in ref_val.lower()
                    or "hs37" in ref_val
                ):
                    return "GRCh37"
    return None


# Discover VCF files


def discover_vcfs(
    test_dir: str,
) -> list[tuple[str, str, str]]:
    """Find all VCF files. Returns [(path, assembly, category), ...]."""
    results = []
    for assembly_dir in ("grch37", "grch38"):
        assembly = assembly_dir.upper().replace("GRCH", "GRCh")
        asm_path = os.path.join(test_dir, assembly_dir)
        if not os.path.isdir(asm_path):
            continue

        # Synthetic files
        for f in sorted(os.listdir(asm_path)):
            if f.endswith(".vcf.gz") and not f.endswith(".tbi"):
                results.append((os.path.join(asm_path, f), assembly, "synthetic"))

        # Real-world files
        rw_path = os.path.join(asm_path, "real_world")
        if os.path.isdir(rw_path):
            for f in sorted(os.listdir(rw_path)):
                if f.endswith(".vcf.gz") and not f.endswith(".tbi"):
                    results.append((os.path.join(rw_path, f), assembly, "real_world"))

    return results


# Main validation


def validate_all(
    test_dir: str,
    fasta_grch37: Optional[str],
    fasta_grch38: Optional[str],
) -> tuple[bool, list[FileResult], dict]:
    """Run all validation checks. Returns (all_passed, results, summary)."""

    # Load FASTAs if provided
    fastas: dict[str, Optional[RefFasta]] = {"GRCh37": None, "GRCh38": None}
    if fasta_grch37:
        print("Loading GRCh37 reference...", file=sys.stderr)
        fastas["GRCh37"] = RefFasta(fasta_grch37, "21")
    if fasta_grch38:
        print("Loading GRCh38 reference...", file=sys.stderr)
        fastas["GRCh38"] = RefFasta(fasta_grch38, "21")

    # Discover files
    vcf_files = discover_vcfs(test_dir)
    if not vcf_files:
        print(f"ERROR: No VCF files found in {test_dir}", file=sys.stderr)
        return False, [], {}

    print(
        f"\nValidating {len(vcf_files)} VCF files in {test_dir}\n",
        file=sys.stderr,
    )

    results: list[FileResult] = []
    all_passed = True

    for vcf_path, assembly, category in vcf_files:
        rel_path = os.path.relpath(vcf_path, test_dir)
        print(f"  Checking {rel_path} ...", file=sys.stderr, end=" ")

        fr = FileResult(path=rel_path, assembly=assembly, category=category)

        # Check 4+6: variant count + version
        fr.variant_count, fr.vcf_version = count_variants_and_version(vcf_path)

        # Check 2: bcftools
        fr.bcftools_ok, fr.bcftools_error = check_bcftools(vcf_path)

        # Check 3: contig lengths
        fr.contig_ok, fr.contig_errors = check_contig_lengths(vcf_path, assembly)

        # Check 1: REF alleles (if FASTA available)
        fasta = fastas.get(assembly)
        if fasta:
            fr.ref_checked = True
            fr.ref_total, fr.ref_mismatches, fr.ref_errors = check_ref_alleles(
                vcf_path, fasta
            )

        status = "PASS" if fr.passed else "FAIL"
        if not fr.passed:
            all_passed = False
        print(
            f"{status}  (v={fr.vcf_version}, n={fr.variant_count:,})",
            file=sys.stderr,
        )

        results.append(fr)

    # Check 5: real-world symmetry
    sym_ok, sym_errors, rw_files = check_real_world_symmetry(test_dir)
    if not sym_ok:
        all_passed = False

    # Build summary
    summary = build_summary(results, sym_ok, sym_errors, rw_files, fastas)

    return all_passed, results, summary


# Summary and reporting


def build_summary(
    results: list[FileResult],
    sym_ok: bool,
    sym_errors: list[str],
    rw_files: dict[str, list[str]],
    fastas: dict[str, Optional[RefFasta]],
) -> dict:
    """Build JSON-serializable summary."""

    # Per-assembly counts
    assembly_counts: dict[str, dict[str, int]] = {}
    for fr in results:
        asm = fr.assembly
        if asm not in assembly_counts:
            assembly_counts[asm] = {
                "synthetic": 0,
                "real_world": 0,
                "total_variants": 0,
                "files": 0,
            }
        assembly_counts[asm][fr.category] += fr.variant_count
        assembly_counts[asm]["total_variants"] += fr.variant_count
        assembly_counts[asm]["files"] += 1

    total_variants = sum(ac["total_variants"] for ac in assembly_counts.values())
    total_files = sum(ac["files"] for ac in assembly_counts.values())

    # Failures
    failures = []
    for fr in results:
        if not fr.passed:
            failure = {"file": fr.path, "assembly": fr.assembly, "issues": []}
            if not fr.bcftools_ok:
                failure["issues"].append(f"bcftools: {fr.bcftools_error}")
            if not fr.contig_ok:
                failure["issues"].extend(fr.contig_errors)
            if fr.ref_checked and fr.ref_mismatches > 0:
                failure["issues"].append(
                    f"REF mismatches: {fr.ref_mismatches}/{fr.ref_total}"
                )
                failure["issues"].extend(fr.ref_errors)
            failures.append(failure)

    # File details
    file_details = []
    for fr in results:
        detail = {
            "file": fr.path,
            "assembly": fr.assembly,
            "category": fr.category,
            "vcf_version": fr.vcf_version,
            "variants": fr.variant_count,
            "bcftools_ok": fr.bcftools_ok,
            "contig_ok": fr.contig_ok,
            "ref_checked": fr.ref_checked,
            "ref_total": fr.ref_total,
            "ref_mismatches": fr.ref_mismatches,
            "passed": fr.passed,
        }
        file_details.append(detail)

    return {
        "overall_passed": len(failures) == 0 and sym_ok,
        "total_files": total_files,
        "total_variants": total_variants,
        "assembly_counts": assembly_counts,
        "ref_validation": {
            "grch37_fasta": fastas.get("GRCh37") is not None,
            "grch38_fasta": fastas.get("GRCh38") is not None,
        },
        "real_world_symmetry": {
            "passed": sym_ok,
            "errors": sym_errors,
            "files": rw_files,
        },
        "failures": failures,
        "files": file_details,
    }


def print_human_report(
    results: list[FileResult],
    summary: dict,
    all_passed: bool,
) -> None:
    """Print human-readable report to stderr."""

    print("\n" + "=" * 78, file=sys.stderr)
    print("VCF VALIDATION REPORT", file=sys.stderr)
    print("=" * 78, file=sys.stderr)

    # Per-file table
    print(
        f"\n{'File':<45} {'Ver':<8} {'Variants':>9} "
        f"{'bcf':>4} {'ctg':>4} {'REF':>10} {'Status':>7}",
        file=sys.stderr,
    )
    print("-" * 78 + "  -------", file=sys.stderr)

    for fr in results:
        bcf_sym = "ok" if fr.bcftools_ok else "FAIL"
        ctg_sym = "ok" if fr.contig_ok else "FAIL"
        if fr.ref_checked:
            if fr.ref_mismatches == 0:
                ref_sym = f"{fr.ref_total:,} ok"
            else:
                ref_sym = f"{fr.ref_mismatches}/{fr.ref_total}"
        else:
            ref_sym = "skip"
        status = "PASS" if fr.passed else "FAIL"

        print(
            f"{fr.path:<45} {fr.vcf_version:<8} {fr.variant_count:>9,} "
            f"{bcf_sym:>4} {ctg_sym:>4} {ref_sym:>10} {status:>7}",
            file=sys.stderr,
        )

    # Assembly totals
    print("\nVariant counts by assembly:", file=sys.stderr)
    for asm, counts in sorted(summary["assembly_counts"].items()):
        print(
            f"  {asm}: {counts['total_variants']:,} total "
            f"({counts['synthetic']:,} synthetic, "
            f"{counts['real_world']:,} real-world) "
            f"across {counts['files']} files",
            file=sys.stderr,
        )

    cross_total = summary["total_variants"]
    print(
        f"  Cross-assembly total: {cross_total:,} variants "
        f"across {summary['total_files']} files",
        file=sys.stderr,
    )

    # REF validation status
    ref_info = summary["ref_validation"]
    print("\nREF allele validation:", file=sys.stderr)
    for asm in ("grch37", "grch38"):
        key = f"{asm}_fasta"
        if ref_info[key]:
            asm_results = [
                r for r in results if r.assembly == asm.upper().replace("GRCH", "GRCh")
            ]
            total_checked = sum(r.ref_total for r in asm_results)
            total_mismatches = sum(r.ref_mismatches for r in asm_results)
            print(
                f"  {asm.upper()}: {total_checked:,} REF alleles checked, "
                f"{total_mismatches:,} mismatches",
                file=sys.stderr,
            )
        else:
            print(
                f"  {asm.upper()}: SKIPPED (no FASTA provided)",
                file=sys.stderr,
            )

    # Symmetry check
    sym = summary["real_world_symmetry"]
    print(
        f"\nReal-world symmetry: {'PASS' if sym['passed'] else 'FAIL'}",
        file=sys.stderr,
    )
    if not sym["passed"]:
        for e in sym["errors"]:
            print(f"  {e}", file=sys.stderr)

    # Failures
    failures = summary["failures"]
    if failures:
        print(f"\nFAILURES ({len(failures)}):", file=sys.stderr)
        for fail in failures:
            print(f"\n  {fail['file']} ({fail['assembly']}):", file=sys.stderr)
            for issue in fail["issues"]:
                print(f"    {issue}", file=sys.stderr)

    # Final verdict
    print("\n" + "=" * 78, file=sys.stderr)
    if all_passed:
        print("RESULT: ALL CHECKS PASSED", file=sys.stderr)
    else:
        print("RESULT: VALIDATION FAILED", file=sys.stderr)
    print("=" * 78 + "\n", file=sys.stderr)


# CLI


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate all test VCFs in the vep-rs SV validation suite."
    )
    parser.add_argument(
        "--fasta-grch37",
        help="Path to GRCh37 reference FASTA (enables REF validation)",
    )
    parser.add_argument(
        "--fasta-grch38",
        help="Path to GRCh38 reference FASTA (enables REF validation)",
    )
    parser.add_argument(
        "--test-dir",
        default="tests/sv_validation",
        help="Root directory of test VCFs (default: tests/sv_validation)",
    )
    args = parser.parse_args()

    # Validate FASTA paths if provided
    for label, path in [
        ("--fasta-grch37", args.fasta_grch37),
        ("--fasta-grch38", args.fasta_grch38),
    ]:
        if path and not os.path.isfile(path):
            print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
            return 1

    # Validate test dir
    if not os.path.isdir(args.test_dir):
        print(
            f"ERROR: test directory not found: {args.test_dir}",
            file=sys.stderr,
        )
        return 1

    # Run validation
    all_passed, results, summary = validate_all(
        test_dir=args.test_dir,
        fasta_grch37=args.fasta_grch37,
        fasta_grch38=args.fasta_grch38,
    )

    # Human-readable report to stderr
    print_human_report(results, summary, all_passed)

    # JSON report to stdout
    print(json.dumps(summary, indent=2))

    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
