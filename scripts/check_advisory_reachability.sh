#!/usr/bin/env bash
#
# Fail when a known third-party advisory is reachable from the production binary.
#
# `cargo audit` answers "does this tree depend on anything with an advisory?".
# It does not answer "is that thing in the binary users run?", and the two are
# very different: an advisory can be reachable only from `vep-cache-builder`, the
# optional offline cache-build utility that talks to Ensembl's servers and
# processes nothing a user supplies, and never from the annotator itself.
#
# The release condition is that no third-party vulnerability is reachable from
# the production binary via untrusted input. Argued in prose and enforced by
# nothing, that condition breaks silently: an advisory confined to a dev-only or
# cache-builder-only crate sits in the tree behind a paragraph explaining why
# its defect cannot be triggered, and a refactor that pulls that crate into
# `vep-cli` invalidates the paragraph without failing anything.
#
# Definitions:
#   production binary = `vep-cli` (the `vep` binary, the only one that consumes
#                       user-supplied variant files). `vep-cache-builder` and
#                       `vep-cache-converter` are the utilities the condition's
#                       own reachability argument excludes, so advisories
#                       confined to them warn instead of failing.
#
# Two flags carry the correctness of the reachability set:
#
#   --target all   `openssl` does not appear in the host-target graph AT ALL,
#                  because `native-tls` uses Security.framework on macOS and
#                  OpenSSL on Linux. Without this flag the gate reports a clean
#                  result for every openssl advisory on the released Linux
#                  artifact.
#
#   --edges normal Dev-dependencies (the benchmark harness among them) are not
#                  in the production binary; `--edges all` would flag advisories
#                  that ship to nobody.
#
# Allowlist (see the file named below): entries suppress a REACHABLE advisory and
# each one must carry an expiry date. An expired entry that is still suppressing
# a live finding FAILS rather than lapsing quietly, and an entry that suppresses
# nothing warns so it gets deleted. A permanent silence is the failure mode this
# design exists to prevent.
#
# Usage: scripts/check_advisory_reachability.sh
# Exit: 0 pass (possibly with warnings), 1 a reachable advisory, 2 cannot check.

set -euo pipefail

PROD_PACKAGE="vep-cli"
ALLOWLIST="${ALLOWLIST:-scripts/advisory-reachability-allow.txt}"

if ! command -v cargo >/dev/null 2>&1; then
    echo "ERROR: [check_advisory_reachability] cargo not on PATH" >&2
    exit 2
fi

# An absent scanner exits 2 ("could not check"), never 0. A gate that silently
# passes when its tool is missing is worse than no gate, because CI stays green.
if ! cargo audit --version >/dev/null 2>&1; then
    echo "ERROR: [check_advisory_reachability] cargo-audit not installed" >&2
    echo "Install: cargo install cargo-audit --locked" >&2
    exit 2
fi

if ! cargo tree -p "$PROD_PACKAGE" --help >/dev/null 2>&1; then
    echo "ERROR: [check_advisory_reachability] cargo tree unavailable" >&2
    exit 2
fi

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT

# The set of crates the production binary actually contains. See the flag notes
# above before touching this command.
if ! cargo tree -p "$PROD_PACKAGE" --target all --edges normal \
    --prefix none --no-dedupe >"$workdir/tree.txt" 2>"$workdir/tree.err"; then
    echo "ERROR: [check_advisory_reachability] cargo tree failed for $PROD_PACKAGE:" >&2
    cat "$workdir/tree.err" >&2
    exit 2
fi
sed 's/ v[0-9].*//' "$workdir/tree.txt" | sed '/^$/d' | sort -u >"$workdir/prod.txt"

reachable_count=$(grep -c . "$workdir/prod.txt") || reachable_count=0
if [[ $reachable_count -lt 10 ]]; then
    # A near-empty reachability set means the tree command silently changed
    # shape, not that the binary depends on nothing. Refuse rather than pass
    # everything.
    echo "ERROR: [check_advisory_reachability] implausible reachability set ($reachable_count crates)" >&2
    exit 2
fi

# `cargo audit` exits non-zero when it finds anything, which is the normal case
# here, so its status is captured rather than allowed to kill the script.
audit_status=0
cargo audit --json >"$workdir/audit.json" 2>"$workdir/audit.err" || audit_status=$?
if [[ ! -s $workdir/audit.json ]]; then
    echo "ERROR: [check_advisory_reachability] cargo audit produced no JSON (exit $audit_status):" >&2
    cat "$workdir/audit.err" >&2
    exit 2
fi

if [[ -f $ALLOWLIST ]]; then
    cp "$ALLOWLIST" "$workdir/allow.txt"
else
    : >"$workdir/allow.txt"
fi

# Python does the JSON walk and the date arithmetic; both are miserable in the
# bash 3.2 that ships with macOS, where this gate is also run by hand.
python3 - "$workdir/audit.json" "$workdir/prod.txt" "$workdir/allow.txt" <<'PYEOF'
import datetime, json, sys

audit_path, prod_path, allow_path = sys.argv[1], sys.argv[2], sys.argv[3]

with open(audit_path) as fh:
    audit = json.load(fh)

prod = {line.strip() for line in open(prod_path) if line.strip()}

allow = {}
for lineno, raw in enumerate(open(allow_path), 1):
    line = raw.split("#", 1)[0].strip()
    if not line:
        continue
    parts = [p.strip() for p in line.split("|")]
    if len(parts) != 3 or not all(parts):
        print("ERROR: [check_advisory_reachability] malformed allowlist line %d: %s"
              % (lineno, raw.rstrip()), file=sys.stderr)
        print("Expected: RUSTSEC-ID | YYYY-MM-DD | reason", file=sys.stderr)
        sys.exit(2)
    rid, expiry, reason = parts
    try:
        expiry_date = datetime.date.fromisoformat(expiry)
    except ValueError:
        print("ERROR: [check_advisory_reachability] allowlist line %d has a bad date: %s"
              % (lineno, expiry), file=sys.stderr)
        sys.exit(2)
    allow[rid] = (expiry_date, reason)

# Vulnerabilities and warnings are shaped alike but live under different keys,
# and warnings are a dict of kind -> list. Both are walked: cargo audit files
# an unsoundness advisory as a "warning", and one of those can be reachable
# from vep-cli just as a vulnerability can.
findings = []
for item in (audit.get("vulnerabilities") or {}).get("list") or []:
    findings.append((item["advisory"]["id"], item["package"]["name"], "vulnerability"))
warnings = audit.get("warnings") or {}
for kind, items in warnings.items():
    for item in items or []:
        findings.append((item["advisory"]["id"], item["package"]["name"], kind))

today = datetime.date.today()
blocking, suppressed, confined = [], [], []

for rid, crate, kind in sorted(set(findings)):
    if crate not in prod:
        confined.append((rid, crate, kind))
        continue
    if rid in allow:
        expiry_date, reason = allow[rid]
        if expiry_date < today:
            blocking.append((rid, crate, kind,
                             "allowlist entry EXPIRED %s" % expiry_date.isoformat()))
        else:
            suppressed.append((rid, crate, kind, expiry_date, reason))
    else:
        blocking.append((rid, crate, kind, "not allowlisted"))

used = {s[0] for s in suppressed} | {b[0] for b in blocking if "EXPIRED" in b[3]}
stale = sorted(set(allow) - used)

print("[check_advisory_reachability] %s reachable crates; %d advisories"
      % (len(prod), len(set(findings))))

if confined:
    print()
    print("Not reachable from the production binary (%d), reported for the record:"
          % len(confined))
    for rid, crate, kind in confined:
        print("  ok    %-20s %-22s %s" % (rid, crate, kind))

if suppressed:
    print()
    print("Reachable but allowlisted (%d):" % len(suppressed))
    for rid, crate, kind, expiry_date, reason in suppressed:
        print("  allow %-20s %-22s expires %s -- %s"
              % (rid, crate, expiry_date.isoformat(), reason))

if stale:
    print()
    print("WARNING: allowlist entries matching nothing (%d) -- delete them:" % len(stale))
    for rid in stale:
        print("  stale %s" % rid)

if blocking:
    print(file=sys.stderr)
    print("ERROR: [check_advisory_reachability] %d advisory/advisories are reachable "
          "from %s:" % (len(blocking), "vep-cli"), file=sys.stderr)
    for rid, crate, kind, why in blocking:
        print("  FAIL  %-20s %-22s %-14s %s" % (rid, crate, kind, why), file=sys.stderr)
    print(file=sys.stderr)
    print("The release requirement is that no third-party vulnerability is reachable "
          "from the", file=sys.stderr)
    print("production binary via untrusted input. Fix by upgrading the crate,",
          file=sys.stderr)
    print("by removing the dependency edge into vep-cli, or -- if the defect genuinely "
          "cannot be", file=sys.stderr)
    print("triggered -- by adding an EXPIRING entry to the allowlist with the reason:",
          file=sys.stderr)
    print("  RUSTSEC-XXXX-NNNN | YYYY-MM-DD | why this cannot be triggered",
          file=sys.stderr)
    sys.exit(1)

print()
print("PASS: no advisory is reachable from the production binary.")
PYEOF
