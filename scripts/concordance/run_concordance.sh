#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=scripts/data/plugin_data_layout.sh
source "${SCRIPT_DIR}/../data/plugin_data_layout.sh"
VEP_RS_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# Third-party checkouts live outside this repository. The default below is a
# sibling of the repo root purely as a convention; set the matching environment
# variable or flag to put it anywhere.
WORKSPACE_DIR="$(cd "${VEP_RS_DIR}/.." && pwd)"

# Benchmark VCFs, under the same ./.vep/ root the other scripts use. Not
# committed: seed them with scripts/data/download_real_world_vcfs.sh, which
# fetches from public sources, or point --benchmark-dir at your own directory.
BENCHMARK_DIR="${VEP_BENCHMARK_DIR:-${VEP_RS_DIR}/.vep/benchmark}"
# Upstream Perl VEP checkout, needed only for --perl-vep-mode local.
ENSEMBL_VEP_DIR="${VEP_ENSEMBL_VEP_DIR:-${WORKSPACE_DIR}/ensembl-vep}"
WORK_DIR="${VEP_RS_DIR}/tmp/concordance"

MODE="smoke"
SMOKE_VARIANTS=5000
SPECIES="homo_sapiens"
ASSEMBLY="GRCh37"
PERL_VEP_MODE="docker"
PERL_VEP_DOCKER_IMAGE="ensemblorg/ensembl-vep:release_115.2"
FIXTURES_DIR=""
RUST_RELEASE=0

PERL_CACHE_DIR=""
JSON_CACHE_DIR=""
FASTA=""
NO_FASTA=0
CANONICAL_CONTIGS=1
NO_CANONICAL_CONTIGS=0

PERL_OUTPUT_CACHE_DIR=""
PERL_OUTPUT_CACHE_ENABLED=1
PERL_OUTPUT_CACHE_REFRESH=0

PERF_ENABLED=1
PERF_DIR=""
PERF_SAMPLE_SECONDS=2

# Plugin testing
PLUGINS=""
PLUGIN_DATA_DIR=""
PLUGIN_TIER="quick"
# Per-engine floor of annotated records below which a plugin FAILS. Zero
# annotations is the signature of a broken plugin (contig-name mismatch,
# unparsed header, wrong position column), not of perfect agreement.
PLUGIN_MIN_ANNOTATED=1
# Upstream Perl VEP_plugins checkout, mounted into the Perl container when it
# exists. Optional: absent means the Perl side runs without Perl-side plugins.
VEP_PLUGINS_DIR="${VEP_PLUGINS_DIR:-${WORKSPACE_DIR}/VEP_plugins}"

# Final process exit status. A verdict that must fail the run sets this, and the
# script exits on it after printing its report paths, so a failure is both visible
# and machine-detectable.
RUN_STATUS=0

# Provenance labels (optional; useful when aggregating runs across machines)
RUN_ID=""
HOST_ID=""

usage() {
    cat <<'EOF'
Run semantic concordance evaluation between ensembl-vep and vep-rs.

Usage:
  scripts/concordance/run_concordance.sh --perl-cache-dir <path> --json-cache-dir <path> \
      --fasta <path> [options]

Required:
  --perl-cache-dir <path>   Perl VEP cache root, the directory above homo_sapiens/ (passed to VEP as --dir_cache)
  --json-cache-dir <path>   Converted JSON cache for vep-rs (must contain info.json)
  --fasta <path>            Reference FASTA, or --no-fasta to run every engine without one

Optional:
  --perl-vep-mode <local|docker|fixtures>  Perl VEP execution mode (default: docker)
                                    fixtures: use pre-existing Perl outputs from --fixtures-dir
                                              (skips Docker/Perl entirely)
  --perl-vep-docker-image <image>  Docker image for perl VEP (default: ensemblorg/ensembl-vep:release_115.2)
  --fixtures-dir <path>            Cached Perl VEP outputs (default:
                                   tests/concordance_fixtures/sv_perl, which is
                                   not committed; produce it with
                                   --perl-vep-mode docker first)
  --rust-release            Run vep-rs using a release build (cargo --release)
  --perl-output-cache-dir <path>   Cache dir for perl VEP outputs (default: vep-rs/tmp/perl_output_cache)
  --no-perl-output-cache   Disable perl output caching (always rerun perl VEP)
  --refresh-perl-output-cache  Ignore cached perl outputs; rerun and update cache
  --perf-dir <path>         Directory for performance reports (default: <work-dir>/reports/perf)
  --no-perf                 Disable performance monitoring
  --perf-sample-seconds <int>  Docker stats sampling interval (default: 2)
  --benchmark-dir <path>    Directory of .vcf / .vcf.gz inputs. Not committed;
                            seed it with scripts/data/download_real_world_vcfs.sh.
                            Default: ./.vep/benchmark/ under the repo root, or
                            $VEP_BENCHMARK_DIR.
  --ensembl-vep-dir <path>  Upstream Perl VEP checkout, for --perl-vep-mode local
                            only. Default: an ensembl-vep/ sibling of the repo, or
                            $VEP_ENSEMBL_VEP_DIR.
  --work-dir <path>         Working output directory (default: vep-rs/tmp/concordance)
  --species <name>          Species (default: homo_sapiens)
  --assembly <name>         Assembly (default: GRCh37)
  --fasta <path>            Reference FASTA. Enables sequence-derived UTR predicates on Perl +
                            vep-rs + fastVEP, plus HGVS 3' normalization. Both engines MUST
                            receive --fasta or both MUST be invoked without it; an asymmetric
                            run reports sequence-dependent differences as false discordants.
  --no-fasta                Run every engine without a FASTA. Without this flag, an empty
                            --fasta is a hard-fail.
  --no-canonical-contigs    Opt out of the canonical-contigs input filter and annotate every
                            contig the input carries. DEFAULT: input VCFs are filtered to the
                            canonical contigs 1-22,X,Y,MT (Ensembl style, for both assemblies;
                            a UCSC chr prefix is stripped first) before being fed to the
                            engines, so that asymmetric alt-contig handling between engines
                            cannot inflate F1 divergences.
  --mode <smoke|full>       Run mode (default: smoke)
  --smoke-variants <int>    Variants per file in smoke mode (default: 5000)
  --plugins <list>          Comma-separated plugins to test (e.g. CADD,REVEL,gnomADc)
  --plugin-data-dir <path>  Plugin data file directory (default: ~/.vep/plugin_data/<assembly>/<tier>)
  --plugin-tier <quick|full>  Plugin data tier (default: quick)
  --run-id <id>             Identifier for this measurement campaign; recorded
                            into provenance.json and every wall-time CSV row's
                            notes column. Optional; useful for grouping runs.
  --host-id <id>            Identifier for the machine running this invocation;
                            recorded into provenance.json and every wall-time CSV
                            row. Optional (e.g. a hostname, when aggregating runs
                            from more than one machine).
  --help                    Print this help

Outputs:
  <work-dir>/inputs     prepared VCFs
  <work-dir>/perl       perl VEP outputs
  <work-dir>/rust       rust VEP outputs
  <work-dir>/reports    concordance reports
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    --perl-cache-dir)
        PERL_CACHE_DIR="$2"
        shift 2
        ;;
    --json-cache-dir)
        JSON_CACHE_DIR="$2"
        shift 2
        ;;
    --perl-vep-mode)
        PERL_VEP_MODE="$2"
        shift 2
        ;;
    --perl-vep-docker-image)
        PERL_VEP_DOCKER_IMAGE="$2"
        shift 2
        ;;
    --fixtures-dir)
        FIXTURES_DIR="$2"
        shift 2
        ;;
    --rust-release)
        RUST_RELEASE=1
        shift 1
        ;;
    --perl-output-cache-dir)
        PERL_OUTPUT_CACHE_DIR="$2"
        shift 2
        ;;
    --no-perl-output-cache)
        PERL_OUTPUT_CACHE_ENABLED=0
        shift 1
        ;;
    --refresh-perl-output-cache)
        PERL_OUTPUT_CACHE_REFRESH=1
        shift 1
        ;;
    --perf-dir)
        PERF_DIR="$2"
        shift 2
        ;;
    --no-perf)
        PERF_ENABLED=0
        shift 1
        ;;
    --perf-sample-seconds)
        PERF_SAMPLE_SECONDS="$2"
        shift 2
        ;;
    --benchmark-dir)
        BENCHMARK_DIR="$2"
        shift 2
        ;;
    --ensembl-vep-dir)
        ENSEMBL_VEP_DIR="$2"
        shift 2
        ;;
    --work-dir)
        WORK_DIR="$2"
        shift 2
        ;;
    --species)
        SPECIES="$2"
        shift 2
        ;;
    --assembly)
        ASSEMBLY="$2"
        shift 2
        ;;
    --fasta)
        FASTA="$2"
        shift 2
        ;;
    --no-fasta)
        NO_FASTA=1
        shift 1
        ;;
    --no-canonical-contigs)
        NO_CANONICAL_CONTIGS=1
        CANONICAL_CONTIGS=0
        shift 1
        ;;
    --mode)
        MODE="$2"
        shift 2
        ;;
    --smoke-variants)
        SMOKE_VARIANTS="$2"
        shift 2
        ;;
    --plugins)
        PLUGINS="$2"
        shift 2
        ;;
    --plugin-data-dir)
        PLUGIN_DATA_DIR="$2"
        shift 2
        ;;
    --plugin-tier)
        PLUGIN_TIER="$2"
        shift 2
        ;;
    --run-id)
        RUN_ID="$2"
        shift 2
        ;;
    --host-id)
        HOST_ID="$2"
        shift 2
        ;;
    --help | -h)
        usage
        exit 0
        ;;
    *)
        echo "Unknown argument: $1" >&2
        usage >&2
        exit 1
        ;;
    esac
done

if [[ -z "${PERL_CACHE_DIR}" ]]; then
    echo "ERROR: --perl-cache-dir is required" >&2
    exit 1
fi

if [[ "${MODE}" != "smoke" && "${MODE}" != "full" ]]; then
    echo "ERROR: --mode must be smoke or full" >&2
    exit 1
fi

if [[ "${PERL_VEP_MODE}" != "local" && "${PERL_VEP_MODE}" != "docker" && "${PERL_VEP_MODE}" != "fixtures" ]]; then
    echo "ERROR: --perl-vep-mode must be local, docker, or fixtures" >&2
    exit 1
fi

# FASTA-required gate: asymmetric --fasta between engines reports sequence-dependent
# differences as false discordants. Either ALL three engines (Perl, vep-rs, fastVEP)
# get --fasta, or none do; the harness enforces "all". --no-fasta is the explicit
# opt-out that runs every engine without one.
if [[ -z "${FASTA}" && "${NO_FASTA}" -ne 1 ]]; then
    echo "ERROR: --fasta is required (or pass --no-fasta to opt out explicitly)" >&2
    echo "  Reason: sequence-derived UTR predicates (_ins_del_start_altered," >&2
    echo "    _ins_del_stop_altered) and HGVS 3' normalization require the reference" >&2
    echo "    FASTA. Asymmetric FASTA between Perl and vep-rs reports sequence-dependent" >&2
    echo "    differences as false discordants." >&2
    echo "  Fix: pass --fasta /path/to/Homo_sapiens.GRCh3{7,8}.dna.primary_assembly.fa" >&2
    echo "    (with .fai index alongside)." >&2
    echo "  Opt-out: pass --no-fasta to run every engine without a FASTA. The output" >&2
    echo "    is tagged with the no-fasta caveat in provenance.json." >&2
    exit 1
fi
if [[ -n "${FASTA}" && ! -f "${FASTA}" ]]; then
    echo "ERROR: --fasta path does not exist: ${FASTA}" >&2
    exit 1
fi
if [[ -n "${FASTA}" && ! -f "${FASTA}.fai" ]]; then
    echo "WARNING: --fasta index missing: ${FASTA}.fai (engines will build it on first read)" >&2
fi

# Canonical-contigs gate: real-world clinical/research
# VCFs annotate against 1-22,X,Y,MT only. Unplaced contigs (NT_*, NW_*, decoys) are
# sequencing artifacts that don't belong in F1. Each engine handles alt contigs differently
# (Perl includes some, vep-rs includes everything in cache, fastVEP drops them entirely),
# producing F1 divergence on multi-chromosome cells purely from contig-handling
# asymmetry, NOT from real consequence-logic differences.
#
# Default: filter input VCFs to canonical contigs before feeding to engines. The filtered
# VCF SHAs are captured in provenance.json under canonical_contigs_filter_active +
# filtered_vcf_sha256s. --no-canonical-contigs is the explicit opt-out for comparison
# against a "with alt contigs" baseline.
if [[ "${CANONICAL_CONTIGS}" -eq 1 ]]; then
    for cmd in bcftools tabix bgzip; do
        if ! command -v "${cmd}" >/dev/null 2>&1; then
            echo "ERROR: canonical-contigs filter requires ${cmd} (install htslib + bcftools)" >&2
            echo "  Or pass --no-canonical-contigs to opt out and annotate every contig" >&2
            exit 1
        fi
    done
fi

# Canonical contig set: universal
# Ensembl-style for both GRCh37 and GRCh38. The assembly difference is in the
# genome-sequence content, NOT contig naming. UCSC-style chr-prefix inputs are
# auto-stripped at filter time so the output is always Ensembl-style. Must match
# scripts/concordance/prepare_benchmark_vcfs.py CANONICAL_CONTIGS frozenset and
# scripts/concordance/compare_vep_outputs.py CANONICAL_CONTIGS frozenset.
#
# The variable is CANONICAL_CONTIGS_LIST, not CANONICAL_CONTIGS: the latter is the
# 0|1 flag among the defaults above that tracks whether the filter is active.
CANONICAL_CONTIGS_LIST="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT,M"

if [[ "${PERL_VEP_MODE}" == "fixtures" ]]; then
    if [[ -z "${FIXTURES_DIR}" ]]; then
        FIXTURES_DIR="${VEP_FIXTURES_DIR:-${VEP_RS_DIR}/tests/concordance_fixtures/sv_perl}"
    fi
    if [[ ! -d "${FIXTURES_DIR}" ]]; then
        echo "ERROR: [run_concordance] no Perl VEP fixtures at ${FIXTURES_DIR}" >&2
        echo "  Fixtures are cached Perl VEP outputs. They are NOT committed: they are" >&2
        echo "  large and specific to one cache release. Produce them with" >&2
        echo "    --perl-vep-mode docker" >&2
        echo "  which runs Perl VEP itself, then point --fixtures-dir at the saved" >&2
        echo "  outputs on later runs. See scripts/README.md." >&2
        exit 1
    fi
fi

for cmd in cargo python3; do
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        echo "ERROR: missing required tool: ${cmd}" >&2
        exit 1
    fi
done

if [[ "${PERL_VEP_MODE}" == "local" ]]; then
    if ! command -v perl >/dev/null 2>&1; then
        echo "ERROR: missing required tool for local mode: perl" >&2
        exit 1
    fi
fi

if [[ "${PERL_VEP_MODE}" == "docker" ]]; then
    if ! command -v docker >/dev/null 2>&1; then
        echo "ERROR: missing required tool for docker mode: docker" >&2
        exit 1
    fi
    if ! docker image inspect "${PERL_VEP_DOCKER_IMAGE}" >/dev/null 2>&1; then
        echo "pulling docker image: ${PERL_VEP_DOCKER_IMAGE}"
        docker pull "${PERL_VEP_DOCKER_IMAGE}"
    fi
fi

if [[ ! -d "${BENCHMARK_DIR}" ]]; then
    echo "ERROR: [run_concordance] no benchmark directory at ${BENCHMARK_DIR}" >&2
    echo "  Benchmark VCFs are not committed. Seed them from public sources with" >&2
    echo "    scripts/data/download_real_world_vcfs.sh" >&2
    echo "  or point --benchmark-dir at a directory of .vcf / .vcf.gz inputs." >&2
    exit 1
fi

shopt -s nullglob
benchmark_inputs=("${BENCHMARK_DIR}"/*.vcf "${BENCHMARK_DIR}"/*.vcf.gz)
shopt -u nullglob
if ((${#benchmark_inputs[@]} == 0)); then
    echo "ERROR: [run_concordance] no .vcf or .vcf.gz files in ${BENCHMARK_DIR}" >&2
    echo "  Seed them with scripts/data/download_real_world_vcfs.sh." >&2
    exit 1
fi

if [[ "${PERL_VEP_MODE}" == "local" ]]; then
    if [[ ! -f "${ENSEMBL_VEP_DIR}/vep" ]]; then
        echo "ERROR: [run_concordance] no ensembl-vep executable at ${ENSEMBL_VEP_DIR}/vep" >&2
        echo "  --perl-vep-mode local needs an upstream Perl VEP checkout. Clone it from" >&2
        echo "    https://github.com/Ensembl/ensembl-vep" >&2
        echo "  and pass --ensembl-vep-dir, or use --perl-vep-mode docker instead." >&2
        exit 1
    fi
fi

if [[ ! -d "${PERL_CACHE_DIR}" ]]; then
    echo "ERROR: perl cache directory not found: ${PERL_CACHE_DIR}" >&2
    exit 1
fi

mkdir -p "${WORK_DIR}"/{inputs,perl,rust,reports}

if [[ -z "${PERL_OUTPUT_CACHE_DIR}" ]]; then
    PERL_OUTPUT_CACHE_DIR="${VEP_RS_DIR}/tmp/perl_output_cache"
fi
if [[ "${PERL_OUTPUT_CACHE_ENABLED}" -eq 1 ]]; then
    mkdir -p "${PERL_OUTPUT_CACHE_DIR}"
    PERL_OUTPUT_CACHE_DIR="$(cd "${PERL_OUTPUT_CACHE_DIR}" && pwd)"
fi

if [[ "${PERF_ENABLED}" -eq 1 ]]; then
    if [[ -z "${PERF_DIR}" ]]; then
        PERF_DIR="${WORK_DIR}/reports/perf"
    fi
    mkdir -p "${PERF_DIR}"
    PERF_DIR="$(cd "${PERF_DIR}" && pwd)"
    # Ensure sampling interval is sane.
    if ! [[ "${PERF_SAMPLE_SECONDS}" =~ ^[0-9]+$ ]] || [[ "${PERF_SAMPLE_SECONDS}" -lt 1 ]]; then
        echo "ERROR: --perf-sample-seconds must be an integer >= 1" >&2
        exit 1
    fi
fi

sha256_file() {
    python3 - "$1" <<'PY'
import hashlib, sys
path = sys.argv[1]
h = hashlib.sha256()
with open(path, "rb") as f:
    for chunk in iter(lambda: f.read(1024 * 1024), b""):
        h.update(chunk)
print(h.hexdigest())
PY
}

sha256_string() {
    python3 - "$1" <<'PY'
import hashlib, sys
print(hashlib.sha256(sys.argv[1].encode("utf-8")).hexdigest())
PY
}

TIME_MODE="none"
if [[ "${PERF_ENABLED}" -eq 1 ]]; then
    if [[ -x "/usr/bin/time" ]] && /usr/bin/time -l true >/dev/null 2>&1; then
        TIME_MODE="bsd"
    elif [[ -x "/usr/bin/time" ]] && /usr/bin/time -v true >/dev/null 2>&1; then
        TIME_MODE="gnu"
    elif [[ -x "/usr/bin/time" ]]; then
        TIME_MODE="posix"
    fi
fi

perf_write_host_info() {
    if [[ "${PERF_ENABLED}" -ne 1 ]]; then
        return
    fi
    python3 - <<PY
import json, platform, subprocess
from datetime import datetime, timezone
from pathlib import Path

def run(cmd):
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
        return {"cmd": cmd, "returncode": p.returncode, "output": p.stdout}
    except Exception as e:
        return {"cmd": cmd, "error": str(e)}

info = {
  "generated_at_utc": datetime.now(timezone.utc).isoformat(),
  "platform": platform.platform(),
  "python": platform.python_version(),
  "uname": run(["uname", "-a"]),
}
if "${PERL_VEP_MODE}" == "docker":
  info["docker_version"] = run(["docker", "version"])
  info["docker_info"] = run(["docker", "info"])

out = Path("${PERF_DIR}") / "host.json"
out.write_text(json.dumps(info, indent=2), encoding="utf-8")
print(f"wrote perf host info: {out}")
PY
}

perf_parse_time_to_json() {
    if [[ "${PERF_ENABLED}" -ne 1 ]]; then
        return
    fi
    local label="$1"
    local time_raw="$2"
    local wall_seconds="$3"
    local exit_code="$4"
    local out_json="${PERF_DIR}/${label}.time.json"
    python3 - "${time_raw}" "${out_json}" "${TIME_MODE}" "${wall_seconds}" "${exit_code}" <<'PY'
import json, re, sys
from pathlib import Path

time_raw, out_json, mode, wall_s, rc = sys.argv[1:]
wall_s = int(wall_s)
rc = int(rc)
text = Path(time_raw).read_text(encoding="utf-8", errors="replace") if Path(time_raw).exists() else ""

metrics = {
  "time_mode": mode,
  "wall_seconds": wall_s,
  "exit_code": rc,
}

if mode == "bsd":
  m = re.search(r"^\s*([0-9.]+)\s+real\s+([0-9.]+)\s+user\s+([0-9.]+)\s+sys\s*$", text, re.M)
  if m:
    metrics["time_real_seconds"] = float(m.group(1))
    metrics["time_user_seconds"] = float(m.group(2))
    metrics["time_sys_seconds"] = float(m.group(3))
  m = re.search(r"^\s*([0-9.]+)%\s+cpu\s*$", text, re.M)
  if m:
    metrics["cpu_percent"] = float(m.group(1))
  m = re.search(r"^\s*(\d+)\s+maximum resident set size\s*$", text, re.M)
  if m:
    # BSD time (-l) reports bytes on macOS.
    metrics["max_rss_bytes"] = int(m.group(1))
elif mode == "gnu":
  def kv(key):
    mm = re.search(rf"^{re.escape(key)}\s*:\s*(.+)$", text, re.M)
    return mm.group(1).strip() if mm else None
  if kv("User time (seconds)") is not None:
    metrics["time_user_seconds"] = float(kv("User time (seconds)"))
  if kv("System time (seconds)") is not None:
    metrics["time_sys_seconds"] = float(kv("System time (seconds)"))
  if kv("Percent of CPU this job got") is not None:
    val = kv("Percent of CPU this job got").rstrip("%")
    metrics["cpu_percent"] = float(val)
  if kv("Maximum resident set size (kbytes)") is not None:
    metrics["max_rss_bytes"] = int(kv("Maximum resident set size (kbytes)")) * 1024
elif mode == "posix":
  # time -p: real/user/sys, no RSS.
  m = re.search(r"^real\s+([0-9.]+)\s*$", text, re.M)
  if m:
    metrics["time_real_seconds"] = float(m.group(1))
  m = re.search(r"^user\s+([0-9.]+)\s*$", text, re.M)
  if m:
    metrics["time_user_seconds"] = float(m.group(1))
  m = re.search(r"^sys\s+([0-9.]+)\s*$", text, re.M)
  if m:
    metrics["time_sys_seconds"] = float(m.group(1))

if metrics.get("cpu_percent") is None:
  real = metrics.get("time_real_seconds")
  user = metrics.get("time_user_seconds", 0.0) or 0.0
  sys_t = metrics.get("time_sys_seconds", 0.0) or 0.0
  if isinstance(real, (int, float)) and real and real > 0:
    metrics["cpu_percent"] = (user + sys_t) / real * 100.0

Path(out_json).write_text(json.dumps(metrics, indent=2), encoding="utf-8")
PY
}

run_timed() {
    local label="$1"
    shift
    if [[ "${PERF_ENABLED}" -ne 1 ]]; then
        "$@"
        return
    fi

    local time_raw="${PERF_DIR}/${label}.time.txt"
    local start_epoch end_epoch rc
    start_epoch="$(date +%s)"

    set +e
    case "${TIME_MODE}" in
    bsd)
        { /usr/bin/time -l "$@"; } 2> >(tee "${time_raw}" >&2)
        rc=$?
        ;;
    gnu)
        { /usr/bin/time -v "$@"; } 2> >(tee "${time_raw}" >&2)
        rc=$?
        ;;
    posix)
        { /usr/bin/time -p "$@"; } 2> >(tee "${time_raw}" >&2)
        rc=$?
        ;;
    *)
        "$@"
        rc=$?
        ;;
    esac
    set -e

    end_epoch="$(date +%s)"
    perf_parse_time_to_json "${label}" "${time_raw}" "$((end_epoch - start_epoch))" "${rc}"

    # Split engine stderr away from /usr/bin/time summary so warning audits see
    # only the engine's output. Engine stderr comes first; /usr/bin/time appends
    # its block at the end. Per-engine stderr lands in
    # ${WORK_DIR}/{perl,rust}/<base>.stderr.log when label matches a known engine
    # prefix; otherwise no split (e.g., compare_outputs).
    if [[ -f "${time_raw}" ]]; then
        local engine_dir="" engine_base=""
        case "${label}" in
        perl_vep_*)
            engine_dir="${WORK_DIR}/perl"
            engine_base="${label#perl_vep_}"
            ;;
        rust_vep_*)
            engine_dir="${WORK_DIR}/rust"
            engine_base="${label#rust_vep_}"
            ;;
        esac
        if [[ -n "${engine_dir}" && -n "${engine_base}" ]]; then
            mkdir -p "${engine_dir}"
            local stderr_log="${engine_dir}/${engine_base}.stderr.log"
            # Strip the trailing /usr/bin/time block (lines after the time-summary marker).
            # - bsd: starts with whitespace + numeric "real" line
            # - gnu: "Command being timed:" line
            # - posix: standalone "real <num>" line
            python3 - "${time_raw}" "${stderr_log}" "${TIME_MODE}" <<'PY' 2>/dev/null || true
import sys, re, pathlib

src, dst, mode = sys.argv[1:4]
text = pathlib.Path(src).read_text(encoding="utf-8", errors="replace")
lines = text.splitlines(keepends=True)

# Find first line of /usr/bin/time summary; everything before it is engine stderr.
markers = {
    "gnu": re.compile(r"^\s*Command being timed:"),
    "posix": re.compile(r"^real\s+\d"),
    "bsd": re.compile(r"^\s+\d+\.\d+ real"),
}
pat = markers.get(mode)
cut = len(lines)
if pat is not None:
    for i, ln in enumerate(lines):
        if pat.match(ln):
            cut = i
            break
pathlib.Path(dst).write_text("".join(lines[:cut]), encoding="utf-8")
PY
        fi
    fi
    return "${rc}"
}

perf_parse_docker_stats_to_json() {
    if [[ "${PERF_ENABLED}" -ne 1 ]]; then
        return
    fi
    local stats_tsv="$1"
    local out_json="$2"
    python3 - "${stats_tsv}" "${out_json}" <<'PY'
import json, re, sys
from pathlib import Path

stats_tsv, out_json = sys.argv[1:]
p = Path(stats_tsv)
if not p.exists():
  Path(out_json).write_text(json.dumps({"error": "missing stats file"}, indent=2), encoding="utf-8")
  raise SystemExit(0)

def parse_percent(s):
  s = (s or "").strip()
  if s.endswith("%"):
    s = s[:-1]
  try:
    return float(s)
  except Exception:
    return None

UNIT = {
  "B": 1,
  "KB": 1000,
  "MB": 1000**2,
  "GB": 1000**3,
  "TB": 1000**4,
  "KiB": 1024,
  "MiB": 1024**2,
  "GiB": 1024**3,
  "TiB": 1024**4,
}

def parse_bytes(s):
  s = (s or "").strip()
  m = re.match(r"^([0-9.]+)\s*([A-Za-z]+)$", s)
  if not m:
    return None
  val = float(m.group(1))
  unit = m.group(2)
  mult = UNIT.get(unit)
  if mult is None:
    return None
  return int(val * mult)

max_cpu = None
max_mem = None
samples = 0
for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
  if not line or line.startswith("epoch\t"):
    continue
  parts = line.split("\t")
  if len(parts) < 3:
    continue
  cpu = parse_percent(parts[1])
  mem_usage = parts[2]
  # docker stats prints "X / Y"
  mem_left = mem_usage.split("/", 1)[0].strip()
  mem_b = parse_bytes(mem_left)
  if cpu is not None:
    max_cpu = cpu if max_cpu is None else max(max_cpu, cpu)
  if mem_b is not None:
    max_mem = mem_b if max_mem is None else max(max_mem, mem_b)
  samples += 1

out = {
  "samples": samples,
  "max_cpu_percent": max_cpu,
  "max_mem_bytes": max_mem,
}
Path(out_json).write_text(json.dumps(out, indent=2), encoding="utf-8")
PY
}

perf_write_host_info

RUST_VERSION="$(
    python3 - <<PY
import re, pathlib
text = pathlib.Path("${VEP_RS_DIR}/crates/vep-core/src/lib.rs").read_text(encoding="utf-8")
m = re.search(r"VEP_VERSION:\\s*u32\\s*=\\s*(\\d+)", text)
print(m.group(1) if m else "unknown")
PY
)"

PERL_CACHE_VERSION="$(
    python3 - <<PY
from pathlib import Path
# The cache root holds homo_sapiens/<version>_<assembly>/info.txt; VEP resolves
# --dir_cache the same way, so the version is read from that leaf.
root = Path("${PERL_CACHE_DIR}")
leaves = sorted(root.glob("homo_sapiens/*_${ASSEMBLY}/info.txt"))
info = leaves[-1] if leaves else root / "info.txt"
version = "unknown"
if info.exists():
    for raw in info.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw or raw.startswith("#") or "\\t" not in raw:
            continue
        key, value = raw.split("\\t", 1)
        lk = key.lower().strip()
        if lk in {"version", "cache_version", "api_version", "variation_api_version"} and value.strip().isdigit():
            version = value.strip()
            break
print(version)
PY
)"

echo "preflight:"
echo "  rust vep version constant: ${RUST_VERSION}"
echo "  perl cache version (from info.txt): ${PERL_CACHE_VERSION}"
echo "  rust release build: ${RUST_RELEASE}"
echo "  perl vep mode: ${PERL_VEP_MODE}"
if [[ "${PERL_VEP_MODE}" == "docker" ]]; then
    echo "  perl vep docker image: ${PERL_VEP_DOCKER_IMAGE}"
fi
if [[ "${PERL_VEP_MODE}" == "fixtures" ]]; then
    echo "  fixtures dir: ${FIXTURES_DIR}"
fi
if [[ "${PERL_OUTPUT_CACHE_ENABLED}" -eq 1 ]]; then
    echo "  perl output cache dir: ${PERL_OUTPUT_CACHE_DIR}"
    echo "  perl output cache refresh: ${PERL_OUTPUT_CACHE_REFRESH}"
fi
if [[ "${RUST_VERSION}" != "unknown" && "${PERL_CACHE_VERSION}" != "unknown" && "${RUST_VERSION}" != "${PERL_CACHE_VERSION}" ]]; then
    echo "WARNING: rust version and perl cache version differ; concordance may be degraded" >&2
fi

# Plugin setup

# Resolve the plugin data directory from the ASSEMBLY, not a hardcoded grch37:
# defaulted independently, `--assembly GRCh38` would read GRCh37 plugin data and
# every plugin would annotate nothing, or from the wrong genome build, while the
# comparator reported a pass.
ASSEMBLY_LC="$(echo "${ASSEMBLY}" | tr '[:upper:]' '[:lower:]')"
if [[ -n "${PLUGINS}" && -z "${PLUGIN_DATA_DIR}" ]]; then
    PLUGIN_DATA_DIR="${HOME}/.vep/plugin_data/${ASSEMBLY_LC}/${PLUGIN_TIER}"
fi

# Interlock: an explicitly-passed --plugin-data-dir naming the OTHER assembly is
# almost certainly a mistake, and the failure mode is silent, so refuse it.
if [[ -n "${PLUGINS}" && -n "${PLUGIN_DATA_DIR}" ]]; then
    case "${ASSEMBLY_LC}" in
    grch37)
        if [[ "${PLUGIN_DATA_DIR}" == *grch38* || "${PLUGIN_DATA_DIR}" == *hg38* ]]; then
            echo "ERROR: --assembly ${ASSEMBLY} but --plugin-data-dir looks like GRCh38: ${PLUGIN_DATA_DIR}" >&2
            exit 1
        fi
        ;;
    grch38)
        if [[ "${PLUGIN_DATA_DIR}" == *grch37* || "${PLUGIN_DATA_DIR}" == *hg19* ]]; then
            echo "ERROR: --assembly ${ASSEMBLY} but --plugin-data-dir looks like GRCh37: ${PLUGIN_DATA_DIR}" >&2
            exit 1
        fi
        ;;
    esac
fi

# Build plugin argument strings for Perl and Rust VEP
PERL_PLUGIN_ARGS=()
RUST_PLUGIN_ARGS=()
PLUGIN_NAMES=()

# Emit the `-v` flags the Perl container needs for the requested plugins. Defined
# once for both call sites: build_plugin_args emits the container paths
# /work/loftee_perl and /work/loftee_data for LoFTEE, and a mount list that
# drifts from it leaves a docker-mode LoFTEE comparison unable to run.
plugin_docker_mounts() {
    local -a vols=()
    if [[ ${#PERL_PLUGIN_ARGS[@]} -eq 0 ]]; then
        return 0
    fi

    vols+=(-v "${PLUGIN_DATA_DIR}:/work/plugin_data:ro")
    # Only mount VEP_plugins when it exists: the ensembl-vep image already installs
    # every .pm into /plugins (ENV VEP_DIR_PLUGINS=/plugins), so the mount is a
    # version pin, not a requirement. Mounting a nonexistent host path makes Docker
    # create an empty dir that SHADOWS the image's own plugins.
    if [[ -d "${VEP_PLUGINS_DIR}" ]]; then
        vols+=(-v "${VEP_PLUGINS_DIR}:/opt/vep/.vep/Plugins:ro")
    fi

    # LoFTEE is not in VEP_plugins (the image builds with --skip_plugins LoF), so it
    # needs its own checkout plus its data bundle.
    local p
    for p in "${PLUGIN_NAMES[@]+"${PLUGIN_NAMES[@]}"}"; do
        if [[ "${p}" == "LoFTEE" ]]; then
            local ld="${LOFTEE_DATA_DIR:-${PLUGIN_DATA_DIR}/loftee/${ASSEMBLY_LC}}"
            local lp="${LOFTEE_PERL_PATH:?set LOFTEE_PERL_PATH to a konradjk/loftee checkout}"
            if [[ ! -d "${lp}" ]]; then
                echo "ERROR: LoFTEE requested but the Perl plugin checkout is missing: ${lp}" >&2
                echo "  Clone konradjk/loftee and set LOFTEE_PERL_PATH (branch: master for GRCh37, grch38 for GRCh38)." >&2
                return 1
            fi
            if [[ ! -d "${ld}" ]]; then
                echo "ERROR: LoFTEE requested but the data bundle is missing: ${ld}" >&2
                return 1
            fi
            vols+=(-v "${lp}:/work/loftee_perl:ro")
            vols+=(-v "${ld}:/work/loftee_data:ro")
            break
        fi
    done

    printf '%s\n' "${vols[@]}"
}

build_plugin_args() {
    # Map plugin name to --plugin argument for a given data directory.
    # $1 = plugin name, $2 = data directory (host path), $3 = "perl" or "rust"
    # For Docker Perl VEP, paths are remapped to /work/plugin_data/
    # File names come from scripts/data/plugin_data_layout.sh, the table
    # setup_plugin_data.sh stages from, so the harness reads what staging wrote.
    local plugin="$1"
    local data_dir="$2"
    local target="$3"
    local docker_data="/work/plugin_data"

    local d="${data_dir}"
    if [[ "${target}" == "perl" && "${PERL_VEP_MODE}" == "docker" ]]; then
        d="${docker_data}"
    fi
    local f

    case "${plugin}" in
    CADD)
        echo "CADD,snv=${d}/$(plugin_data_file CADD snv "${ASSEMBLY_LC}"),indels=${d}/$(plugin_data_file CADD indel "${ASSEMBLY_LC}")"
        ;;
    REVEL)
        echo "REVEL,file=${d}/$(plugin_data_file REVEL file "${ASSEMBLY_LC}")"
        ;;
    SpliceAI)
        # The indel file is passed only when staged (it is not freely distributed);
        # Perl's SpliceAI.pm adds both files and fails on an absent one, vep-rs accepts
        # either alone, so a run without it is a vep-rs-only run of this plugin.
        f="SpliceAI,snv=${d}/$(plugin_data_file SpliceAI snv "${ASSEMBLY_LC}")"
        if [[ -e "${data_dir}/$(plugin_data_file SpliceAI indel "${ASSEMBLY_LC}")" ]]; then
            f="${f},indel=${d}/$(plugin_data_file SpliceAI indel "${ASSEMBLY_LC}")"
        fi
        echo "${f}"
        ;;
    gnomADc)
        echo "gnomADc,${d}/$(plugin_data_file gnomADc file "${ASSEMBLY_LC}")"
        ;;
    AlphaMissense)
        echo "AlphaMissense,${d}/$(plugin_data_file AlphaMissense file "${ASSEMBLY_LC}")"
        ;;
    dbNSFP)
        echo "dbNSFP,${d}/$(plugin_data_file dbNSFP file "${ASSEMBLY_LC}"),SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred"
        ;;
    dbscSNV)
        echo "dbscSNV,${d}/$(plugin_data_file dbscSNV file "${ASSEMBLY_LC}")"
        ;;
    LoFtool)
        echo "LoFtool,${d}/$(plugin_data_file LoFtool file "${ASSEMBLY_LC}")"
        ;;
    pLI)
        echo "pLI,${d}/$(plugin_data_file pLI file "${ASSEMBLY_LC}")"
        ;;
    GWAS)
        # NHGRI-EBI GWAS Catalog. GRCh38 ONLY: the catalog publishes CHR_POS in
        # GRCh38 coordinates, so a GRCh37 run would compare against the wrong
        # genome. vep-rs reads the TSV (gzip or plain); Perl's GWAS.pm wants
        # `file=` plus `type=curated`.
        if [[ "${ASSEMBLY_LC}" != "grch38" ]]; then
            echo "ERROR: GWAS is GRCh38-only (the catalog publishes GRCh38 coordinates); got ${ASSEMBLY}" >&2
            return 1
        fi
        if [[ "${target}" == "perl" ]]; then
            echo "GWAS,file=${d}/$(plugin_data_file GWAS file "${ASSEMBLY_LC}"),type=curated"
        else
            echo "GWAS,file=${d}/$(plugin_data_file GWAS file "${ASSEMBLY_LC}")"
        fi
        ;;
    LoFTEE)
        # LoFTEE data lives under a dedicated dir (LOFTEE_DATA_DIR, default
        # ${d}/loftee/<assembly>), not the tabix tier layout. Perl plugin is
        # "LoF" (key:value args + loftee_path); vep-rs is "LoFTEE" (key=value).
        # Both are fed the SAME data so HC/LC uses identical inputs.
        #
        # The two upstream bundles differ STRUCTURALLY: GRCh38 ships GERP as a
        # bigWig, GRCh37 as a tabix per-base TSV. Hardcoding the bigWig made a
        # GRCh37 invocation impossible to build.
        local ld="${LOFTEE_DATA_DIR:-${d}/loftee/${ASSEMBLY_LC}}"
        local gerp_file gerp_key ancestor_plain ancestor_bgzf
        gerp_file="$(basename "$(plugin_data_file LoFTEE gerp "${ASSEMBLY_LC}")")"
        # The Perl plugin reads the bgzipped ancestral FASTA (with .fai and .gzi); vep-rs
        # reads the plain one. Both are staged, so both names come from the layout table.
        ancestor_plain="$(basename "$(plugin_data_file LoFTEE ancestor "${ASSEMBLY_LC}")")"
        ancestor_bgzf="$(basename "$(plugin_data_file LoFTEE ancestor_bgzf "${ASSEMBLY_LC}")")"
        if [[ "${ASSEMBLY_LC}" == "grch37" ]]; then
            gerp_key="gerp_tabix"
        else
            gerp_key="gerp_bigwig"
        fi
        if [[ "${target}" == "perl" ]]; then
            local lp="${LOFTEE_PERL_PATH:-./loftee_perl}"
            local pld="${ld}"
            if [[ "${PERL_VEP_MODE}" == "docker" ]]; then
                lp="/work/loftee_perl"
                pld="/work/loftee_data"
            fi
            # Perl reads gerp_bigwig on the grch38 branch, gerp_file on master.
            if [[ "${ASSEMBLY_LC}" == "grch37" ]]; then
                echo "LoF,loftee_path:${lp},human_ancestor_fa:${pld}/${ancestor_bgzf},gerp_file:${pld}/${gerp_file}"
            else
                echo "LoF,loftee_path:${lp},human_ancestor_fa:${pld}/${ancestor_bgzf},gerp_bigwig:${pld}/${gerp_file}"
            fi
        else
            echo "LoFTEE,human_ancestor_fa=${ld}/${ancestor_plain},${gerp_key}=${ld}/${gerp_file}"
        fi
        ;;
    *)
        echo "ERROR: unknown plugin: ${plugin}" >&2
        return 1
        ;;
    esac
}

if [[ -n "${PLUGINS}" ]]; then
    IFS=',' read -ra PLUGIN_NAMES <<<"${PLUGINS}"

    if [[ ! -d "${PLUGIN_DATA_DIR}" ]]; then
        echo "ERROR: plugin data directory not found: ${PLUGIN_DATA_DIR}" >&2
        echo "  Run: scripts/data/setup_plugin_data.sh --assembly ${ASSEMBLY} --tier ${PLUGIN_TIER} --data-dir ${PLUGIN_DATA_DIR}" >&2
        exit 1
    fi

    echo "plugins:"
    echo "  plugins:   ${PLUGIN_NAMES[*]}"
    echo "  data dir:  ${PLUGIN_DATA_DIR}"
    echo "  tier:      ${PLUGIN_TIER}"

    for plugin in "${PLUGIN_NAMES[@]}"; do
        perl_arg="$(build_plugin_args "${plugin}" "${PLUGIN_DATA_DIR}" "perl")"
        rust_arg="$(build_plugin_args "${plugin}" "${PLUGIN_DATA_DIR}" "rust")"
        PERL_PLUGIN_ARGS+=(--plugin "${perl_arg}")
        RUST_PLUGIN_ARGS+=(--plugin "${rust_arg}")
        echo "  ${plugin}: perl='${perl_arg}'"
    done
fi

if [[ -z "${JSON_CACHE_DIR}" ]]; then
    JSON_CACHE_DIR="${WORK_DIR}/cache_json"
    echo "ERROR: --json-cache-dir is required; no JSON cache was given (for example ${JSON_CACHE_DIR})." >&2
    echo "  Build one with vep-cache-builder (no Perl required):" >&2
    echo "    vep-cache-builder --release 115 --assembly GRCh38 \\" >&2
    echo "      --genome-fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \\" >&2
    echo "      --include-predictions --output-dir ${JSON_CACHE_DIR}" >&2
    echo "  Or provide an existing cache with --json-cache-dir." >&2
    exit 1
else
    JSON_CACHE_DIR="$(cd "${JSON_CACHE_DIR}" && pwd)"
    if [[ ! -f "${JSON_CACHE_DIR}/info.json" ]]; then
        echo "ERROR: --json-cache-dir must contain info.json: ${JSON_CACHE_DIR}" >&2
        exit 1
    fi
fi

PREP_ARGS=(
    --input-dir "${BENCHMARK_DIR}"
    --output-dir "${WORK_DIR}/inputs"
    --chrom-mode auto
    --cache-dir "${JSON_CACHE_DIR}"
    --assembly "${ASSEMBLY}"
)
if [[ "${MODE}" == "smoke" ]]; then
    PREP_ARGS+=(--max-variants "${SMOKE_VARIANTS}")
fi
if [[ "${CANONICAL_CONTIGS}" -eq 1 ]]; then
    PREP_ARGS+=(--canonical-contigs)
else
    PREP_ARGS+=(--no-canonical-contigs)
fi

echo "preparing benchmark vcfs (${MODE})"
run_timed "prepare_benchmark_vcfs" python3 "${VEP_RS_DIR}/scripts/concordance/prepare_benchmark_vcfs.py" "${PREP_ARGS[@]}"

echo "running annotations"
shopt -s nullglob
prepared_inputs=("${WORK_DIR}"/inputs/*.vcf)
shopt -u nullglob
if ((${#prepared_inputs[@]} == 0)); then
    echo "ERROR: no prepared VCFs found in ${WORK_DIR}/inputs" >&2
    exit 1
fi

# ---- Provenance pinning ------------------------------------------------------
# Pin the exact bits this run used, so the result can be reproduced.
echo "capturing provenance pins (Docker digest, cache SHAs, FASTA SHA, input SHAs)..."
PERL_DOCKER_DIGEST=""
if [[ "${PERL_VEP_MODE}" == "docker" ]] && command -v docker >/dev/null 2>&1; then
    docker pull "${PERL_VEP_DOCKER_IMAGE}" >/dev/null 2>&1 || true
    PERL_DOCKER_DIGEST="$(docker inspect --format='{{index .RepoDigests 0}}' "${PERL_VEP_DOCKER_IMAGE}" 2>/dev/null || echo "")"
fi

dir_sha256() {
    # sha256 over the contents of a directory; deterministic across reruns.
    local dir="$1"
    if [[ -d "${dir}" ]]; then
        (cd "${dir}" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum 2>/dev/null | sha256sum | awk '{print $1}')
    fi
}

PERL_CACHE_SHA=""
[[ -n "${PERL_CACHE_DIR}" && -d "${PERL_CACHE_DIR}" ]] && PERL_CACHE_SHA="$(dir_sha256 "${PERL_CACHE_DIR}")"
JSON_CACHE_SHA=""
[[ -n "${JSON_CACHE_DIR}" && -d "${JSON_CACHE_DIR}" ]] && JSON_CACHE_SHA="$(dir_sha256 "${JSON_CACHE_DIR}")"
FASTA_SHA=""
[[ -n "${FASTA}" && -f "${FASTA}" ]] && FASTA_SHA="$(sha256sum "${FASTA}" 2>/dev/null | awk '{print $1}')"

# Per-VCF SHAs as a TSV side-channel; the provenance.json writer below reads it.
# The side-channel keeps escaping/quoting in Python, not bash, and a plain loop
# needs no associative array (bash 3.2, which macOS ships, has none).
input_shas_tsv="${WORK_DIR}/reports/_input_shas.tsv"
: >"${input_shas_tsv}"
for vcf in "${prepared_inputs[@]}"; do
    printf '%s\t%s\n' "$(basename "${vcf}")" "$(sha256sum "${vcf}" 2>/dev/null | awk '{print $1}')" >>"${input_shas_tsv}"
done
export PERL_DOCKER_DIGEST PERL_CACHE_SHA JSON_CACHE_SHA FASTA_SHA

# Cache key for reusable Perl output.
#
# `plugins=` alone is not enough: it records WHICH plugins ran, not which DATA
# they read. Two runs naming the same plugins but pointing at different files
# (GRCh37 vs GRCh38 CADD, or a re-prepared REVEL index) collided, and the second
# silently reused the first's output.
PLUGIN_DATA_FINGERPRINT=""
if [[ -n "${PLUGINS}" && -d "${PLUGIN_DATA_DIR}" ]]; then
    # Size+mtime per file, not a content hash: these are multi-GB tabix files and
    # hashing them would dominate the run. `stat` is not portable (BSD -f vs
    # GNU -c), so probe once and use whichever the host accepts.
    if stat -f '%z' "${PLUGIN_DATA_DIR}" >/dev/null 2>&1; then
        _stat_fmt=(-f '%N:%z:%m')
    else
        _stat_fmt=(-c '%n:%s:%Y')
    fi
    if command -v sha256sum >/dev/null 2>&1; then
        _sha_cmd=(sha256sum)
    else
        _sha_cmd=(shasum -a 256)
    fi
    PLUGIN_DATA_FINGERPRINT="$(
        find "${PLUGIN_DATA_DIR}" -type f \
            \( -name '*.gz' -o -name '*.bgz' -o -name '*.txt' -o -name '*.tsv' -o -name '*.bw' \) \
            -exec stat "${_stat_fmt[@]}" {} \; 2>/dev/null | sort | "${_sha_cmd[@]}" | cut -d' ' -f1
    )"
fi
PERL_VEP_SIGNATURE="mode=${PERL_VEP_MODE};docker_image=${PERL_VEP_DOCKER_IMAGE};perl_cache_dir=${PERL_CACHE_DIR};species=${SPECIES};assembly=${ASSEMBLY};format=vcf;buffer_size=5000;fork=1;no_headers=1;no_stats=1;offline=1;cache=1;plugins=${PLUGINS};plugin_data_dir=${PLUGIN_DATA_DIR};plugin_data=${PLUGIN_DATA_FINGERPRINT};fasta=${FASTA:-}"
PERL_VEP_SIGNATURE_SHA="$(sha256_string "${PERL_VEP_SIGNATURE}")"

for input_vcf in "${prepared_inputs[@]}"; do
    base_name="$(basename "${input_vcf}" .vcf)"
    perl_out="${WORK_DIR}/perl/${base_name}.txt"
    rust_out="${WORK_DIR}/rust/${base_name}.txt"

    echo "  perl vep: ${base_name}"

    if [[ "${PERL_VEP_MODE}" == "fixtures" ]]; then
        fixture_file="${FIXTURES_DIR}/${base_name}.txt"
        if [[ ! -f "${fixture_file}" ]]; then
            echo "ERROR: fixture file not found: ${fixture_file}" >&2
            exit 1
        fi
        echo "    -> fixture (${fixture_file})"
        cp "${fixture_file}" "${perl_out}"
    else

        perl_cache_hit=0
        cache_key=""
        cache_txt=""
        if [[ "${PERL_OUTPUT_CACHE_ENABLED}" -eq 1 ]]; then
            input_sha="$(sha256_file "${input_vcf}")"
            cache_key="$(sha256_string "${input_sha}::${PERL_VEP_SIGNATURE_SHA}")"
            cache_txt="${PERL_OUTPUT_CACHE_DIR}/${cache_key}.txt"
            if [[ "${PERL_OUTPUT_CACHE_REFRESH}" -eq 0 && -s "${cache_txt}" ]]; then
                perl_cache_hit=1
                echo "    -> cached (${cache_key:0:12})"
                if [[ "${PERF_ENABLED}" -eq 1 ]]; then
                    run_timed "perl_copy_${base_name}" cp "${cache_txt}" "${perl_out}"
                else
                    cp "${cache_txt}" "${perl_out}"
                fi
            fi
        fi

        if [[ "${perl_cache_hit}" -eq 0 ]]; then
            if [[ "${PERL_VEP_MODE}" == "docker" && "${PERF_ENABLED}" -eq 1 ]]; then
                # A named container, so docker stats can sample it while it runs.
                container_name_raw="vep_concordance_${base_name}_$RANDOM"
                container_name="${container_name_raw//[^a-zA-Z0-9_.-]/_}"
                stats_tsv="${PERF_DIR}/perl_${base_name}.docker_stats.tsv"
                stats_json="${PERF_DIR}/perl_${base_name}.docker_stats.json"
                echo -e "epoch\tcpu_perc\tmem_usage\tmem_perc\tnet_io\tblock_io\tpids" >"${stats_tsv}"

                cleanup_container() {
                    docker rm -f "${container_name}" >/dev/null 2>&1 || true
                }
                trap cleanup_container EXIT

                start_epoch="$(date +%s)"
                # Build Docker volume mounts for plugins
                plugin_docker_vols=()
                if [[ ${#PERL_PLUGIN_ARGS[@]} -gt 0 ]]; then
                    # Captured first: a failure inside a process substitution is invisible
                    # to the loop, so a missing LoFTEE checkout would start the container
                    # with no plugin mounts instead of stopping the run.
                    _mounts="$(plugin_docker_mounts)" || exit 1
                    while IFS= read -r _v; do
                        [[ -n "${_v}" ]] && plugin_docker_vols+=("${_v}")
                    done <<<"${_mounts}"
                fi

                # FASTA parity: when FASTA is set, both engines must receive --fasta.
                # Asymmetry hides 5'UTR start_lost and 3'UTR stop_retained DNA-level divergences.
                perl_fasta_vols=()
                perl_fasta_args=()
                if [[ -n "${FASTA}" ]]; then
                    perl_fasta_vols+=(-v "${FASTA}:/work/fasta.fa:ro")
                    if [[ -f "${FASTA}.fai" ]]; then
                        perl_fasta_vols+=(-v "${FASTA}.fai:/work/fasta.fa.fai:ro")
                    fi
                    perl_fasta_args+=(--fasta /work/fasta.fa)
                fi

                container_id="$(docker run -d \
                    --name "${container_name}" \
                    --user "$(id -u):$(id -g)" \
                    -v "${input_vcf}:/work/input.vcf:ro" \
                    -v "${WORK_DIR}/perl:/work/perl_out" \
                    -v "${PERL_CACHE_DIR}:/work/cache:ro" \
                    "${perl_fasta_vols[@]+"${perl_fasta_vols[@]}"}" \
                    "${plugin_docker_vols[@]+"${plugin_docker_vols[@]}"}" \
                    "${PERL_VEP_DOCKER_IMAGE}" \
                    vep \
                    -i /work/input.vcf \
                    -o "/work/perl_out/${base_name}.txt" \
                    --offline \
                    --cache \
                    --dir_cache /work/cache \
                    --species "${SPECIES}" \
                    --assembly "${ASSEMBLY}" \
                    --format vcf \
                    --buffer_size 5000 \
                    --fork 1 \
                    --no_headers \
                    --no_stats \
                    --quiet \
                    --force_overwrite \
                    "${perl_fasta_args[@]+"${perl_fasta_args[@]}"}" \
                    "${PERL_PLUGIN_ARGS[@]+"${PERL_PLUGIN_ARGS[@]}"}")"

                # Take an immediate snapshot; fast runs can finish before the sampler loop
                # gets scheduled, leaving an empty stats file.
                ts="$(date +%s)"
                stats="$(docker stats --no-stream --format '{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}\t{{.NetIO}}\t{{.BlockIO}}\t{{.PIDs}}' "${container_name}" 2>/dev/null || true)"
                echo -e "${ts}\t${stats}" >>"${stats_tsv}"

                (
                    while docker inspect -f '{{.State.Running}}' "${container_name}" >/dev/null 2>&1; do
                        running="$(docker inspect -f '{{.State.Running}}' "${container_name}" 2>/dev/null || echo "false")"
                        ts="$(date +%s)"
                        stats="$(docker stats --no-stream --format '{{.CPUPerc}}\t{{.MemUsage}}\t{{.MemPerc}}\t{{.NetIO}}\t{{.BlockIO}}\t{{.PIDs}}' "${container_name}" 2>/dev/null || true)"
                        echo -e "${ts}\t${stats}" >>"${stats_tsv}"
                        if [[ "${running}" != "true" ]]; then
                            break
                        fi
                        sleep "${PERF_SAMPLE_SECONDS}"
                    done
                ) &
                stats_pid=$!

                set +e
                perl_rc="$(docker wait "${container_name}")"
                rc=$?
                set -e
                end_epoch="$(date +%s)"

                kill "${stats_pid}" >/dev/null 2>&1 || true
                wait "${stats_pid}" >/dev/null 2>&1 || true
                perf_parse_docker_stats_to_json "${stats_tsv}" "${stats_json}"

                # Always capture container stderr/logs for warning audits.
                docker logs "${container_name}" >"${WORK_DIR}/perl/${base_name}.stderr.log" 2>&1 || true
                if [[ "${rc}" -ne 0 || "${perl_rc}" -ne 0 ]]; then
                    cp "${WORK_DIR}/perl/${base_name}.stderr.log" "${PERF_DIR}/perl_${base_name}.docker.log" 2>/dev/null || true
                fi
                cleanup_container
                trap - EXIT

                python3 - <<PY
import json
from pathlib import Path
out = Path("${PERF_DIR}") / "perl_${base_name}.run.json"
out.write_text(json.dumps({
  "base_name": "${base_name}",
  "container_id": "${container_id}",
  "container_name": "${container_name}",
  "exit_code": int("${perl_rc}") if "${perl_rc}".isdigit() else None,
  "wall_seconds": int("${end_epoch}") - int("${start_epoch}"),
  "stats_tsv": "${stats_tsv}",
  "stats_json": "${stats_json}",
}, indent=2), encoding="utf-8")
PY
                # Ensure docker wait succeeded.
                if [[ "${rc}" -ne 0 || "${perl_rc}" -ne 0 ]]; then
                    echo "ERROR: perl vep docker run failed for ${base_name}" >&2
                    exit 1
                fi
            else
                if [[ "${PERL_VEP_MODE}" == "docker" ]]; then
                    # Build Docker volume mounts for plugins
                    plugin_docker_vols=()
                    if [[ ${#PERL_PLUGIN_ARGS[@]} -gt 0 ]]; then
                        # Captured first: a failure inside a process substitution is invisible
                        # to the loop, so a missing LoFTEE checkout would start the container
                        # with no plugin mounts instead of stopping the run.
                        _mounts="$(plugin_docker_mounts)" || exit 1
                        while IFS= read -r _v; do
                            [[ -n "${_v}" ]] && plugin_docker_vols+=("${_v}")
                        done <<<"${_mounts}"
                    fi

                    # FASTA parity: when FASTA is set, both engines must receive --fasta.
                    perl_fasta_vols=()
                    perl_fasta_args=()
                    if [[ -n "${FASTA}" ]]; then
                        perl_fasta_vols+=(-v "${FASTA}:/work/fasta.fa:ro")
                        if [[ -f "${FASTA}.fai" ]]; then
                            perl_fasta_vols+=(-v "${FASTA}.fai:/work/fasta.fa.fai:ro")
                        fi
                        perl_fasta_args+=(--fasta /work/fasta.fa)
                    fi

                    run_timed "perl_vep_${base_name}" docker run --rm \
                        --user "$(id -u):$(id -g)" \
                        -v "${input_vcf}:/work/input.vcf:ro" \
                        -v "${WORK_DIR}/perl:/work/perl_out" \
                        -v "${PERL_CACHE_DIR}:/work/cache:ro" \
                        "${perl_fasta_vols[@]+"${perl_fasta_vols[@]}"}" \
                        "${plugin_docker_vols[@]+"${plugin_docker_vols[@]}"}" \
                        "${PERL_VEP_DOCKER_IMAGE}" \
                        vep \
                        -i /work/input.vcf \
                        -o "/work/perl_out/${base_name}.txt" \
                        --offline \
                        --cache \
                        --dir_cache /work/cache \
                        --species "${SPECIES}" \
                        --assembly "${ASSEMBLY}" \
                        --format vcf \
                        --buffer_size 5000 \
                        --fork 1 \
                        --no_headers \
                        --no_stats \
                        --quiet \
                        --force_overwrite \
                        "${perl_fasta_args[@]+"${perl_fasta_args[@]}"}" \
                        "${PERL_PLUGIN_ARGS[@]+"${PERL_PLUGIN_ARGS[@]}"}"
                else
                    # FASTA parity: when FASTA is set, both engines must receive --fasta.
                    perl_fasta_args=()
                    if [[ -n "${FASTA}" ]]; then
                        perl_fasta_args+=(--fasta "${FASTA}")
                    fi

                    run_timed "perl_vep_${base_name}" perl "${ENSEMBL_VEP_DIR}/vep" \
                        -i "${input_vcf}" \
                        -o "${perl_out}" \
                        --offline \
                        --cache \
                        --dir_cache "${PERL_CACHE_DIR}" \
                        --species "${SPECIES}" \
                        --assembly "${ASSEMBLY}" \
                        --format vcf \
                        --buffer_size 5000 \
                        --fork 1 \
                        --no_headers \
                        --no_stats \
                        --quiet \
                        --force_overwrite \
                        --dir_plugins "${VEP_PLUGINS_DIR}" \
                        "${perl_fasta_args[@]+"${perl_fasta_args[@]}"}" \
                        "${PERL_PLUGIN_ARGS[@]+"${PERL_PLUGIN_ARGS[@]}"}"
                fi
            fi

            if [[ "${PERL_OUTPUT_CACHE_ENABLED}" -eq 1 ]]; then
                # Compute the cache key if the earlier pass did not.
                if [[ -z "${cache_key}" ]]; then
                    input_sha="$(sha256_file "${input_vcf}")"
                    cache_key="$(sha256_string "${input_sha}::${PERL_VEP_SIGNATURE_SHA}")"
                    cache_txt="${PERL_OUTPUT_CACHE_DIR}/${cache_key}.txt"
                fi
                cache_tmp="${cache_txt}.tmp.$$"
                cp "${perl_out}" "${cache_tmp}"
                mv "${cache_tmp}" "${cache_txt}"
                python3 - <<PY
import json
from datetime import datetime, timezone
from pathlib import Path

meta = {
  "generated_at_utc": datetime.now(timezone.utc).isoformat(),
  "cache_key": "${cache_key}",
  "cache_key_prefix": "${cache_key:0:12}",
  "input_vcf_sha256": "${input_sha}",
  "base_name": "${base_name}",
  "perl_vep_signature": "${PERL_VEP_SIGNATURE}",
}
out = Path("${PERL_OUTPUT_CACHE_DIR}") / f"{meta['cache_key']}.json"
out.write_text(json.dumps(meta, indent=2), encoding="utf-8")
PY
                echo "    -> stored (${cache_key:0:12})"
            fi
        fi
    fi # end fixtures else

    echo "  rust vep: ${base_name}"
    (
        cd "${VEP_RS_DIR}"
        rust_extra_args=()
        if [[ -n "${FASTA}" ]]; then
            rust_extra_args+=(--fasta "${FASTA}")
        fi
        cargo_profile_args=()
        if [[ "${RUST_RELEASE}" -eq 1 ]]; then
            cargo_profile_args+=(--release)
        fi
        run_timed "rust_vep_${base_name}" cargo run -q "${cargo_profile_args[@]}" -p vep-cli -- \
            -i "${input_vcf}" \
            -o "${rust_out}" \
            --offline \
            "${rust_extra_args[@]}" \
            --json_cache "${JSON_CACHE_DIR}" \
            --species "${SPECIES}" \
            --assembly "${ASSEMBLY}" \
            --format vcf \
            --buffer_size 5000 \
            --fork 1 \
            --no_headers \
            --no_stats \
            --force \
            --quiet \
            "${RUST_PLUGIN_ARGS[@]+"${RUST_PLUGIN_ARGS[@]}"}"
    )
done

python3 - <<PY
import json, os
from datetime import datetime, timezone
from pathlib import Path

provenance = {
    "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    "benchmark_dir": "${BENCHMARK_DIR}",
    "ensembl_vep_dir": "${ENSEMBL_VEP_DIR}",
    "vep_rs_dir": "${VEP_RS_DIR}",
    "perl_cache_dir": "${PERL_CACHE_DIR}",
    "json_cache_dir": "${JSON_CACHE_DIR}",
    "fasta": "${FASTA}" if "${FASTA}" else None,
    "perl_vep_mode": "${PERL_VEP_MODE}",
    "perl_vep_docker_image": "${PERL_VEP_DOCKER_IMAGE}" if "${PERL_VEP_MODE}" == "docker" else None,
    "fixtures_dir": "${FIXTURES_DIR}" if "${PERL_VEP_MODE}" == "fixtures" else None,
    "mode": "${MODE}",
    "smoke_variants": int("${SMOKE_VARIANTS}"),
    "species": "${SPECIES}",
    "assembly": "${ASSEMBLY}",
    "rust_version_constant": "${RUST_VERSION}",
    "rust_release": bool(int("${RUST_RELEASE}")),
    "perl_cache_version": "${PERL_CACHE_VERSION}",
    "plugins": "${PLUGINS}" if "${PLUGINS}" else None,
    "plugin_data_dir": "${PLUGIN_DATA_DIR}" if "${PLUGIN_DATA_DIR}" else None,
    "plugin_tier": "${PLUGIN_TIER}" if "${PLUGINS}" else None,
    "perl_docker_digest": os.environ.get("PERL_DOCKER_DIGEST") or None,
    "perl_cache_sha256": os.environ.get("PERL_CACHE_SHA") or None,
    "json_cache_sha256": os.environ.get("JSON_CACHE_SHA") or None,
    "fasta_sha256": os.environ.get("FASTA_SHA") or None,
    "fasta_required_gate": True,
    "no_fasta_opt_out": bool(int("${NO_FASTA}")),
    "canonical_contigs_filter_active": bool(int("${CANONICAL_CONTIGS}")),
    "no_canonical_contigs_opt_out": bool(int("${NO_CANONICAL_CONTIGS}")),
    "run_id": "${RUN_ID}" or None,
    "host_id": "${HOST_ID}" or None,
    "input_vcf_sha256s": dict(
        line.rstrip("\n").split("\t", 1)
        for line in Path("${WORK_DIR}/reports/_input_shas.tsv").read_text(encoding="utf-8").splitlines()
        if line
    ),
}
# Best-effort: pull the input-prep manifest if it exists and merge canonical-contigs
# audit fields into provenance.json (per_file dropped counts + filtered sha256s).
prep_manifest_path = Path("${WORK_DIR}/inputs/manifest.json")
if prep_manifest_path.exists():
    try:
        prep_manifest = json.loads(prep_manifest_path.read_text(encoding="utf-8"))
        provenance["input_prep_manifest"] = {
            "canonical_contigs_filter_active": prep_manifest.get("canonical_contigs_filter_active"),
            "canonical_contigs_set": prep_manifest.get("canonical_contigs_set"),
            "assembly": prep_manifest.get("assembly"),
            "files": [
                {
                    "input_file": f.get("input_file"),
                    "output_file": f.get("output_file"),
                    "variants_written": f.get("variants_written"),
                    "variants_dropped_non_canonical": f.get("variants_dropped_non_canonical", 0),
                    "output_sha256": f.get("output_sha256", ""),
                }
                for f in prep_manifest.get("files", [])
            ],
        }
    except Exception as exc:
        provenance["input_prep_manifest_error"] = str(exc)
out = Path("${WORK_DIR}/reports/provenance.json")
out.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
print(f"wrote provenance: {out}")
PY

echo "comparing outputs"
COMPARE_ARGS=(
    --perl-dir "${WORK_DIR}/perl"
    --rust-dir "${WORK_DIR}/rust"
    --report-dir "${WORK_DIR}/reports"
    --assembly "${ASSEMBLY}"
)
if [[ "${CANONICAL_CONTIGS}" -eq 1 ]]; then
    COMPARE_ARGS+=(--canonical-contigs)
else
    COMPARE_ARGS+=(--no-canonical-contigs)
fi
run_timed "compare_outputs" python3 "${VEP_RS_DIR}/scripts/concordance/compare_vep_outputs.py" "${COMPARE_ARGS[@]}"

# Plugin concordance comparison
if [[ ${#PLUGIN_NAMES[@]} -gt 0 ]]; then
    echo "comparing plugin outputs"
    PLUGIN_LIST_CSV="$(
        IFS=','
        echo "${PLUGIN_NAMES[*]}"
    )"
    # No `|| true`: the comparator's exit status IS the plugin verdict. Swallowing
    # it meant a plugin FAIL -- including one that annotated nothing at all --
    # could never fail the run.
    PLUGIN_COMPARE_STATUS=0
    run_timed "compare_plugin_outputs" python3 "${VEP_RS_DIR}/scripts/concordance/compare_plugin_outputs.py" \
        --perl-dir "${WORK_DIR}/perl" \
        --rust-dir "${WORK_DIR}/rust" \
        --plugins "${PLUGIN_LIST_CSV}" \
        --min-annotated "${PLUGIN_MIN_ANNOTATED}" \
        --json-output "${WORK_DIR}/reports/plugin_concordance.json" || PLUGIN_COMPARE_STATUS=$?

    if [[ "${PLUGIN_COMPARE_STATUS}" -ne 0 ]]; then
        echo "ERROR: [plugin_concordance] comparison FAILED (exit ${PLUGIN_COMPARE_STATUS}); see ${WORK_DIR}/reports/plugin_concordance.json" >&2
        # Carry it to the process exit. Printing the ERROR is not enough: this
        # script's last statement is an `echo`, so a status that stops here leaves
        # the run reporting success and a plugin that annotated nothing reading
        # as a pass -- the same false green the `|| true` produced, one layer up.
        RUN_STATUS="${PLUGIN_COMPARE_STATUS}"
    fi

    # Plugin performance benchmark: baseline (no plugins) vs with plugins
    if [[ "${PERF_ENABLED}" -eq 1 ]]; then
        echo "running plugin performance benchmark (baseline vs with-plugins)"
        mkdir -p "${WORK_DIR}/rust_baseline"

        # Pick the first prepared input for the benchmark
        benchmark_vcf="${prepared_inputs[0]}"
        benchmark_base="$(basename "${benchmark_vcf}" .vcf)"

        (
            cd "${VEP_RS_DIR}"
            cargo_profile_args=()
            if [[ "${RUST_RELEASE}" -eq 1 ]]; then
                cargo_profile_args+=(--release)
            fi
            rust_extra_args=()
            if [[ -n "${FASTA}" ]]; then
                rust_extra_args+=(--fasta "${FASTA}")
            fi

            # Baseline run (no plugins)
            run_timed "rust_baseline_${benchmark_base}" cargo run -q "${cargo_profile_args[@]}" -p vep-cli -- \
                -i "${benchmark_vcf}" \
                -o "${WORK_DIR}/rust_baseline/${benchmark_base}.txt" \
                --offline \
                "${rust_extra_args[@]}" \
                --json_cache "${JSON_CACHE_DIR}" \
                --species "${SPECIES}" \
                --assembly "${ASSEMBLY}" \
                --format vcf \
                --buffer_size 5000 \
                --fork 1 \
                --no_headers \
                --no_stats \
                --force \
                --quiet
        )

        # Generate plugin timing report
        python3 - <<PY
import json
from pathlib import Path

perf_dir = Path("${PERF_DIR}")
report = {}

# Parse baseline timing
baseline_file = perf_dir / "rust_baseline_${benchmark_base}.time.json"
if baseline_file.exists():
    baseline = json.loads(baseline_file.read_text())
    report["baseline_rust_seconds"] = baseline.get("wall_seconds", 0)

# Parse with-plugins timing
plugins_file = perf_dir / "rust_vep_${benchmark_base}.time.json"
if plugins_file.exists():
    with_plugins = json.loads(plugins_file.read_text())
    report["with_plugins_rust_seconds"] = with_plugins.get("wall_seconds", 0)

# Parse Perl timing
perl_file = perf_dir / "perl_vep_${benchmark_base}.time.json"
if perl_file.exists():
    perl_data = json.loads(perl_file.read_text())
    report["perl_with_plugins_seconds"] = perl_data.get("wall_seconds", 0)

# Compute derived metrics
b = report.get("baseline_rust_seconds", 0)
p = report.get("with_plugins_rust_seconds", 0)
if b > 0 and p > 0:
    report["plugin_overhead_seconds"] = round(p - b, 2)
    report["plugin_overhead_pct"] = round((p - b) / b * 100, 1)

perl_s = report.get("perl_with_plugins_seconds", 0)
if perl_s > 0 and p > 0:
    report["speedup_vs_perl"] = f"{perl_s / p:.1f}x"

report["benchmark_file"] = "${benchmark_base}"
report["plugins"] = "${PLUGIN_LIST_CSV}"

out = perf_dir / "plugin_timing.json"
out.write_text(json.dumps(report, indent=2), encoding="utf-8")
print(f"  plugin timing: {out}")
PY
    fi
fi

echo "done"
echo "  summary:    ${WORK_DIR}/reports/summary.json"
echo "  markdown:   ${WORK_DIR}/reports/summary.md"
echo "  discordant: ${WORK_DIR}/reports/discordant.tsv"
echo "  provenance: ${WORK_DIR}/reports/provenance.json"
if [[ ${#PLUGIN_NAMES[@]} -gt 0 ]]; then
    echo "  plugins:    ${WORK_DIR}/reports/plugin_concordance.json"
    if [[ "${PERF_ENABLED}" -eq 1 ]]; then
        echo "  perf:       ${PERF_DIR}/plugin_timing.json"
    fi
fi

exit "${RUN_STATUS}"
