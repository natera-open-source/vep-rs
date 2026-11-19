#!/usr/bin/env bash
# Prints the recipe .github/workflows/release.yml builds one Linux archive by,
# as the `key=value` lines of a step's outputs. That workflow is the one
# statement of the recipe: stage_tree.sh and the pull-request job of
# .github/workflows/docker.yml read it here and state no constant of their own.
#   toolchain               RUST_TOOLCHAIN
#   glibc_floor             GLIBC_FLOOR
#   zig_version             ZIG_VERSION
#   zig_sha256              the ZIG_SHA256_* constant of the row's `zig`
#   cargo_zigbuild_version  CARGO_ZIGBUILD_VERSION
#   target, cpu, zig, openssl, openssl_sha256
#                           the keys of the build matrix row whose `name` is
#                           <archive> (the archive's suffix after `vep-X.Y.Z-`)
# A constant the workflow does not state once, in the shape read here, is
# refused by its name, and a refused run prints no output line.
#
# Usage: scripts/docker/release_recipe.sh <archive>
# WORKFLOW overrides the workflow file read.
set -euo pipefail

usage="usage: release_recipe.sh <archive>"
archive="${1:?$usage}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORKFLOW="${WORKFLOW:-$ROOT/.github/workflows/release.yml}"

fail() {
    echo "ERROR: [release_recipe] $*" >&2
    exit 1
}

[[ -f "$WORKFLOW" ]] || fail "no workflow at $WORKFLOW"

constant() { # constant <NAME> <ERE>: the value of the workflow-level `env` line `  NAME: value`
    local value
    value="$(sed -nE "s/^  $1: \"?([^\" ]+)\"?\$/\1/p" "$WORKFLOW")"
    [[ "$value" =~ $2 ]] || fail "$WORKFLOW does not state $1 once in the shape the recipe reads: '$value'"
    printf '%s\n' "$value"
}

# A matrix row is a `- key: value` line and the `key: value` lines under it.
first='^          - ([a-z_0-9]+): (.+)$'
more='^            ([a-z_0-9]+): (.+)$'
rows=0
row=""
name=""
lines=""
close_row() {
    if [[ "$name" == "$archive" ]]; then
        rows=$((rows + 1))
        row="$lines"
    fi
    name=""
    lines=""
}
while IFS= read -r line; do
    if [[ "$line" =~ $first ]]; then
        close_row
    elif [[ ! "$line" =~ $more ]]; then
        close_row
        continue
    fi
    [[ "${BASH_REMATCH[1]}" != name ]] || name="${BASH_REMATCH[2]}"
    lines+="${BASH_REMATCH[1]} ${BASH_REMATCH[2]}"$'\n'
done <"$WORKFLOW"
close_row
[[ "$rows" -eq 1 ]] || fail "$WORKFLOW builds $rows archives named '$archive', one expected"

key() { # key <name> <ERE>: the value of the row's `name: value` line
    local value
    value="$(sed -nE "s/^$1 (.+)\$/\1/p" <<<"$row")"
    [[ "$value" =~ $2 ]] || fail "the row of '$archive' in $WORKFLOW does not state $1 once in the shape the recipe reads: '$value'"
    printf '%s\n' "$value"
}

semver='^[0-9]+\.[0-9]+\.[0-9]+$'
sha256='^[0-9a-f]{64}$'
toolchain="$(constant RUST_TOOLCHAIN "$semver")"
glibc_floor="$(constant GLIBC_FLOOR '^[0-9]+\.[0-9]+$')"
zig_version="$(constant ZIG_VERSION "$semver")"
cargo_zigbuild_version="$(constant CARGO_ZIGBUILD_VERSION "$semver")"
target="$(key target '^(x86_64|aarch64)-unknown-linux-gnu$')"
cpu="$(key cpu '^[a-z0-9-]+$')"
zig="$(key zig '^(x86_64|aarch64)-linux$')"
openssl="$(key openssl '^linux-[a-z0-9]+/openssl-[0-9][A-Za-z0-9._]*-[a-z0-9_]+\.conda$')"
openssl_sha256="$(key openssl_sha256 "$sha256")"
zig_sha256="$(constant "ZIG_SHA256_$(tr 'a-z-' 'A-Z_' <<<"$zig")" "$sha256")"

printf '%s\n' "toolchain=$toolchain" "glibc_floor=$glibc_floor" "zig_version=$zig_version" "zig_sha256=$zig_sha256" \
    "cargo_zigbuild_version=$cargo_zigbuild_version" "target=$target" "cpu=$cpu" "zig=$zig" \
    "openssl=$openssl" "openssl_sha256=$openssl_sha256"
