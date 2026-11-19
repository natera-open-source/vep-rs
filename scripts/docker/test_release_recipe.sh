#!/usr/bin/env bash
# Tests for release_recipe.sh: the recipe of a row of a fabricated release
# workflow is read whole, a constant the workflow does not state once in its
# shape is refused by its name, and on the tree's own workflows every constant
# read is the one .github/workflows/release.yml states while
# .github/workflows/docker.yml states none of them.
#
# Pure bash. No network, no GitHub.
# Run: bash scripts/docker/test_release_recipe.sh
# shellcheck disable=SC2016  # the workflow's text is matched as written, unexpanded
set -uo pipefail

PASS=0
FAIL=0
ok() {
    PASS=$((PASS + 1))
    echo "  ok   $1"
}
bad() {
    FAIL=$((FAIL + 1))
    echo "  FAIL $1"
}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RECIPE="$HERE/release_recipe.sh"
ROOT=$(cd "$HERE/../.." && pwd)
RELEASE_YML="$ROOT/.github/workflows/release.yml"
DOCKER_YML="$ROOT/.github/workflows/docker.yml"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

A=$(printf 'a%.0s' {1..64})
B=$(printf 'b%.0s' {1..64})
C=$(printf 'c%.0s' {1..64})
D=$(printf 'd%.0s' {1..64})
cat >"$TMP/release.yml" <<EOF
name: Release

env:
  RUST_TOOLCHAIN: 1.2.3
  # The oldest glibc the archives run on.
  GLIBC_FLOOR: "2.99"
  ZIG_VERSION: 0.1.2
  ZIG_SHA256_X86_64_LINUX: $A
  ZIG_SHA256_AARCH64_LINUX: $B
  CARGO_ZIGBUILD_VERSION: 4.5.6

jobs:
  build:
    strategy:
      matrix:
        include:
          - target: x86_64-unknown-linux-gnu
            name: x86_64-unknown-linux-gnu-v2
            cpu: x86-64-v2
            os: ubuntu-latest
            zig: x86_64-linux
            openssl: linux-64/openssl-9.9.9-habc_0.conda
            openssl_sha256: $C
          - target: aarch64-unknown-linux-gnu
            name: aarch64-unknown-linux-gnu-generic
            cpu: generic
            os: ubuntu-24.04-arm
            zig: aarch64-linux
            openssl: linux-aarch64/openssl-9.9.9-hdef_1.conda
            openssl_sha256: $D
          - target: aarch64-apple-darwin
            name: aarch64-apple-darwin
            os: macos-latest
    steps:
      - uses: actions/checkout@v1
        with:
          toolchain: \${{ env.RUST_TOOLCHAIN }}
EOF

run() { # run <workflow> <archive>
    env -i PATH="$PATH" WORKFLOW="$1" bash "$RECIPE" "$2"
}

edited() { # edited <sed expression>: a copy of the fabricated workflow with one edit
    sed -E "$1" "$TMP/release.yml" >"$TMP/edited.yml"
    echo "$TMP/edited.yml"
}

expect_fail() { # expect_fail <name> <stderr substring> <workflow> <archive>: exits non-zero with no output line
    local name="$1" want="$2" out err rc
    out=$(run "$3" "$4" 2>"$TMP/stderr")
    rc=$?
    err=$(<"$TMP/stderr")
    if [[ "$rc" -ne 0 && "$err" == "ERROR: [release_recipe] "*"$want"* && -z "$out" ]]; then ok "$name"; else bad "$name (rc=$rc, stdout: $out, stderr: $err)"; fi
}

bash -n "$RECIPE" 2>/dev/null && ok "release_recipe.sh parses" || bad "release_recipe.sh is missing or has a syntax error"

want="toolchain=1.2.3
glibc_floor=2.99
zig_version=0.1.2
zig_sha256=$A
cargo_zigbuild_version=4.5.6
target=x86_64-unknown-linux-gnu
cpu=x86-64-v2
zig=x86_64-linux
openssl=linux-64/openssl-9.9.9-habc_0.conda
openssl_sha256=$C"
got=$(run "$TMP/release.yml" x86_64-unknown-linux-gnu-v2 2>&1)
[[ "$got" == "$want" ]] && ok "the recipe of a row is its ten constants, the workflow's and the row's" ||
    bad "the x86_64 row's recipe: got [$got]"
got=$(run "$TMP/release.yml" aarch64-unknown-linux-gnu-generic 2>&1)
[[ "$got" == *"zig_sha256=$B"* && "$got" == *"openssl_sha256=$D"* && "$got" == *"cpu=generic"* ]] &&
    ok "the zig checksum is the one of the row's zig, the rest the row's own" ||
    bad "the arm64 row's recipe: got [$got]"

expect_fail "an archive the workflow does not build is refused" "builds 0 archives named 'x86_64-unknown-linux-gnu-v4'" \
    "$TMP/release.yml" x86_64-unknown-linux-gnu-v4
expect_fail "an archive the workflow builds twice is refused" "builds 2 archives named 'aarch64-unknown-linux-gnu-generic'" \
    "$(edited 's/^(            name: )x86_64-unknown-linux-gnu-v2$/\1aarch64-unknown-linux-gnu-generic/')" aarch64-unknown-linux-gnu-generic
expect_fail "a macOS archive, which has no zig recipe, is refused" "does not state target once in the shape the recipe reads: 'aarch64-apple-darwin'" \
    "$TMP/release.yml" aarch64-apple-darwin
expect_fail "a missing workflow is refused" "no workflow at $TMP/absent.yml" "$TMP/absent.yml" x86_64-unknown-linux-gnu-v2

# Each constant, taken out of the workflow or out of its shape, is refused by its name.
refused=0
wrong=""
refuse() { # refuse <name the error carries> <sed expression> [archive]
    local out rc
    out=$(run "$(edited "$2")" "${3:-aarch64-unknown-linux-gnu-generic}" 2>&1)
    rc=$?
    if [[ "$rc" -ne 0 && "$out" == "ERROR: [release_recipe] "*"does not state $1 once in the shape the recipe reads"* ]]; then
        refused=$((refused + 1))
    else
        wrong+=" [$1: $2]"
    fi
}
refuse RUST_TOOLCHAIN '/^  RUST_TOOLCHAIN: /d'
refuse RUST_TOOLCHAIN 's/^(  RUST_TOOLCHAIN: ).*/\1stable/'
refuse GLIBC_FLOOR '/^  GLIBC_FLOOR: /d'
refuse GLIBC_FLOOR 's/^(  GLIBC_FLOOR: .*)$/\1\n  GLIBC_FLOOR: "2.98"/'
refuse ZIG_VERSION 's/^(  ZIG_VERSION: ).*/\1master/'
refuse CARGO_ZIGBUILD_VERSION '/^  CARGO_ZIGBUILD_VERSION: /d'
refuse ZIG_SHA256_AARCH64_LINUX '/^  ZIG_SHA256_AARCH64_LINUX: /d'
refuse ZIG_SHA256_X86_64_LINUX 's/^(  ZIG_SHA256_X86_64_LINUX: ).*/\1abc/' x86_64-unknown-linux-gnu-v2
refuse cpu '/^            cpu: generic$/d'
refuse zig 's/^(            zig: )aarch64-linux$/\1aarch64-macos/'
refuse openssl 's/^(            openssl: )linux-aarch64.*/\1openssl/'
refuse openssl_sha256 "s/^(            openssl_sha256: )$D\$/\\1${D:0:63}/"
[[ "$refused" -eq 12 && -z "$wrong" ]] && ok "a constant the workflow lacks, repeats or misshapes is refused by its name" ||
    bad "$refused of 12 edits refused by the constant's name; not refused:$wrong"

# The tree's own workflows. release.yml is read here a second way, line by
# line, and each constant the recipe answers must be the one it states.
archives=$(sed -nE 's/^  ARCHIVES: (.+)$/\1/p' "$DOCKER_YML")
rows=$(awk '
    /^          - target: / { target = $3; name = ""; cpu = ""; zig = ""; openssl = "" }
    /^            name: / { name = $2 }
    /^            cpu: / { cpu = $2 }
    /^            zig: / { zig = $2 }
    /^            openssl: / { openssl = $2 }
    /^            openssl_sha256: / && name != "" { print name, target, cpu, zig, openssl, $2 }' "$RELEASE_YML")
stated() { sed -nE "s/^  $1: \"?([^\" ]+)\"?\$/\1/p" "$RELEASE_YML"; }
answer() { sed -n "s/^$1=//p" <<<"$2"; }
count=0
differs=""
for name in $archives; do
    count=$((count + 1))
    got=$(env -i PATH="$PATH" bash "$RECIPE" "$name" 2>/dev/null) || { differs+=" $name:unread"; continue; }
    read -r _ target cpu zig openssl openssl_sha256 <<<"$(grep "^$name " <<<"$rows")"
    [[ "$(answer toolchain "$got")" == "$(stated RUST_TOOLCHAIN)" ]] || differs+=" $name:RUST_TOOLCHAIN"
    [[ "$(answer glibc_floor "$got")" == "$(stated GLIBC_FLOOR)" ]] || differs+=" $name:GLIBC_FLOOR"
    [[ "$(answer zig_version "$got")" == "$(stated ZIG_VERSION)" ]] || differs+=" $name:ZIG_VERSION"
    [[ "$(answer cargo_zigbuild_version "$got")" == "$(stated CARGO_ZIGBUILD_VERSION)" ]] || differs+=" $name:CARGO_ZIGBUILD_VERSION"
    case "$zig" in
        x86_64-linux) sum=$(stated ZIG_SHA256_X86_64_LINUX) ;;
        aarch64-linux) sum=$(stated ZIG_SHA256_AARCH64_LINUX) ;;
        *) sum="" ;;
    esac
    [[ -n "$sum" && "$(answer zig_sha256 "$got")" == "$sum" ]] || differs+=" $name:ZIG_SHA256"
    [[ -n "$target" && "$(answer target "$got")" == "$target" ]] || differs+=" $name:target"
    [[ -n "$cpu" && "$(answer cpu "$got")" == "$cpu" ]] || differs+=" $name:cpu"
    [[ -n "$zig" && "$(answer zig "$got")" == "$zig" ]] || differs+=" $name:zig"
    [[ -n "$openssl" && "$(answer openssl "$got")" == "$openssl" ]] || differs+=" $name:openssl"
    [[ -n "$openssl_sha256" && "$(answer openssl_sha256 "$got")" == "$openssl_sha256" ]] || differs+=" $name:openssl_sha256"
done
[[ "$count" -eq 4 && -z "$differs" ]] && ok "every archive docker.yml packages has the recipe release.yml states for it" ||
    bad "$count archives read from docker.yml, 4 expected; the recipe differs from release.yml in:$differs"

# docker.yml holds no second statement of a recipe constant: its toolchain
# step and its build take them from the recipe.
second=""
for constant in RUST_TOOLCHAIN GLIBC_FLOOR ZIG_VERSION ZIG_SHA256_X86_64_LINUX ZIG_SHA256_AARCH64_LINUX CARGO_ZIGBUILD_VERSION; do
    value=$(stated "$constant")
    [[ -n "$value" ]] || { second+=" $constant:unread"; continue; }
    ! grep -qF -e "$constant" -e "$value" "$DOCKER_YML" || second+=" $constant"
done
while read -r name _ _ _ openssl openssl_sha256; do
    [[ -n "$name" ]] || continue
    ! grep -qF -e "$openssl" -e "$openssl_sha256" "$DOCKER_YML" || second+=" $name:openssl"
done <<<"$rows"
[[ -n "$rows" && -z "$second" ]] && ok "docker.yml states no constant of the release recipe" ||
    bad "docker.yml states a constant release.yml states:$second"
grep -qxF '        run: scripts/docker/release_recipe.sh "${{ matrix.archive }}" >> "$GITHUB_OUTPUT"' "$DOCKER_YML" &&
    grep -qxF '          toolchain: ${{ steps.recipe.outputs.toolchain }}' "$DOCKER_YML" &&
    grep -qxF '          targets: ${{ steps.recipe.outputs.target }}' "$DOCKER_YML" &&
    [[ "$(grep -c '^          toolchain: ' "$DOCKER_YML")" -eq 1 ]] &&
    ok "docker.yml installs the toolchain and target the recipe names" ||
    bad "docker.yml's toolchain step does not take its toolchain and target from release_recipe.sh"

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
