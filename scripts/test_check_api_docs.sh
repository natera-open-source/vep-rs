#!/usr/bin/env bash
# Tests for check_api_docs.sh: a crate whose doc comments rustdoc warns about
# is refused by name with its count of diagnostics, the same crate with
# well-formed comments passes, and the arguments and trees the gate refuses.
#
# Every scratch crate is free of dependencies, so the tests need cargo, rustdoc
# and python3 and no network. Run: bash scripts/test_check_api_docs.sh
#
# The scratch sources and the expected lines hold backticks as text.
# shellcheck disable=SC2016
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
GATE="$HERE/check_api_docs.sh"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
# cargo writes every scratch build under $WORK, whatever target directory the environment or a
# cargo config names.
export CARGO_TARGET_DIR="$WORK/target"

crate() { # crate <dir> <package name> <lib.rs text> [further manifest text]
    mkdir -p "$1/src"
    printf '[package]\nname = "%s"\nversion = "0.0.0"\nedition = "2021"\n%s' "$2" "${4:-}" > "$1/Cargo.toml"
    printf '%s\n' "$3" > "$1/src/lib.rs"
}
lock() { # lock <manifest>: the lockfile `--locked` needs
    cargo generate-lockfile --offline --manifest-path "$1"
}
run() { # run <tree> [argument...]: the gate's output and status, run from the tree's root
    local tree=$1
    shift
    OUT=$(cd "$tree" && bash "$GATE" "$@" 2>&1)
    RC=$?
}
expect() { # expect <name> <status> [substring the output holds]...
    local name=$1 want=$2 text
    shift 2
    if [[ "$RC" -ne "$want" ]]; then
        bad "$name (rc=$RC, want $want): $OUT"
        return
    fi
    for text in "$@"; do
        if [[ "$OUT" != *"$text"* ]]; then
            bad "$name: the output lacks [$text]: $OUT"
            return
        fi
    done
    ok "$name"
}
refuse() { # refuse <name> <substring the output must not hold>
    if [[ "$OUT" != *"$2"* ]]; then
        ok "$1"
    else
        bad "$1: the output holds [$2]: $OUT"
    fi
}

CLEAN='//! A crate whose documentation is well-formed.

/// Does nothing.
pub fn f() {}'

if bash -n "$GATE"; then
    ok "check_api_docs.sh parses"
else
    bad "check_api_docs.sh has a syntax error"
fi

# One Cargo root under tools/ with one broken intra-doc link.
T="$WORK/one-tool"
crate "$T/tools/linked" linked '//! A crate whose documentation names [`NoSuchItem`].'
lock "$T/tools/linked/Cargo.toml"
run "$T" tools
expect "a tool with a broken intra-doc link is refused by name with its count" 1 \
    'ERROR: [api-docs] linked: 1 rustdoc diagnostics'
expect "rustdoc's own diagnostic is printed with the error line" 1 'unresolved link to `NoSuchItem`'
refuse "a refused run does not print OK" 'check_api_docs: OK'
crate "$T/tools/linked" linked '//! A crate whose documentation names `NoSuchItem` in a code span.'
run "$T" tools
expect "the same tool with the link fixed passes" 0 'check_api_docs: OK'
refuse "a passing run prints no error line" 'ERROR: [api-docs]'

# cargo reads CARGO_ENCODED_RUSTDOCFLAGS in place of RUSTDOCFLAGS when the environment sets it.
T="$WORK/encoded"
crate "$T/tools/linked" linked '//! A crate whose documentation names [`NoSuchItem`].'
lock "$T/tools/linked/Cargo.toml"
OUT=$(cd "$T" && CARGO_ENCODED_RUSTDOCFLAGS=$'--cfg\037encoded' bash "$GATE" tools 2>&1)
RC=$?
expect "a tool is refused when the environment sets cargo's encoded rustdoc flags" 1 \
    'ERROR: [api-docs] linked: 1 rustdoc diagnostics' 'unresolved link to `NoSuchItem`'

# Two roots under tools/: a failure in the first does not hide the second.
T="$WORK/two-tools"
crate "$T/tools/first" first '//! Fetches from https://example.org when asked.'
crate "$T/tools/second" second '//! One record per line as <name>.

/// Reads as [`hidden`] does.
pub fn shown() {}

fn hidden() {}'
crate "$T/tools/third" third "$CLEAN"
for root in first second third; do lock "$T/tools/$root/Cargo.toml"; done
run "$T" tools
expect "every failing tool root is named, each with its own count" 1 \
    'ERROR: [api-docs] first: 1 rustdoc diagnostics' 'ERROR: [api-docs] second: 2 rustdoc diagnostics'
refuse "a clean tool root beside failing ones is not named" 'ERROR: [api-docs] third'

# A workspace: the gate's workspace form names the `bench` feature of the
# `vep-cli` package, so the scratch workspace has a member of that name with it.
T="$WORK/workspace"
mkdir -p "$T"
printf '[workspace]\nmembers = ["crates/vep-cli", "crates/other", "crates/clean"]\nresolver = "2"\n' > "$T/Cargo.toml"
crate "$T/crates/vep-cli" vep-cli '//! A crate with two doc comments rustdoc warns about.

/// Reads as [`hidden`] does.
pub fn shown() {}

/// Fetches from https://example.org when asked.
pub fn fetch() {}

fn hidden() {}' $'\n[features]\nbench = []\n'
crate "$T/crates/other" other '//! One record per line as <name>.'
crate "$T/crates/clean" clean "$CLEAN"
lock "$T/Cargo.toml"
# One job at a time: cargo then starts no crate after the first failure unless told to keep going.
OUT=$(cd "$T" && CARGO_BUILD_JOBS=1 bash "$GATE" workspace 2>&1)
RC=$?
expect "each failing workspace crate is named with its own count" 1 \
    'ERROR: [api-docs] vep-cli: 2 rustdoc diagnostics' 'ERROR: [api-docs] other: 1 rustdoc diagnostics'
refuse "a clean workspace crate beside failing ones is not named" 'ERROR: [api-docs] clean'
crate "$T/crates/vep-cli" vep-cli '//! A crate whose doc comments are well-formed.

/// Reads as `hidden` does.
pub fn shown() {}

/// Fetches from <https://example.org> when asked.
pub fn fetch() {}

fn hidden() {}' $'\n[features]\nbench = []\n'
crate "$T/crates/other" other '//! One record per line as `<name>`.'
run "$T" workspace
expect "the same workspace with the comments fixed passes" 0 'check_api_docs: OK'

# A crate another one depends on is compiled by rustc, whose warnings are not rustdoc's.
T="$WORK/checked"
mkdir -p "$T"
printf '[workspace]\nmembers = ["crates/vep-cli", "crates/lower"]\nresolver = "2"\n' > "$T/Cargo.toml"
crate "$T/crates/lower" lower '//! A crate whose documentation is well-formed.

fn unused() {}'
crate "$T/crates/vep-cli" vep-cli "$CLEAN" $'\n[features]\nbench = []\n\n[dependencies]\nlower = { path = "../lower" }\n'
lock "$T/Cargo.toml"
run "$T" workspace
expect "a rustc warning on a crate compiled for another is printed and the run passes" 0 \
    'is never used' 'check_api_docs: OK'
refuse "a rustc warning is not counted as a rustdoc diagnostic" 'ERROR: [api-docs]'

# A build that fails before rustdoc runs is also refused: no lockfile under --locked.
T="$WORK/unlocked"
crate "$T/tools/unlocked" unlocked "$CLEAN"
run "$T" tools
expect "a build that stops before rustdoc is refused" 1 \
    'ERROR: [api-docs] tools/unlocked/Cargo.toml: cargo doc exited'

# Arguments and trees the gate does not read.
run "$WORK/one-tool" docs
expect "an unknown argument exits 2 with the usage line" 2 'Usage: scripts/check_api_docs.sh workspace|tools'
run "$WORK/one-tool"
expect "no argument exits 2" 2 'Usage:'
run "$WORK/one-tool" workspace tools
expect "two arguments exit 2" 2 'Usage:'
run "$WORK/one-tool" workspace
expect "the workspace form outside a Cargo root exits 2" 2 'ERROR: [api-docs] no Cargo.toml here'
run "$WORK/workspace" tools
expect "the tools form in a tree without tools/ exits 2" 2 'ERROR: [api-docs] no Cargo root under tools/'

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
