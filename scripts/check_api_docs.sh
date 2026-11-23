#!/usr/bin/env bash
#
# Build the API documentation with rustdoc's warnings denied.
#
# rustdoc reports a doc comment it cannot render as written as a warning: a link
# to an item that does not exist or that the built page cannot reach, a bare
# URL, an angle-bracketed word read as an HTML tag. `cargo doc` exits 0 on a
# warning, and the test, clippy and fmt commands never read one, so this gate
# builds with `RUSTDOCFLAGS="-D warnings"`, under which each is an error.
#
# Usage: scripts/check_api_docs.sh workspace|tools    (from the repository root)
#   workspace  the root workspace's crates, with the feature that builds the
#              `bench_plugins` binary
#   tools      every Cargo root under tools/, which no `--workspace` command
#              reaches
#
# Prints rustdoc's diagnostics, then one
# `ERROR: [api-docs] <crate>: <n> rustdoc diagnostics` line per failing crate,
# and exits 1. Exits 2 on a usage error.

set -euo pipefail

usage() {
    echo "Usage: scripts/check_api_docs.sh workspace|tools" >&2
    exit 2
}

[[ $# -eq 1 ]] || usage

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

status=0

# cargo reads CARGO_ENCODED_RUSTDOCFLAGS, its flags divided by the unit separator, in place
# of RUSTDOCFLAGS when the environment sets it, so that variable takes the denial as well.
if [[ -n ${CARGO_ENCODED_RUSTDOCFLAGS+set} ]]; then
    sep=$'\037'
    encoded=${CARGO_ENCODED_RUSTDOCFLAGS:+$CARGO_ENCODED_RUSTDOCFLAGS$sep}
    export CARGO_ENCODED_RUSTDOCFLAGS="${encoded}-D${sep}warnings"
fi

# document <manifest> [cargo doc arguments]
#
# `--keep-going` documents every crate after the first failure, so one run names
# every failing crate. The diagnostics of crates documented in parallel
# interleave in cargo's rendered output, so the count per crate is taken from
# its JSON messages, each of which carries its package.
document() {
    local manifest=$1
    shift
    local rc=0
    RUSTDOCFLAGS="${RUSTDOCFLAGS:+$RUSTDOCFLAGS }-D warnings" cargo doc \
        --manifest-path "$manifest" --no-deps --locked --keep-going \
        --message-format=json "$@" > "$scratch/messages.jsonl" || rc=$?
    if ! cargo metadata --manifest-path "$manifest" --no-deps --format-version 1 \
        > "$scratch/metadata.json"; then
        echo "ERROR: [api-docs] $manifest: cargo metadata failed" >&2
        status=1
        return
    fi
    if ! python3 - "$scratch/messages.jsonl" "$scratch/metadata.json" <<'PYEOF'; then
import json
import sys

messages, metadata = sys.argv[1], sys.argv[2]
with open(metadata, encoding="utf-8") as handle:
    names = {package["id"]: package["name"] for package in json.load(handle)["packages"]}
counts = {}
with open(messages, encoding="utf-8") as handle:
    for line in handle:
        if not line.startswith("{"):
            continue
        event = json.loads(line)
        if event.get("reason") != "compiler-message":
            continue
        diagnostic = event["message"]
        sys.stderr.write(diagnostic.get("rendered") or "")
        # rustc closes a failed compilation with a summary that is no diagnostic of its own.
        if diagnostic.get("level") != "error" or diagnostic.get("message", "").startswith("aborting due to"):
            continue
        name = names.get(event["package_id"], event["package_id"])
        counts[name] = counts.get(name, 0) + 1
for name in sorted(counts):
    print(f"ERROR: [api-docs] {name}: {counts[name]} rustdoc diagnostics", file=sys.stderr)
sys.exit(1 if counts else 0)
PYEOF
        status=1
    elif [[ $rc -ne 0 ]]; then
        echo "ERROR: [api-docs] $manifest: cargo doc exited $rc with no rustdoc diagnostic" >&2
        status=1
    fi
}

case $1 in
    workspace)
        if [[ ! -f Cargo.toml ]]; then
            echo "ERROR: [api-docs] no Cargo.toml here: run from the repository root" >&2
            exit 2
        fi
        document Cargo.toml --workspace --features vep-cli/bench
        ;;
    tools)
        shopt -s nullglob
        manifests=(tools/*/Cargo.toml)
        if [[ ${#manifests[@]} -eq 0 ]]; then
            echo "ERROR: [api-docs] no Cargo root under tools/: run from the repository root" >&2
            exit 2
        fi
        for manifest in "${manifests[@]}"; do
            document "$manifest"
        done
        ;;
    *)
        usage
        ;;
esac

if [[ $status -eq 0 ]]; then
    echo "check_api_docs: OK"
fi

exit "$status"
