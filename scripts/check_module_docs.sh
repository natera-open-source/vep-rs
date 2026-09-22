#!/usr/bin/env bash
#
# Verify every `crates/**/*.rs` file's module documentation survived the license
# header.
#
# Each source file carries a two-line Apache header. When the header's last line
# and the module doc share one physical line, the result is:
#
#     // Licensed under the Apache License, Version 2.0 (see LICENSE).//! Module doc.
#
# which is a valid line comment, so the file compiles, `cargo fmt --check`
# passes, and every test stays green -- while the module's first rustdoc line is
# silently demoted to an ordinary comment and vanishes from the docs. Nothing
# else in CI can see it, hence this gate.
#
# Checks, per file:
#   1. No `LICENSE).//` fusion anywhere.
#   2. When a license header is present, it is followed by a blank line.
#
# Usage: scripts/check_module_docs.sh [crates_dir]
# Exits non-zero on the first category of failure, listing every offending file.

set -euo pipefail

CRATES_DIR="${1:-crates}"

if [[ ! -d $CRATES_DIR ]]; then
    echo "ERROR: [check_module_docs] no such directory: $CRATES_DIR" >&2
    exit 2
fi

status=0

# 1. Fused header + doc comment on one physical line.
fused=$(grep -rn 'LICENSE)\.//' "$CRATES_DIR" --include='*.rs' || true)
if [[ -n $fused ]]; then
    echo "ERROR: [check_module_docs] license header fused to the following line:" >&2
    echo "$fused" >&2
    echo "Fix: end the license header with a newline so the module doc starts its own line." >&2
    status=1
fi

# 2. A license header must be followed by a blank line.
while IFS= read -r file; do
    if [[ $(sed -n '1p' "$file") == "// Copyright"* ]]; then
        third=$(sed -n '3p' "$file")
        if [[ -n $third ]]; then
            echo "ERROR: [check_module_docs] $file: line 3 must be blank after the 2-line license header, got: $third" >&2
            status=1
        fi
    fi
done < <(find "$CRATES_DIR" -name '*.rs' -type f)

if [[ $status -eq 0 ]]; then
    echo "check_module_docs: OK"
fi

exit "$status"
