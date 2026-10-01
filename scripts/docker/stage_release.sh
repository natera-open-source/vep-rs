#!/usr/bin/env bash
# Stages the Docker build context for a published release: downloads its
# x86_64 Linux archive, checks it against the release's SHA256SUMS and its
# build-provenance attestation from release.yml at that tag, then moves the
# three binaries plus LICENSE and NOTICE into <context-dir>. Any failed check
# exits 1 with nothing at <context-dir>.
#
# Usage: scripts/docker/stage_release.sh <tag> <context-dir>
# Needs an authenticated `gh` (GH_TOKEN in CI). REPO overrides the repository.
set -euo pipefail

usage="usage: stage_release.sh <tag> <context-dir>"
tag="${1:?$usage}"
dest="${2:?$usage}"
REPO="${REPO:-natera-open-source/vep-rs}"
TARGET=x86_64-unknown-linux-gnu
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

fail() {
    echo "ERROR: [stage_release] $*" >&2
    exit 1
}

[[ "$tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || fail "tag '$tag' is not vX.Y.Z"
[[ ! -e "$dest" ]] || fail "$dest already exists"
version="${tag#v}"
archive="vep-$version-$TARGET.tar.gz"

dl="$(mktemp -d)"
trap 'rm -rf "$dl"' EXIT
gh release download "$tag" --repo "$REPO" --dir "$dl" --pattern "$archive" --pattern SHA256SUMS ||
    fail "could not download $archive and SHA256SUMS from release $tag"

line="$(awk -v f="$archive" '$2 == f' "$dl/SHA256SUMS")"
[[ -n "$line" ]] || fail "SHA256SUMS has no line for $archive"
(cd "$dl" && printf '%s\n' "$line" | sha256sum --check --strict --status -) ||
    fail "$archive does not match SHA256SUMS"
gh attestation verify "$dl/$archive" --repo "$REPO" \
    --signer-workflow "$REPO/.github/workflows/release.yml" \
    --source-ref "refs/tags/$tag" >/dev/null ||
    fail "$archive has no attestation from $REPO's release workflow at $tag"

tar -xzf "$dl/$archive" -C "$dl"
for bin in vep vep-cache-builder vep-cache-converter; do
    [[ -f "$dl/vep-$version-$TARGET/$bin" ]] || fail "$archive has no $bin"
done
mkdir -p "$dl/context"
for bin in vep vep-cache-builder vep-cache-converter; do
    install -m 0755 "$dl/vep-$version-$TARGET/$bin" "$dl/context/$bin"
done
cp "$ROOT/LICENSE" "$ROOT/NOTICE" "$dl/context/"
mkdir -p "$(dirname "$dest")"
mv "$dl/context" "$dest"
echo "staged $tag ($archive, sha256 ${line%% *}) in $dest"
