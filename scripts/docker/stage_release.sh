#!/usr/bin/env bash
# Stages the Docker build context for one archive of a published release:
# downloads the Linux archive named by <target> (the archive's suffix after
# `vep-X.Y.Z-`: a Linux target triple, with a generic build's CPU suffix when
# it has one, such as `aarch64-unknown-linux-gnu-generic`), checks it against
# the release's SHA256SUMS and its build-provenance attestation from release.yml
# at that tag, then moves the three binaries plus LICENSE and NOTICE into
# <context-dir>. Any failed check exits 1 with nothing at <context-dir>.
#
# Usage: scripts/docker/stage_release.sh <tag> <target> <context-dir>
# Needs an authenticated `gh` (GH_TOKEN in CI). REPO overrides the repository.
set -euo pipefail

usage="usage: stage_release.sh <tag> <target> <context-dir>"
tag="${1:?$usage}"
target="${2:?$usage}"
dest="${3:?$usage}"
REPO="${REPO:-natera-open-source/vep-rs}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

fail() {
    echo "ERROR: [stage_release] $*" >&2
    exit 1
}

[[ "$tag" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || fail "tag '$tag' is not vX.Y.Z"
[[ "$target" =~ ^(x86_64|aarch64)-unknown-linux-gnu(-[a-z0-9]+)?$ ]] || fail "target '$target' is not a Linux archive name"
[[ ! -e "$dest" ]] || fail "$dest already exists"
version="${tag#v}"
archive="vep-$version-$target.tar.gz"

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
    [[ -f "$dl/vep-$version-$target/$bin" ]] || fail "$archive has no $bin"
done
mkdir -p "$dl/context"
for bin in vep vep-cache-builder vep-cache-converter; do
    install -m 0755 "$dl/vep-$version-$target/$bin" "$dl/context/$bin"
done
cp "$ROOT/LICENSE" "$ROOT/NOTICE" "$dl/context/"
mkdir -p "$(dirname "$dest")"
mv "$dl/context" "$dest"
echo "staged $tag ($archive, sha256 ${line%% *}) in $dest"
