#!/usr/bin/env bash
# Stages the Docker build context for one archive from the checked-out tree,
# built as .github/workflows/release.yml builds that archive (the archive's
# suffix after `vep-X.Y.Z-`, as in stage_release.sh). The recipe is read from
# that workflow by release_recipe.sh: zig and cargo-zigbuild at its versions
# and the target's conda-forge OpenSSL are installed under <tools-dir>, each
# download held to the checksum the workflow records, and the three binaries
# are linked against its glibc floor with the archive's `-C target-cpu`, the
# build clock at the commit's time as the workflow sets it. The
# release's ELF check (scripts/release/check_glibc_floor.py) holds them to
# that floor, so they do not depend on the build host's glibc, and they move
# with LICENSE and NOTICE into <context-dir>. Any failed step exits 1 with
# nothing at <context-dir>.
#
# Usage: scripts/docker/stage_tree.sh <archive> <context-dir> <tools-dir>
# Needs the recipe's Rust toolchain, git, curl, unzip, zstd, python3 and a
# host of the archive's architecture.
set -euo pipefail

usage="usage: stage_tree.sh <archive> <context-dir> <tools-dir>"
archive="${1:?$usage}"
dest="${2:?$usage}"
tools="${3:?$usage}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

fail() {
    echo "ERROR: [stage_tree] $*" >&2
    exit 1
}

recipe="$(bash "$HERE/release_recipe.sh" "$archive")" || fail "release.yml states no build recipe for '$archive'"
value() { sed -n "s/^$1=//p" <<<"$recipe"; }
GLIBC_FLOOR="$(value glibc_floor)"
ZIG_VERSION="$(value zig_version)"
ZIG_SHA256="$(value zig_sha256)"
CARGO_ZIGBUILD_VERSION="$(value cargo_zigbuild_version)"
TARGET="$(value target)"
CPU="$(value cpu)"
ZIG_ARCH="$(value zig)"
OPENSSL_PACKAGE="$(value openssl)"
OPENSSL_SHA256="$(value openssl_sha256)"

[[ ! -e "$dest" ]] || fail "$dest already exists"
[[ ! -e "$tools" ]] || fail "$tools already exists"
host="$(uname -m)"
[[ "$host" == "${TARGET%%-*}" ]] || fail "$archive needs a host of architecture ${TARGET%%-*}, not $host"

# The build clock is the time of the tree's own commit: a tree unpacked inside
# another checkout would otherwise read that checkout's. Under
# log.showSignature git prints a signed commit's signature ahead of its time.
clock="$(git -C "$ROOT" -c log.showSignature=false log -1 --format=%ct 2>/dev/null)" || clock=""
[[ "$(git -C "$ROOT" rev-parse --show-toplevel 2>/dev/null)" -ef "$ROOT" && "$clock" =~ ^[0-9]+$ ]] ||
    fail "the tree at $ROOT has no commit of its own to take the build clock from"

mkdir -p "$tools/zig" "$tools/openssl"
tools="$(cd "$tools" && pwd)"

# zig links the binaries against glibc stubs of the floor.
curl -fsSLo "$tools/zig.tar.xz" "https://ziglang.org/download/$ZIG_VERSION/zig-$ZIG_ARCH-$ZIG_VERSION.tar.xz" ||
    fail "could not download zig $ZIG_VERSION for $ZIG_ARCH"
echo "$ZIG_SHA256  $tools/zig.tar.xz" | sha256sum --check --strict --status - ||
    fail "the zig $ZIG_VERSION tarball does not match the checksum release.yml records"
tar -xJf "$tools/zig.tar.xz" -C "$tools/zig" --strip-components=1
(cd "$ROOT" && cargo install cargo-zigbuild --locked --version "$CARGO_ZIGBUILD_VERSION" --root "$tools/cargo") ||
    fail "cargo could not install cargo-zigbuild $CARGO_ZIGBUILD_VERSION"

# Only the link needs the package: the image's own libssl.so.3 loads at run
# time. openssl-sys reads the target-prefixed OPENSSL_DIR.
curl -fsSLo "$tools/openssl.conda" "https://conda.anaconda.org/conda-forge/$OPENSSL_PACKAGE" ||
    fail "could not download $OPENSSL_PACKAGE"
echo "$OPENSSL_SHA256  $tools/openssl.conda" | sha256sum --check --strict --status - ||
    fail "$OPENSSL_PACKAGE does not match the checksum release.yml records"
unzip -q "$tools/openssl.conda" 'pkg-*.tar.zst' -d "$tools"
tar --zstd -xf "$tools"/pkg-*.tar.zst -C "$tools/openssl" lib include
openssl_dir="$(echo "$TARGET" | tr 'a-z-' 'A-Z_')_OPENSSL_DIR=$tools/openssl"

# RUSTFLAGS replaces the target rustflags of .cargo/config.toml, as in
# release.yml. The installed cargo-zigbuild is run by its path: `cargo
# zigbuild` runs the first one cargo finds, its own home's before PATH's.
(cd "$ROOT" && env PATH="$tools/zig:$PATH" "$openssl_dir" RUSTFLAGS="-C target-cpu=$CPU" SOURCE_DATE_EPOCH="$clock" \
    "$tools/cargo/bin/cargo-zigbuild" zigbuild --release --locked --target "$TARGET.$GLIBC_FLOOR" \
    --bin vep --bin vep-cache-builder --bin vep-cache-converter) ||
    fail "cargo could not build the three binaries for $archive"

# cargo resolves a relative CARGO_TARGET_DIR from the directory it ran in.
out="${CARGO_TARGET_DIR:-target}"
[[ "$out" == /* ]] || out="$ROOT/$out"
out="$out/$TARGET/release"
for bin in vep vep-cache-builder vep-cache-converter; do
    [[ -f "$out/$bin" ]] || fail "the build left no $bin in $out"
done
python3 "$ROOT/scripts/release/check_glibc_floor.py" --floor "$GLIBC_FLOOR" \
    "$out/vep" "$out/vep-cache-builder" "$out/vep-cache-converter" ||
    fail "a binary built for $archive needs a glibc newer than the floor $GLIBC_FLOOR"

stage="$(mktemp -d)"
trap 'rm -rf "$stage"' EXIT
mkdir "$stage/context"
for bin in vep vep-cache-builder vep-cache-converter; do
    install -m 0755 "$out/$bin" "$stage/context/$bin"
done
cp "$ROOT/LICENSE" "$ROOT/NOTICE" "$stage/context/"
mkdir -p "$(dirname "$dest")"
mv "$stage/context" "$dest"
echo "staged the tree's build for $archive (target-cpu $CPU, glibc floor $GLIBC_FLOOR) in $dest"
