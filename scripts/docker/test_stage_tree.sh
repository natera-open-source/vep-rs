#!/usr/bin/env bash
# Tests for stage_tree.sh: in a fabricated tree whose release workflow states
# a recipe, stub `curl`, `tar`, `unzip`, `cargo`, the `cargo-zigbuild` it
# installs and `python3` record the downloads, the install, the build and the
# ELF check the script asks for and leave three fabricated binaries, the
# archive argument picks the row whose
# recipe is followed, the build's clock is the time of the tree's one commit,
# and each fail-closed gate (archive, existing
# destination or tools directory, host architecture, a tree with no commit of
# its own, failed download,
# checksum, failed install, failed build, missing binary, glibc floor) is
# shown to stop the run with nothing staged. The checksums are
# checked by the real `sha256sum` against the fabricated downloads.
#
# Bash and git. No compiler, no network. The script runs from a fabricated
# tree, so nothing is written into this one.
# Run: bash scripts/docker/test_stage_tree.sh
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
ROOT=$(cd "$HERE/../.." && pwd)
RELEASE_YML="$ROOT/.github/workflows/release.yml"
DOCKER_YML="$ROOT/.github/workflows/docker.yml"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

TREE="$TMP/tree"
STAGE="$TREE/scripts/docker/stage_tree.sh"
STUBS="$TMP/bin"
LOG="$TMP/calls.log"
mkdir -p "$TREE/scripts/docker" "$TREE/scripts/release" "$TREE/.github/workflows" "$STUBS" "$TMP/sums"
cp "$HERE/stage_tree.sh" "$HERE/release_recipe.sh" "$TREE/scripts/docker/" 2>/dev/null
echo license >"$TREE/LICENSE"
echo notice >"$TREE/NOTICE"
: >"$TREE/scripts/release/check_glibc_floor.py"

# The downloads the stub curl hands out, and their checksums for the
# fabricated workflow.
echo "zig tarball" >"$TMP/sums/zig"
echo "openssl package" >"$TMP/sums/openssl"
ZIG_SUM=$(sha256sum <"$TMP/sums/zig" | cut -d' ' -f1)
OPENSSL_SUM=$(sha256sum <"$TMP/sums/openssl" | cut -d' ' -f1)
{
    printf '%s\n' 'env:' '  RUST_TOOLCHAIN: 1.2.3' '  GLIBC_FLOOR: "2.99"' '  ZIG_VERSION: 0.1.2' \
        "  ZIG_SHA256_X86_64_LINUX: $ZIG_SUM" "  ZIG_SHA256_AARCH64_LINUX: $ZIG_SUM" '  CARGO_ZIGBUILD_VERSION: 4.5.6' \
        'jobs:' '  build:' '    strategy:' '      matrix:' '        include:'
    while read -r name target cpu zig package; do
        printf '%s\n' "          - target: $target" "            name: $name" "            cpu: $cpu" "            zig: $zig" \
            "            openssl: $package/openssl-9.9.9-habc_0.conda" "            openssl_sha256: $OPENSSL_SUM"
    done <<'ROWS'
x86_64-unknown-linux-gnu x86_64-unknown-linux-gnu x86-64-v3 x86_64-linux linux-64
x86_64-unknown-linux-gnu-v2 x86_64-unknown-linux-gnu x86-64-v2 x86_64-linux linux-64
aarch64-unknown-linux-gnu aarch64-unknown-linux-gnu neoverse-v1 aarch64-linux linux-aarch64
aarch64-unknown-linux-gnu-generic aarch64-unknown-linux-gnu generic aarch64-linux linux-aarch64
ROWS
    printf '%s\n' '          - target: aarch64-apple-darwin' '            name: aarch64-apple-darwin'
} >"$TREE/.github/workflows/release.yml"

# A fabricated tree is a checkout with one commit at a constant time, the time
# stage_tree.sh reads. The blank LFS filters keep `git add` storing the files
# themselves on a host where git-lfs is installed.
COMMIT_TIME=1700000000
init_tree() { # init_tree <dir>: a checkout with no commit
    git -C "$1" init -q >/dev/null 2>&1 &&
        git -C "$1" config user.name "stage tree test" && git -C "$1" config user.email stage-tree-test@example.invalid &&
        git -C "$1" config commit.gpgsign false && git -C "$1" config filter.lfs.clean "" &&
        git -C "$1" config filter.lfs.process "" && git -C "$1" config filter.lfs.required false
}
commit_tree() { # commit_tree <dir>: every file of <dir> in one commit at COMMIT_TIME
    init_tree "$1" && git -C "$1" add -A &&
        GIT_AUTHOR_DATE="$COMMIT_TIME +0000" GIT_COMMITTER_DATE="$COMMIT_TIME +0000" \
            git -C "$1" commit -q -m "fabricated tree" >/dev/null 2>&1
}
commit_tree "$TREE" || bad "git could not commit the fabricated tree"

# curl: logs the URL and writes the fabricated download for it; STUB_CURL_FAIL
# names a host whose download fails, STUB_CURL_BYTES replaces the zig tarball.
cat >"$STUBS/curl" <<'STUB'
#!/usr/bin/env bash
echo "curl $*" >>"$STUB_LOG"
[[ "$1" == -fsSLo ]] || exit 2
[[ -z "${STUB_CURL_FAIL:-}" || "$3" != *"$STUB_CURL_FAIL"* ]] || exit 22
case "$3" in
    https://ziglang.org/*) echo "${STUB_CURL_BYTES:-zig tarball}" >"$2" ;;
    https://conda.anaconda.org/*) echo "openssl package" >"$2" ;;
    *) exit 2 ;;
esac
STUB
cat >"$STUBS/tar" <<'STUB'
#!/usr/bin/env bash
echo "tar $*" >>"$STUB_LOG"
STUB
# unzip: leaves the member the script's glob reads.
cat >"$STUBS/unzip" <<'STUB'
#!/usr/bin/env bash
echo "unzip $*" >>"$STUB_LOG"
while [[ $# -gt 0 ]]; do
    [[ "$1" == -d ]] && : >"$2/pkg-stub.tar.zst"
    shift
done
STUB
# cargo: logs its arguments and, for `install`, leaves the stub cargo-zigbuild
# under --root. A cargo-zigbuild on PATH is a decoy the build must not run.
cat >"$STUBS/cargo" <<'STUB'
#!/usr/bin/env bash
printf '%s\n' "cargo cwd=$PWD" "cargo args=$*" >>"$STUB_LOG"
[[ "$1 ${2:-}" == "install cargo-zigbuild" ]] || exit 2
[[ "${STUB_INSTALL_RC:-0}" -eq 0 ]] || exit "$STUB_INSTALL_RC"
while [[ $# -gt 0 ]]; do
    [[ "$1" == --root ]] && mkdir -p "$2/bin" && cp "$STUB_ZIGBUILD" "$2/bin/cargo-zigbuild"
    shift
done
STUB
cat >"$STUBS/cargo-zigbuild" <<'STUB'
#!/usr/bin/env bash
echo "decoy cargo-zigbuild $*" >>"$STUB_LOG"
exit 2
STUB
# The installed cargo-zigbuild: logs where it ran, the environment the build
# reads and its arguments, then leaves the three binaries where a build for
# --target does; STUB_CARGO_SKIP names one to leave out.
cat >"$TMP/cargo-zigbuild" <<'STUB'
#!/usr/bin/env bash
printf '%s\n' "zigbuild self=$0" "zigbuild cwd=$PWD" "zigbuild rustflags=${RUSTFLAGS:-}" "zigbuild path=$PATH" \
    "zigbuild clock=${SOURCE_DATE_EPOCH:-unset}" \
    "zigbuild openssl=$(env | grep '_OPENSSL_DIR=' | sort | tr '\n' ' ')" "zigbuild args=$*" >>"$STUB_LOG"
[[ "$1 ${2:-}" == "zigbuild --release" ]] || exit 2
[[ "${STUB_CARGO_RC:-0}" -eq 0 ]] || exit "$STUB_CARGO_RC"
triple=""
while [[ $# -gt 0 ]]; do
    [[ "$1" == --target ]] && triple="${2%.*.*}"
    shift
done
out="${CARGO_TARGET_DIR:-target}/$triple/release"
mkdir -p "$out"
for bin in vep vep-cache-builder vep-cache-converter; do
    [[ "$bin" == "${STUB_CARGO_SKIP:-}" ]] || echo "$bin $triple" >"$out/$bin"
done
STUB
cat >"$STUBS/python3" <<'STUB'
#!/usr/bin/env bash
echo "python3 $*" >>"$STUB_LOG"
exit "${STUB_ELF_RC:-0}"
STUB
cat >"$STUBS/uname" <<'STUB'
#!/usr/bin/env bash
[[ "$1" == -m ]] && echo "${STUB_ARCH:-aarch64}"
STUB
# gpg: the line git prints ahead of a signed commit where its configuration
# shows signatures.
cat >"$STUBS/gpg" <<'STUB'
#!/usr/bin/env bash
echo "gpg: fabricated signature" >&2
exit 1
STUB
chmod +x "$STUBS"/* "$TMP/cargo-zigbuild"

run() { # run <dest> <tools> [VAR=value ...]; TARGET in the caller's env overrides the generic arm64 archive
    local dest="$1" tools="$2"
    shift 2
    : >"$LOG"
    env -u CARGO_TARGET_DIR -u RUSTFLAGS PATH="$STUBS:$PATH" STUB_LOG="$LOG" STUB_ZIGBUILD="$TMP/cargo-zigbuild" "$@" \
        bash "$STAGE" "${TARGET:-aarch64-unknown-linux-gnu-generic}" "$dest" "$tools"
}

# The furthest step the log shows: none, zig, install, openssl, build or check.
reached() {
    local step=none
    grep -q '^curl .*ziglang' "$LOG" && step=zig
    grep -q '^cargo args=install ' "$LOG" && step=install
    grep -q '^curl .*conda' "$LOG" && step=openssl
    grep -q '^zigbuild args=zigbuild --release ' "$LOG" && step=build
    grep -q '^python3 ' "$LOG" && step=check
    echo "$step"
}

# The runs below take numbered directories: a name drawn at random can repeat
# one an earlier run left, and stage_tree.sh refuses a directory that exists.
SEQ=0

expect_fail() { # expect_fail <name> <stderr substring> <furthest step reached> [VAR=value ...]
    local name="$1" want="$2" step="$3" dest err rc
    shift 3
    SEQ=$((SEQ + 1))
    dest="$TMP/ctx-$SEQ"
    err=$(run "$dest" "$TMP/tools-$SEQ" "$@" 2>&1 >/dev/null)
    rc=$?
    if [[ "$rc" -ne 0 && "$err" == *"$want"* && ! -e "$dest" && "$(reached)" == "$step" ]]; then
        ok "$name"
    else
        bad "$name (rc=$rc, staged=$([[ -e "$dest" ]] && echo yes || echo no), reached=$(reached), stderr: $err)"
    fi
}

expect_build() { # expect_build <name> <archive> <triple> <cpu>: the build cargo is asked for
    local name="$1" target="$2" triple="$3" cpu="$4" dest
    SEQ=$((SEQ + 1))
    dest="$TMP/ctx-$SEQ"
    if TARGET="$target" run "$dest" "$TMP/tools-$SEQ" STUB_ARCH="${triple%%-*}" >/dev/null 2>&1 &&
        grep -qxF "zigbuild args=zigbuild --release --locked --target $triple.2.99 --bin vep --bin vep-cache-builder --bin vep-cache-converter" "$LOG" &&
        grep -qxF "zigbuild rustflags=-C target-cpu=$cpu" "$LOG" && [[ "$(cat "$dest/vep")" == "vep $triple" ]]; then
        ok "$name"
    else
        bad "$name: $(cat "$LOG" 2>&1)"
    fi
}

bash -n "$STAGE" && ok "stage_tree.sh parses" || bad "stage_tree.sh has a syntax error"

dest="$TMP/ctx-good"
tools="$TMP/tools-good"
if run "$dest" "$tools" >/dev/null 2>&1 &&
    [[ -x "$dest/vep" && -x "$dest/vep-cache-builder" && -x "$dest/vep-cache-converter" &&
        "$(cat "$dest/LICENSE")" == license && "$(cat "$dest/NOTICE")" == notice &&
        "$(cat "$dest/vep-cache-converter")" == "vep-cache-converter aarch64-unknown-linux-gnu" ]]; then
    ok "a build stages three executables, LICENSE and NOTICE"
else
    bad "a build did not stage the context: $(ls -l "$dest" 2>&1)"
fi
grep -qxF "cargo cwd=$TREE" "$LOG" && grep -qxF "zigbuild cwd=$TREE" "$LOG" && ok "cargo runs at the root of the tree" || bad "cargo ran elsewhere: $(cat "$LOG")"
[[ -f "$TREE/target/aarch64-unknown-linux-gnu/release/vep" ]] &&
    ok "the binaries are read from the tree's target directory" ||
    bad "no build under $TREE/target"
grep -qxF "curl -fsSLo $tools/zig.tar.xz https://ziglang.org/download/0.1.2/zig-aarch64-linux-0.1.2.tar.xz" "$LOG" &&
    grep -qxF "tar -xJf $tools/zig.tar.xz -C $tools/zig --strip-components=1" "$LOG" &&
    ok "zig is downloaded at the recipe's version for the row's architecture and unpacked under the tools directory" ||
    bad "the zig download: $(grep -E '^(curl|tar) ' "$LOG")"
grep -qxF "cargo args=install cargo-zigbuild --locked --version 4.5.6 --root $tools/cargo" "$LOG" &&
    ok "cargo-zigbuild is installed at the recipe's version under the tools directory" ||
    bad "the cargo-zigbuild install: $(grep '^cargo args=install' "$LOG")"
grep -qxF "curl -fsSLo $tools/openssl.conda https://conda.anaconda.org/conda-forge/linux-aarch64/openssl-9.9.9-habc_0.conda" "$LOG" &&
    grep -qxF "tar --zstd -xf $tools/pkg-stub.tar.zst -C $tools/openssl lib include" "$LOG" &&
    ok "the row's OpenSSL package is downloaded and unpacked under the tools directory" ||
    bad "the OpenSSL download: $(grep -E '^(curl|unzip|tar) ' "$LOG")"
grep -qxF "zigbuild args=zigbuild --release --locked --target aarch64-unknown-linux-gnu.2.99 --bin vep --bin vep-cache-builder --bin vep-cache-converter" "$LOG" &&
    grep -qxF "zigbuild rustflags=-C target-cpu=generic" "$LOG" &&
    ok "the build links against the recipe's glibc floor with the row's target-cpu" ||
    bad "the build: $(grep '^zigbuild ' "$LOG")"
grep -qxF "zigbuild self=$tools/cargo/bin/cargo-zigbuild" "$LOG" && ! grep -q '^decoy ' "$LOG" &&
    grep -q "^zigbuild path=$tools/zig:" "$LOG" &&
    grep -qxF "zigbuild openssl=AARCH64_UNKNOWN_LINUX_GNU_OPENSSL_DIR=$tools/openssl " "$LOG" &&
    ok "the build runs the cargo-zigbuild installed under the tools directory, with its zig first on PATH and the target's OpenSSL" ||
    bad "the build's environment: $(grep -E '^(zigbuild|decoy) ' "$LOG")"
out="$TREE/target/aarch64-unknown-linux-gnu/release"
grep -qxF "python3 $TREE/scripts/release/check_glibc_floor.py --floor 2.99 $out/vep $out/vep-cache-builder $out/vep-cache-converter" "$LOG" &&
    ok "the release's ELF check holds the three binaries to the recipe's glibc floor" ||
    bad "the ELF check: $(grep '^python3 ' "$LOG")"

dest="$TMP/ctx-clock"
if run "$dest" "$TMP/tools-clock" SOURCE_DATE_EPOCH=1 >/dev/null 2>&1 && grep -qxF "zigbuild clock=$COMMIT_TIME" "$LOG"; then
    ok "the build's clock is the time of the tree's commit, whatever SOURCE_DATE_EPOCH the caller holds"
else
    bad "the build's clock, $COMMIT_TIME expected: $(grep '^zigbuild clock=' "$LOG")"
fi

sign_tree() { # sign_tree <dir>: the commit of <dir> again, signed, in a checkout that shows signatures
    local signed
    signed=$({
        git -C "$1" cat-file commit HEAD | sed -n '1,/^committer /p'
        printf '%s\n' 'gpgsig -----BEGIN PGP SIGNATURE-----' ' fabricated' ' -----END PGP SIGNATURE-----' '' 'fabricated tree'
    } | git -C "$1" hash-object -t commit -w --stdin) && git -C "$1" update-ref HEAD "$signed" &&
        git -C "$1" config gpg.program "$STUBS/gpg" && git -C "$1" config log.showSignature true
}
sign_tree "$TREE" || bad "git could not sign the fabricated tree's commit"
dest="$TMP/ctx-signed"
if run "$dest" "$TMP/tools-signed" >/dev/null 2>&1 && grep -qxF "zigbuild clock=$COMMIT_TIME" "$LOG"; then
    ok "the build's clock is the time of a signed commit where git's configuration shows signatures"
else
    bad "the clock of a signed commit, $COMMIT_TIME expected: $(grep '^zigbuild clock=' "$LOG")"
fi
git -C "$TREE" config --unset log.showSignature

dest="$TMP/ctx-target-dir"
if run "$dest" "$TMP/tools-target-dir" CARGO_TARGET_DIR="$TMP/elsewhere" >/dev/null 2>&1 &&
    [[ -f "$TMP/elsewhere/aarch64-unknown-linux-gnu/release/vep" && -x "$dest/vep" ]]; then
    ok "CARGO_TARGET_DIR moves the build the binaries are read from"
else
    bad "CARGO_TARGET_DIR was not honoured: $(ls -l "$dest" 2>&1)"
fi

expect_build "the generic x86_64 archive builds for x86-64-v2" x86_64-unknown-linux-gnu-v2 x86_64-unknown-linux-gnu x86-64-v2
expect_build "the measured x86_64 archive builds for x86-64-v3" x86_64-unknown-linux-gnu x86_64-unknown-linux-gnu x86-64-v3
expect_build "the generic arm64 archive builds for the generic CPU" aarch64-unknown-linux-gnu-generic aarch64-unknown-linux-gnu generic
expect_build "the measured arm64 archive builds for neoverse-v1" aarch64-unknown-linux-gnu aarch64-unknown-linux-gnu neoverse-v1

# The tree's own release.yml in the fabricated tree, behind a sha256sum that
# records the checksum it is handed: `name target cpu zig openssl
# openssl_sha256` of every Linux row, read here line by line.
REAL="$TMP/real"
mkdir -p "$REAL/scripts/docker" "$REAL/scripts/release" "$REAL/.github/workflows" "$TMP/sumbin"
cp "$HERE/stage_tree.sh" "$HERE/release_recipe.sh" "$REAL/scripts/docker/" 2>/dev/null
cp "$TREE/LICENSE" "$TREE/NOTICE" "$REAL/"
cp "$RELEASE_YML" "$REAL/.github/workflows/release.yml"
commit_tree "$REAL" || bad "git could not commit the tree that holds the real release.yml"
cat >"$TMP/sumbin/sha256sum" <<'STUB'
#!/usr/bin/env bash
echo "sha256sum $(cat)" >>"$STUB_LOG"
STUB
chmod +x "$TMP/sumbin/sha256sum"
rows=$(awk '
    /^          - target: / { target = $3; name = ""; cpu = ""; zig = ""; openssl = "" }
    /^            name: / { name = $2 }
    /^            cpu: / { cpu = $2 }
    /^            zig: / { zig = $2 }
    /^            openssl: / { openssl = $2 }
    /^            openssl_sha256: / && name != "" { print name, target, cpu, zig, openssl, $2 }' "$RELEASE_YML")
stated() { sed -nE "s/^  $1: \"?([^\" ]+)\"?\$/\1/p" "$RELEASE_YML"; }
floor=$(stated GLIBC_FLOOR)
zig_version=$(stated ZIG_VERSION)
zigbuild=$(stated CARGO_ZIGBUILD_VERSION)
real() { # real <archive> <triple> <tools>: stage_tree.sh on the tree's release.yml
    : >"$LOG"
    SEQ=$((SEQ + 1))
    env -u CARGO_TARGET_DIR -u RUSTFLAGS PATH="$TMP/sumbin:$STUBS:$PATH" STUB_LOG="$LOG" STUB_ARCH="${2%%-*}" STUB_ZIGBUILD="$TMP/cargo-zigbuild" \
        bash "$REAL/scripts/docker/stage_tree.sh" "$1" "$TMP/ctx-$SEQ" "$3" >/dev/null 2>&1
}
count=0
wrong=""
while read -r name triple cpu zig openssl openssl_sha256; do
    [[ -n "$name" ]] || continue
    count=$((count + 1))
    tools="$TMP/tools-real-$count"
    case "$zig" in
        x86_64-linux) sum=$(stated ZIG_SHA256_X86_64_LINUX) ;;
        aarch64-linux) sum=$(stated ZIG_SHA256_AARCH64_LINUX) ;;
        *) sum=unread ;;
    esac
    real "$name" "$triple" "$tools" &&
        grep -qF -- "--target $triple.$floor --bin vep " "$LOG" && grep -qxF "zigbuild rustflags=-C target-cpu=$cpu" "$LOG" &&
        grep -qF "/download/$zig_version/zig-$zig-$zig_version.tar.xz" "$LOG" && grep -qxF "sha256sum $sum  $tools/zig.tar.xz" "$LOG" &&
        grep -qF -- "--version $zigbuild --root " "$LOG" &&
        grep -qF "/conda-forge/$openssl" "$LOG" && grep -qxF "sha256sum $openssl_sha256  $tools/openssl.conda" "$LOG" || wrong+=" $name"
done <<<"$rows"
[[ "$count" -eq 4 && -n "$floor" && -z "$wrong" ]] &&
    ok "every Linux archive builds by its release.yml row: target, cpu, glibc floor, zig, cargo-zigbuild and OpenSSL, each checksum the recorded one" ||
    bad "$count Linux rows read from release.yml, 4 expected; built otherwise:$wrong"

archives=$(sed -nE 's/^  ARCHIVES: (.+)$/\1/p' "$DOCKER_YML")
count=0
wrong=""
for name in $archives; do
    count=$((count + 1))
    row=$(grep "^$name " <<<"$rows") || { wrong+=" $name"; continue; }
    read -r _ triple _ <<<"$row"
    real "$name" "$triple" "$TMP/tools-docker-$count" || wrong+=" $name"
done
[[ "$count" -eq 4 && -z "$wrong" ]] && ok "every archive docker.yml packages is one stage_tree.sh builds" ||
    bad "$count archives read from docker.yml, 4 expected; not built, or not in release.yml:$wrong"

# The commands of release.yml's install and build steps, which stage_tree.sh
# repeats with the same variable names.
missing=""
while IFS= read -r fragment; do
    grep -qF -- "$fragment" "$RELEASE_YML" && grep -qF -- "$fragment" "$HERE/stage_tree.sh" || missing+=" [$fragment]"
done <<'FRAGMENTS'
"https://ziglang.org/download/$ZIG_VERSION/zig-$ZIG_ARCH-$ZIG_VERSION.tar.xz"
| sha256sum --check --strict
--strip-components=1
cargo install cargo-zigbuild --locked --version "$CARGO_ZIGBUILD_VERSION"
"https://conda.anaconda.org/conda-forge/$OPENSSL_PACKAGE"
'pkg-*.tar.zst' -d
/pkg-*.tar.zst -C
/openssl" lib include
| tr 'a-z-' 'A-Z_')_OPENSSL_DIR=
-C target-cpu=
zigbuild --release --locked --target
--bin vep --bin vep-cache-builder --bin vep-cache-converter
SOURCE_DATE_EPOCH=
log -1 --format=%ct
check_glibc_floor.py
--floor "$GLIBC_FLOOR"
FRAGMENTS
[[ -z "$missing" ]] && ok "stage_tree.sh downloads, installs, builds and checks with the commands of release.yml" ||
    bad "a command release.yml and stage_tree.sh do not share:$missing"

# release.yml's own clock: one unconditional step of the job `build`, ahead
# of both build steps, writes the commit's time for every later step.
verdict=$(awk '
    /^  [a-z][a-z0-9_-]*:$/ { job = $1 }
    job != "build:" { next }
    /^      - / { step++ }
    /^      (- |  )if: / { conditional[step] = 1 }
    index($0, "echo \"SOURCE_DATE_EPOCH=$(git log -1 --format=%ct)\" >> \"$GITHUB_ENV\"") { clocks++; clock = step }
    /^ +cargo (zig)?build --release --locked --target / { builds++; if (!first) first = step }
    END { print clocks + 0, builds + 0, (clock && first && clock < first && !conditional[clock]) ? "ahead" : "not ahead" }' "$RELEASE_YML")
if [[ "$verdict" == "1 2 ahead" ]]; then
    ok "release.yml writes SOURCE_DATE_EPOCH, the commit's time, to GITHUB_ENV in an unconditional step ahead of both build steps"
else
    bad "release.yml's build clock (clock steps, build steps, order): $verdict"
fi

no_commit() { # no_commit <name> <tree> [init]: stage_tree.sh from a tree with no commit of its own
    local name="$1" tree="$2" err rc
    SEQ=$((SEQ + 1))
    mkdir -p "$tree/scripts/docker" "$tree/.github/workflows"
    cp "$HERE/stage_tree.sh" "$HERE/release_recipe.sh" "$tree/scripts/docker/" 2>/dev/null
    cp "$TREE/.github/workflows/release.yml" "$tree/.github/workflows/"
    [[ "${3:-}" != init ]] || init_tree "$tree"
    : >"$LOG"
    err=$(env -u CARGO_TARGET_DIR -u RUSTFLAGS PATH="$STUBS:$PATH" STUB_LOG="$LOG" STUB_ZIGBUILD="$TMP/cargo-zigbuild" \
        bash "$tree/scripts/docker/stage_tree.sh" aarch64-unknown-linux-gnu-generic "$TMP/ctx-$SEQ" "$TMP/tools-$SEQ" 2>&1 >/dev/null)
    rc=$?
    if [[ "$rc" -ne 0 && "$err" == *"has no commit of its own to take the build clock from"* &&
        ! -e "$TMP/ctx-$SEQ" && ! -e "$TMP/tools-$SEQ" && ! -s "$LOG" ]]; then
        ok "$name"
    else
        bad "$name (rc=$rc, staged=$([[ -e "$TMP/ctx-$SEQ" ]] && echo yes || echo no), calls: $(cat "$LOG"), stderr: $err)"
    fi
}
no_commit "a tree that is no git checkout stops before any download, with nothing staged" "$TMP/unpacked"
no_commit "a checkout with no commit stops before any download, with nothing staged" "$TMP/uncommitted" init
no_commit "a tree unpacked inside another checkout does not take that checkout's commit time" "$TREE/unpacked"

TARGET=aarch64-apple-darwin expect_fail "a macOS target fails" "release.yml states no build recipe for 'aarch64-apple-darwin'" none
TARGET=x86_64-unknown-linux-gnu-v4 expect_fail "an archive the release does not build fails" "release.yml states no build recipe for 'x86_64-unknown-linux-gnu-v4'" none
TARGET=x86_64-unknown-linux-gnu/../x expect_fail "a target with a path in it fails" "release.yml states no build recipe for" none
mkdir -p "$TMP/exists"
err=$(run "$TMP/exists" "$TMP/tools-exists" 2>&1 >/dev/null) && bad "an existing destination was overwritten" ||
    { [[ "$err" == *"$TMP/exists already exists"* && ! -s "$LOG" && ! -e "$TMP/exists/vep" ]] && ok "an existing destination is never overwritten" ||
        bad "existing destination: $err"; }
mkdir -p "$TMP/tools-kept"
err=$(run "$TMP/ctx-kept" "$TMP/tools-kept" 2>&1 >/dev/null) && bad "an existing tools directory was reused" ||
    { [[ "$err" == *"$TMP/tools-kept already exists"* && ! -s "$LOG" && ! -e "$TMP/ctx-kept" ]] && ok "an existing tools directory is never reused" ||
        bad "existing tools directory: $err"; }
expect_fail "a host of another architecture fails" "needs a host of architecture aarch64, not x86_64" none STUB_ARCH=x86_64
expect_fail "a failed zig download fails" "could not download zig 0.1.2 for aarch64-linux" zig STUB_CURL_FAIL=ziglang.org
expect_fail "a zig tarball that is not the recorded one stops the run before anything is installed" \
    "the zig 0.1.2 tarball does not match the checksum release.yml records" zig STUB_CURL_BYTES=other
expect_fail "a failed cargo-zigbuild install fails" "cargo could not install cargo-zigbuild 4.5.6" install STUB_INSTALL_RC=101
expect_fail "a failed OpenSSL download fails" "could not download linux-aarch64/openssl-9.9.9-habc_0.conda" openssl STUB_CURL_FAIL=conda.anaconda.org
sed -E 's/^(            openssl_sha256: ).*/\1'"$ZIG_SUM"'/' "$TREE/.github/workflows/release.yml" >"$TMP/other-sum.yml"
expect_fail "an OpenSSL package that is not the recorded one stops the run before the build" \
    "openssl-9.9.9-habc_0.conda does not match the checksum release.yml records" openssl WORKFLOW="$TMP/other-sum.yml"
expect_fail "a failed build fails" "cargo could not build the three binaries" build STUB_CARGO_RC=101
expect_fail "a build without vep-cache-converter leaves nothing staged" "the build left no vep-cache-converter" build \
    STUB_CARGO_SKIP=vep-cache-converter CARGO_TARGET_DIR="$TMP/partial"
expect_fail "a binary that needs a glibc above the floor leaves nothing staged" \
    "a binary built for aarch64-unknown-linux-gnu-generic needs a glibc newer than the floor 2.99" check STUB_ELF_RC=1

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
