#!/usr/bin/env bash
# Tests for the jobs of .github/workflows/docker.yml, read as text: the job a
# pull request runs can read the repository and nothing else and holds no step
# that logs in, pushes, attests, records or uploads; no job that declares a
# write permission builds the tree, installs a toolchain or runs cargo; the
# job that pushes runs for a release or a dispatch alone; and the two image
# jobs state the same rows and build the image with the same arguments.
#
# Pure bash. No network, no GitHub.
# Run: bash scripts/docker/test_workflow_jobs.sh
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
ROOT=$(cd "$HERE/../.." && pwd)
WORKFLOW="$ROOT/.github/workflows/docker.yml"

jobs() { # the ids of the workflow's jobs
    awk '/^jobs:$/ { inside = 1; next } inside && /^  [a-z_-]+:$/ { sub(/^  /, ""); sub(/:$/, ""); print }' "$WORKFLOW"
}
job() { # job <id>: the lines of the job, comments left out
    awk -v id="$1" '
        /^jobs:$/ { jobs = 1; next }
        jobs && /^  [a-z_-]+:$/ { inside = ($0 == "  " id ":") }
        inside && !/^ *#/' "$WORKFLOW"
}
steps() { # steps <id>: one line per step of the job, `<its if:, or always> | <its name, action or first line>`
    job "$1" | awk '
        function flush() { if (label != "") print cond " | " label; label = ""; cond = "always" }
        BEGIN { cond = "always" }
        /^      - / { flush(); label = $0; sub(/^      - (name: |uses: |run: )?/, "", label) }
        /^        if: / { cond = $0; sub(/^        if: /, "", cond) }
        END { flush() }'
}
permissions() { # permissions <id>: the job's own permissions, one `scope: access` per line
    job "$1" | awk '/^    permissions:$/ { inside = 1; next } inside && /^      [a-z-]+: / { sub(/^      /, ""); print; next } { inside = 0 }'
}
rows() { # rows <id>: `archive platform os` of each matrix row of the job
    job "$1" | awk '
        /^          - archive: / { archive = $3 }
        /^            platform: / { platform = $2 }
        /^            os: / { print archive, platform, $2 }'
}

TREE=$(job tree)
IMAGE=$(job image)
PUBLISH=$(job publish)

[[ "$(grep -c '^    if: ' <<<"$TREE")" -eq 1 ]] && grep -qxF "    if: github.event_name == 'pull_request'" <<<"$TREE" &&
    ok "the job tree runs for a pull request alone" ||
    bad "the job tree has the condition [$(grep '^    if: ' <<<"$TREE")]"

[[ "$(permissions tree)" == "contents: read" ]] &&
    [[ "$(awk '/^permissions:$/ { getline; print }' "$WORKFLOW")" == "  contents: read" ]] &&
    ok "the job tree can read the repository and nothing else" ||
    bad "the job tree declares the permissions [$(permissions tree | paste -sd, -)]"

found=""
for needle in 'docker login' 'push=true' 'push-by-digest' '--push' 'actions/attest' 'upload-artifact' 'download-artifact' \
    'imagetools' 'docker pull' 'stage_release.sh' 'gh release' 'gh attestation'; do
    ! grep -qF -- "$needle" <<<"$TREE" || found+=" [$needle]"
done
[[ -n "$TREE" && -z "$found" ]] && ok "the job tree holds no step that logs in, pushes, attests, records, uploads or stages a release" ||
    bad "the job tree is missing or holds:$found"

# Every job inherits the workflow's own `env:` block, so a token named there
# is in the job tree's environment as one named in the job is.
SHARED_ENV=$(awk '/^env:$/ { inside = 1; next } /^[a-z]/ { inside = 0 } inside && !/^ *#/' "$WORKFLOW")
found=""
for needle in 'github.token' 'secrets.' 'GH_TOKEN' 'GITHUB_TOKEN' 'id-token' 'packages:'; do
    ! grep -qF -- "$needle" <<<"$TREE"$'\n'"$SHARED_ENV" || found+=" [$needle]"
done
checkouts=$(grep -c 'uses: actions/checkout@' <<<"$TREE")
[[ "$checkouts" -ge 1 && "$(grep -c '^          persist-credentials: false$' <<<"$TREE")" -eq "$checkouts" && -z "$found" ]] &&
    ! grep -qE '^    env:$' <<<"$TREE" && grep -q '^  ARCHIVES: ' <<<"$SHARED_ENV" &&
    ok "the job tree sets no token in its environment and its checkout keeps none" ||
    bad "the job tree or the workflow's env names a token, the job sets an environment, its checkout keeps credentials or the workflow's env was not read:$found"

conditional=$(steps tree | grep -vc '^always | ')
[[ "$(steps tree | wc -l)" -ge 8 && "$conditional" -eq 0 ]] && ok "no step of the job tree is left to a condition" ||
    bad "$conditional of the job tree's $(steps tree | wc -l) steps carry an if:"

order=$(grep -n -o -E 'scripts/docker/release_recipe\.sh|dtolnay/rust-toolchain@|scripts/docker/stage_tree\.sh|docker buildx build|/smoke/smoke_test\.sh' <<<"$TREE" | cut -d: -f2 | paste -sd' ' -)
[[ "$order" == "scripts/docker/release_recipe.sh dtolnay/rust-toolchain@ scripts/docker/stage_tree.sh docker buildx build /smoke/smoke_test.sh" ]] &&
    ok "the job tree reads the release recipe, installs its toolchain, builds and stages the tree, builds the image and runs the smoke test" ||
    bad "the job tree runs, in order: [$order]"

grep -qF '| tee "$RUNNER_TEMP/smoke.log"' <<<"$TREE" && grep -qF "if grep '^  skip ' \"\$RUNNER_TEMP/smoke.log\"; then" <<<"$TREE" &&
    ok "a check the smoke test skips fails a pull request's row" ||
    bad "the job tree's smoke test step does not fail on a skip line"

writers=0
found=""
for id in $(jobs); do
    permissions "$id" | grep ': write$' >/dev/null || continue
    writers=$((writers + 1))
    for needle in 'stage_tree.sh' 'release_recipe.sh' 'dtolnay/rust-toolchain' 'cargo' 'rustup' 'zigbuild'; do
        ! job "$id" | grep -F -- "$needle" >/dev/null || found+=" $id:[$needle]"
    done
done
[[ "$writers" -eq 2 && -z "$found" ]] && ok "no job that declares a write permission builds the tree, installs a toolchain or runs cargo" ||
    bad "$writers jobs declare a write permission, 2 expected; one of them runs:$found"

grep -qxF "    if: github.event_name == 'workflow_dispatch' || (github.event_name == 'release' && !github.event.release.prerelease)" <<<"$IMAGE" &&
    [[ "$(grep -c '^    if: ' <<<"$IMAGE")" -eq 1 ]] &&
    ok "the job that pushes runs for a release or a dispatch alone" ||
    bad "the job image has the condition [$(grep '^    if: ' <<<"$IMAGE")]"

# A refused dispatch ends at the resolve step, so every step that reads a
# release or a registry comes after it.
position() { steps image | awk -F' [|] ' -v name="$1" '$2 == name { print NR; exit }'; }
resolve=$(position 'Resolve what to package and whether to push')
early=""
for step in 'Check that the release carries the archives' 'Check out the golden corpora and examples of the packaged tree' \
    'Stage the verified release binaries' 'Create the builder' 'Log in to GHCR' 'Build'; do
    at=$(position "$step")
    [[ -n "$resolve" && -n "$at" && "$resolve" -lt "$at" ]] || early+=" [$step]"
done
! grep -qF 'outputs.source' "$WORKFLOW" && [[ -z "$early" ]] &&
    [[ "$(steps image | grep -F 'Stage the verified release binaries')" == "always | Stage the verified release binaries" ]] &&
    [[ "$(steps image | grep -F 'Check that the release carries the archives')" == "always | Check that the release carries the archives" ]] &&
    ok "the job that pushes resolves first, then checks and stages a release's archives on every run, whatever the source answered" ||
    bad "a step of the workflow is left to the source ($(grep -n 'outputs.source' "$WORKFLOW" | paste -sd' ' -)) or precedes the resolve step:$early"

others=$(steps publish | grep -vF ' | Require every image' | grep -vc "^needs.image.outputs.push == 'true' | ")
grep -qxF '    needs: [tree, image]' <<<"$PUBLISH" && grep -qxF '    name: Build, smoke-test and publish' <<<"$PUBLISH" &&
    [[ "$(steps publish | head -n 1)" == "always | Require every image" && "$others" -eq 0 ]] &&
    ok "the check waits for both image jobs, and each step after its requirement runs only when the job image pushed" ||
    bad "the check's job does not need both image jobs, or $others of its steps run without a push by the job image: $(steps publish | paste -sd';' -)"

archives=$(sed -nE 's/^  ARCHIVES: (.+)$/\1/p' "$WORKFLOW")
[[ "$(rows tree | wc -l)" -eq 4 && "$(rows tree)" == "$(rows image)" && "$(rows tree | cut -d' ' -f1 | paste -sd' ' -)" == "$archives" ]] &&
    ok "both image jobs state the rows of ARCHIVES, each with one platform and runner" ||
    bad "the rows differ: tree [$(rows tree | paste -sd';' -)] image [$(rows image | paste -sd';' -)] ARCHIVES [$archives]"

arguments() { grep -A1 -F 'args=(--platform "$PLATFORM"' <<<"$1"; }
[[ "$(arguments "$TREE" | wc -l)" -eq 2 && "$(arguments "$TREE")" == "$(arguments "$IMAGE")" ]] &&
    grep -qxF '          docker buildx build "${args[@]}" --load --tag vep-rs:smoke' <<<"$TREE" &&
    grep -qxF '            docker buildx build "${args[@]}" --load --tag vep-rs:smoke' <<<"$IMAGE" &&
    ok "both image jobs build the image with the same arguments, the tree's as a release row that does not push" ||
    bad "the build arguments differ: tree [$(arguments "$TREE")] image [$(arguments "$IMAGE")]"

lint='        run: docker run --rm -i "$HADOLINT_IMAGE" < docker/Dockerfile'
builder='        run: docker buildx create --name vep-rs --driver docker-container --driver-opt "image=$BUILDKIT_IMAGE" --use'
grep -qxF "$lint" <<<"$TREE" && grep -qxF "$lint" <<<"$IMAGE" && grep -qxF "$builder" <<<"$TREE" && grep -qxF "$builder" <<<"$IMAGE" &&
    [[ "$(grep -cE '^  (HADOLINT|BUILDKIT)_IMAGE: [a-z/]+:v[0-9.]+@sha256:[0-9a-f]{64}$' "$WORKFLOW")" -eq 2 ]] &&
    ! grep -vE '^  (HADOLINT|BUILDKIT)_IMAGE: ' "$WORKFLOW" | grep -E 'hadolint/hadolint:|moby/buildkit:' >/dev/null &&
    ok "both image jobs lint with one pinned hadolint and build on one pinned BuildKit, each stated once" ||
    bad "the lint or builder step differs between the image jobs, or a pinned image is stated twice"

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
