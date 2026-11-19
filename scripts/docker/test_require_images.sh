#!/usr/bin/env bash
# Tests for require_images.sh: every pair of an event and the two image jobs'
# results, the errors that name the job, and the copy of the script's body
# that the step `Require every image` of .github/workflows/docker.yml runs.
#
# Pure bash. No network, no GitHub.
# Run: bash scripts/docker/test_require_images.sh
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
REQUIRE="$HERE/require_images.sh"
ROOT=$(cd "$HERE/../.." && pwd)
WORKFLOW="$ROOT/.github/workflows/docker.yml"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

run() { # run [VAR=value ...]: require_images.sh in an environment of those variables alone
    env -i PATH="$PATH" "$@" bash "$REQUIRE"
}

bash -n "$REQUIRE" 2>/dev/null && ok "require_images.sh parses" || bad "require_images.sh is missing or has a syntax error"

# Every pair: an event passes on one pair of results alone, its owed job's
# `success` with the other job `skipped`.
RESULTS="success failure cancelled skipped"
table() { # table <name> <event> <owed job> <tree result that passes> <image result that passes>
    local name="$1" event="$2" owed="$3" tree image out rc pairs=0 passed="" wrong=""
    for tree in $RESULTS; do
        for image in $RESULTS; do
            pairs=$((pairs + 1))
            out=$(run EVENT="$event" TREE_RESULT="$tree" IMAGE_RESULT="$image" 2>&1)
            rc=$?
            if [[ "$tree" == "$4" && "$image" == "$5" ]]; then
                [[ "$rc" -eq 0 && "$out" == "every image row of the job '$owed' succeeded" ]] && passed="$tree/$image" || wrong+=" $tree/$image:rc=$rc"
            else
                [[ "$rc" -eq 1 && "$out" == "ERROR: [require_images] the job '"* && "$(wc -l <<<"$out")" -eq 1 ]] || wrong+=" $tree/$image:rc=$rc"
            fi
        done
    done
    [[ "$pairs" -eq 16 && "$passed" == "$4/$5" && -z "$wrong" ]] && ok "$name" ||
        bad "$name: $pairs pairs of 16, passed [$passed], decided otherwise:$wrong"
}
table "a pull request passes on the job tree's success with the job image skipped, and on no other pair" pull_request tree success skipped
table "a release passes on the job image's success with the job tree skipped, and on no other pair" release image skipped success
table "a dispatch passes on the job image's success with the job tree skipped, and on no other pair" workflow_dispatch image skipped success

expect_fail() { # expect_fail <name> <the whole error line> [VAR=value ...]
    local name="$1" want="$2" out rc
    shift 2
    out=$(run "$@" 2>&1)
    rc=$?
    if [[ "$rc" -eq 1 && "$out" == "$want" ]]; then ok "$name"; else bad "$name (rc=$rc, output: $out)"; fi
}
expect_fail "a pull request whose rows failed names the job tree" \
    "ERROR: [require_images] the job 'tree', which a pull_request owes, ended 'failure'" \
    EVENT=pull_request TREE_RESULT=failure IMAGE_RESULT=skipped
expect_fail "a pull request whose rows were skipped is an error, not a pass" \
    "ERROR: [require_images] the job 'tree', which a pull_request owes, ended 'skipped'" \
    EVENT=pull_request TREE_RESULT=skipped IMAGE_RESULT=skipped
expect_fail "a release whose rows were cancelled names the job image" \
    "ERROR: [require_images] the job 'image', which a release owes, ended 'cancelled'" \
    EVENT=release TREE_RESULT=skipped IMAGE_RESULT=cancelled
expect_fail "a dispatch whose rows were skipped names the job image" \
    "ERROR: [require_images] the job 'image', which a workflow_dispatch owes, ended 'skipped'" \
    EVENT=workflow_dispatch TREE_RESULT=skipped IMAGE_RESULT=skipped
expect_fail "a pull request on which the pushing job ran is an error that names it" \
    "ERROR: [require_images] the job 'image', which a pull_request leaves out, ended 'success'" \
    EVENT=pull_request TREE_RESULT=success IMAGE_RESULT=success
expect_fail "a release on which the tree was built is an error that names the job" \
    "ERROR: [require_images] the job 'tree', which a release leaves out, ended 'success'" \
    EVENT=release TREE_RESULT=success IMAGE_RESULT=success
expect_fail "a result the run did not report is an error" \
    "ERROR: [require_images] the job 'tree', which a pull_request owes, ended ''" EVENT=pull_request
expect_fail "an event the workflow builds no image for fails" \
    "ERROR: [require_images] event 'push' is not one the workflow builds an image for" \
    EVENT=push TREE_RESULT=success IMAGE_RESULT=success
expect_fail "a run that names no event fails" \
    "ERROR: [require_images] event '' is not one the workflow builds an image for" TREE_RESULT=success IMAGE_RESULT=skipped

# The step runs the script's body, from `set -euo pipefail` to its last line,
# as its `run:` block: its job checks out no tree before the requirement.
step=$(awk '
    /^      - / || /^  [a-z_-]+:$/ { inside = 0 }
    /^      - name: Require every image$/ { inside = 1 }
    inside' "$WORKFLOW")
body=$(awk '
    /^        run: \|$/ { block = 1; next }
    block && /^          / { print substr($0, 11); next }
    block && /^$/ { print ""; next }
    block { exit }' <<<"$step")
script=$(sed -n '/^set -euo pipefail$/,$p' "$REQUIRE" 2>/dev/null)
[[ -n "$script" && "$body" == "$script" ]] && ok "the step Require every image runs the script's body, line for line" ||
    bad "the step's run block differs from require_images.sh: $(diff <(echo "$body") <(echo "$script") | head -n 6)"
grep -qxF '          EVENT: ${{ github.event_name }}' <<<"$step" &&
    grep -qxF '          TREE_RESULT: ${{ needs.tree.result }}' <<<"$step" &&
    grep -qxF '          IMAGE_RESULT: ${{ needs.image.result }}' <<<"$step" &&
    ok "the step hands the script the event and the result of each image job" ||
    bad "the step's env does not name the event and both jobs' results: $(grep -E '^          [A-Z_]+: ' <<<"$step")"
! grep -qE '^        (uses|if): |scripts/' <<<"$step" && [[ -n "$body" ]] &&
    grep -q '^    needs: \[tree, image\]$' "$WORKFLOW" &&
    ok "the step runs on every event the job runs for, after both image jobs, and reads no file of a tree" ||
    bad "the step is conditional, uses an action, runs a file of the tree, or its job does not need both image jobs"

echo
echo "passed: $PASS   failed: $FAIL"
[[ "$FAIL" -eq 0 ]] || exit 1
