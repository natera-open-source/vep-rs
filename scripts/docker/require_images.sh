#!/usr/bin/env bash
# Decides the check `Build, smoke-test and publish` of
# .github/workflows/docker.yml from the event and the results of its two image
# jobs. Each event owes one of them: a pull request the job `tree`, which
# builds the tree's images and holds no write permission; a release or a
# dispatch the job `image`, which packages a release and can push. The owed
# job must have ended `success` and the other must have been skipped, as its
# `if:` leaves it out of the event; any other pair exits 1 and names the job.
#
# The check's job holds write permissions and checks out no tree on a pull
# request, so its step `Require every image` carries the lines from
# `set -euo pipefail` on in place of a call to this file;
# test_require_images.sh holds the two equal.
#
# Usage: EVENT=<event> TREE_RESULT=<result> IMAGE_RESULT=<result>
#            scripts/docker/require_images.sh
set -euo pipefail

fail() {
    echo "ERROR: [require_images] $*" >&2
    exit 1
}

case "${EVENT:-}" in
    pull_request) owed=tree owed_result="${TREE_RESULT:-}" other=image other_result="${IMAGE_RESULT:-}" ;;
    release | workflow_dispatch) owed=image owed_result="${IMAGE_RESULT:-}" other=tree other_result="${TREE_RESULT:-}" ;;
    *) fail "event '${EVENT:-}' is not one the workflow builds an image for" ;;
esac
[[ "$owed_result" == success ]] || fail "the job '$owed', which a $EVENT owes, ended '$owed_result'"
[[ "$other_result" == skipped ]] || fail "the job '$other', which a $EVENT leaves out, ended '$other_result'"
echo "every image row of the job '$owed' succeeded"
