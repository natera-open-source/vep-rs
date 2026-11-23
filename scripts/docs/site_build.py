#!/usr/bin/env python3
"""Build the documentation site and refuse a warning of the builder.

Runs `mdbook build` with the arguments given, in the directory it is started
from, and passes the builder's lines through as they come. mdBook builds the
site under a warning, so a fault it only warns of reaches the pages: an HTML
tag a page leaves unclosed is written as an element that takes in the text
after it. A placeholder such as `<dir>` is such a tag outside a code span, and
inside one that wraps so that the placeholder opens a line, where an element
name of CommonMark's HTML-block list ends the paragraph and the span with it.

A warning is a line whose first word is the level `WARN`, as mdBook 0.5 writes
it; the lines under it carry on its message. `ADMITTED` lists the warnings
every build of the tree prints, by the words their message opens with, each
with its reason. The builder runs without `MDBOOK_LOG`, which sets the level
it logs from and puts each line's source before its message.

Usage::

    python3 scripts/docs/site_build.py [MDBOOK BUILD ARGUMENT ...]

Exits with the builder's status when the build fails, 1 when it succeeds under
a warning `ADMITTED` does not list (one `ERROR:` line a warning, after the
build), and 2 when `mdbook` is not on PATH. Stdlib only.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys

BUILDER = "mdbook"
# Set to `error` it silences the builder's warnings, and set at all it changes the lines `WARNING` reads.
LOG_VARIABLE = "MDBOOK_LOG"
# The mdBook series whose log lines `WARNING` reads; the Pages workflow pins a version of it.
MDBOOK_SERIES = "0.5"
WARNING = re.compile(r"^\s*WARN\s+(.*?)\s*$")
# the words a warning's message opens with -> why every build of the tree prints it
ADMITTED = {
    "search index is very large": "the search index of the site's pages is larger than the size mdBook warns at",
}


def refused(lines):
    """The messages of the warning lines among `lines` that `ADMITTED` does not list, in order."""
    messages = (match.group(1) for match in map(WARNING.match, lines) if match)
    return [message for message in messages if not message.startswith(tuple(ADMITTED))]


def main(arguments=None):
    arguments = sys.argv[1:] if arguments is None else arguments
    if shutil.which(BUILDER) is None:
        print(f"ERROR: [site_build] {BUILDER} is not on PATH", file=sys.stderr)
        return 2
    build = subprocess.Popen(
        [BUILDER, "build", *arguments],
        env={name: value for name, value in os.environ.items() if name != LOG_VARIABLE},
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        encoding="utf-8",
        errors="replace",
    )
    lines = []
    for line in build.stdout:
        try:
            sys.stdout.write(line)
        except UnicodeEncodeError:
            encoding = sys.stdout.encoding
            sys.stdout.write(line.encode(encoding, errors="backslashreplace").decode(encoding))
        sys.stdout.flush()
        lines.append(line)
    status = build.wait()
    if status != 0:
        print(f"ERROR: [site_build] {BUILDER} build ended with status {status}", file=sys.stderr)
        return status if status > 0 else 1
    warnings = refused(lines)
    for message in warnings:
        print(f"ERROR: [site_build] the builder warned: {message}", file=sys.stderr)
    return 1 if warnings else 0


if __name__ == "__main__":
    sys.exit(main())
