#!/usr/bin/env python3
"""mdBook preprocessor that makes the repository's markdown read correctly as a site.

The book's source is `docs/`, but those pages link outside it (`../README.md`,
`../scripts/README.md`, the provenance JSON), and some chapters show repository
files whose links are relative to the repository root. mdBook copies only
`docs/` into the site, so those links would 404. This preprocessor:

- replaces the content of the `changelog.md` chapter with `CHANGELOG.md` and of
  the `contributing.md` chapter with `CONTRIBUTING.md`, so each has one source;
  it also replaces the `index.md` chapter with `README.md` when that chapter is
  the three-line `# vep-rs documentation` stub, or when `$VEP_RS_DOCS_README_HOME`
  is `1` (`0` disables it), so a tree without a landing page keeps the README as
  its home page. A heading of those files that opens with a link
  (a changelog's `## [0.3.2] - <date>`) keeps its text, and so its id, and the
  link moves to a line of its own under it: mdBook makes every heading a link,
  and a link inside a link has no accessible name;
- replaces a line `<!-- readme-section: <heading prefix> -->` with the README
  section whose `## ` heading starts with the prefix, from that heading to the
  line before the next `## ` heading or the end of the file; exactly one heading
  must match, or the build fails. `| noheading` drops the heading line, and
  `| fragments: <chapter>` points every bare `#fragment` link of the inserted
  section at that chapter (a path under the book source), for a section whose
  headings render on another page. A link of the inserted section labelled by
  its own repository path in backticks is labelled by the file's name instead,
  the target unchanged, so a list of records reads as file names;
- replaces `{{vep-rs-ref}}` with `$VEP_RS_DOCS_REF` (default `main`) and
  `{{vep-rs-release}}` with `$VEP_RS_DOCS_RELEASE`, else the most recent `v*`
  tag `git describe` finds from the checkout's HEAD, else a link to the latest
  release on GitHub (a checkout without tags, such as a shallow clone);
- resolves every relative markdown link against the file it was written in;
- keeps a link to a page in the book relative (mdBook turns `.md` into `.html`),
  points `README.md`, `CHANGELOG.md` and `CONTRIBUTING.md` at their pages when
  the book has them, and sends any other repository path to its GitHub view at
  the ref the site is built from.

Links, directives and tokens inside fenced code blocks are left alone. The
GitHub base URL comes from `output.html.git-repository-url` in `book.toml`; the
ref is `$VEP_RS_DOCS_REF` (a release tag, set by
scripts/docs/build_release_site.sh), else `main`.

mdBook runs it (`[preprocessor.repo-links]` in `book.toml`): `supports <renderer>`
exits 0, and otherwise it reads `[context, book]` JSON on stdin and writes the
book back on stdout. Stdlib only.
"""

import dataclasses
import functools
import json
import os
import re
import subprocess
import sys

LINK = re.compile(r"(\]\()([^)\s]+)((?:\s+\"[^\"]*\")?\))")
FENCE = re.compile(r"^\s*(```|~~~)")
SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
DIRECTIVE = re.compile(r"^\s*<!--\s*readme-section:\s*(.*?)\s*-->\s*$")
SELF_NAMED_LINK = re.compile(r"\[`([^`\s]+)`\](?=\(([^)\s#]+))")
HEADING_LINK = re.compile(r"^(#{1,6} +)\[([^\]]+)\](?!\[)(?:\(([^)\s]+)\))?(.*?)\s*$")
REFERENCE = re.compile(r"^ {0,3}\[([^\]]+)\]:\s*(\S+)", re.MULTILINE)
# Chapters whose content is a repository file outside the book source.
REPO_FILE_CHAPTERS = {"changelog.md": "CHANGELOG.md", "contributing.md": "CONTRIBUTING.md"}
STUB_HEADING = "# vep-rs documentation"


@dataclasses.dataclass(frozen=True)
class Site:
    """The book's layout and the build's variables, shared by every chapter's rewrite."""

    root: str
    src: str
    chapters: frozenset
    blob_url: str
    ref: str
    chapter_of_repo_file: dict
    tokens: dict


@functools.cache
def read_repo_file(root, name):
    try:
        with open(os.path.join(root, name), encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        sys.exit(f"mdbook_repo_links: {name} is missing from {root}")


def link_label(text, url):
    """Name a link moved out of a heading by the GitHub page it opens, else by its text."""
    for marker, label in (("/compare/", "Compare"), ("/releases/tag/", "Release")):
        if marker in url:
            return f"{label} {url.split(marker, 1)[1]}"
    return text


def unlink_headings(content):
    """Move the link that opens a heading to a line of its own under the heading."""
    references = {label.lower(): url for label, url in REFERENCE.findall(content)}
    out, in_fence = [], False
    for line in content.splitlines(keepends=True):
        if FENCE.match(line):
            in_fence = not in_fence
        heading = None if in_fence else HEADING_LINK.match(line)
        url = heading and (heading[3] or references.get(heading[2].lower()))
        if url:
            hashes, text, _, rest = heading.groups()
            line = f"{hashes}{text}{rest}\n\n[{link_label(text, url)}]({url})\n\n"
        out.append(line)
    return "".join(out)


def latest_release(root):
    """Return the most recent `v*` tag reachable from the checkout's HEAD, or None."""
    try:
        described = subprocess.run(
            ["git", "-C", root, "describe", "--tags", "--abbrev=0", "--match", "v[0-9]*"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return described.stdout.strip() or None


def rewrite_target(target, base_dir, chapter_dir, site, fragment_page=None):
    """Return the site-correct form of one link target.

    base_dir is the repository-relative directory the link was written in;
    chapter_dir is the chapter's directory relative to the book source;
    fragment_page, when set, is the page a bare `#fragment` link points at.
    """
    if target.startswith("#"):
        return fragment_page + target if fragment_page else target
    if SCHEME.match(target):
        return target
    path, sep, frag = target.partition("#")
    resolved = os.path.normpath(os.path.join(base_dir, path)).replace(os.sep, "/")
    if resolved.startswith("../"):
        return target
    in_book = None
    if site.chapter_of_repo_file.get(resolved) in site.chapters:
        in_book = site.chapter_of_repo_file[resolved]
    elif resolved.startswith(site.src + "/") and resolved[len(site.src) + 1 :] in site.chapters:
        in_book = resolved[len(site.src) + 1 :]
    if in_book is not None:
        rel = os.path.relpath(in_book, chapter_dir or ".").replace(os.sep, "/")
        return rel + sep + frag
    kind = "tree" if os.path.isdir(os.path.join(site.root, resolved)) else "blob"
    return f"{site.blob_url}/{kind}/{site.ref}/{resolved}" + sep + frag


def parse_directive(spec):
    """Return (heading prefix, fragments chapter or None, keep heading) of a directive body."""
    prefix, *options = (part.strip() for part in spec.split("|"))
    fragments, keep_heading = None, True
    for option in options:
        if option == "noheading":
            keep_heading = False
        elif option.startswith("fragments:") and option[len("fragments:") :].strip():
            fragments = option[len("fragments:") :].strip()
        else:
            sys.exit(
                f"mdbook_repo_links: readme-section option {option!r} is not"
                " `noheading` or `fragments: <chapter>`"
            )
    if not prefix:
        sys.exit("mdbook_repo_links: readme-section needs a heading prefix")
    return prefix, fragments, keep_heading


def select_section(lines, prefix):
    """Return the README lines from the one `## ` heading starting with prefix to the next `## `."""
    headings, in_fence = [], False
    for i, line in enumerate(lines):
        if FENCE.match(line):
            in_fence = not in_fence
        elif not in_fence and line.startswith("## "):
            headings.append(i)
    matches = [i for i in headings if lines[i][3:].lstrip().startswith(prefix)]
    if len(matches) != 1:
        found = ", ".join(lines[i].rstrip() for i in matches) or "none"
        sys.exit(
            f"mdbook_repo_links: readme-section {prefix!r} must match one `## ` heading"
            f" of README.md; matched: {found}"
        )
    following = [i for i in headings if i > matches[0]]
    return lines[matches[0] : following[0] if following else len(lines)]


def readme_section(spec, chapter_dir, site):
    """Return the README section a readme-section directive names, rewritten for the chapter."""
    prefix, fragments, keep_heading = parse_directive(spec)
    readme = read_repo_file(site.root, "README.md").splitlines(keepends=True)
    lines = select_section(readme, prefix)
    if not keep_heading:
        lines = lines[1:]
    fragment_page = None
    if fragments:
        if fragments not in site.chapters:
            sys.exit(f"mdbook_repo_links: readme-section fragments: {fragments!r} is not a chapter")
        fragment_page = os.path.relpath(fragments, chapter_dir or ".").replace(os.sep, "/")
    text = "".join(lines)
    if not text.endswith("\n"):
        text += "\n"
    return rewrite(text, "", chapter_dir, site, fragment_page, included=True)


def name_by_file(m):
    """Label a link written as [`<path>`](<path>) with the path's last component."""
    text, path = m[1], m[2]
    if text == path and "/" in text:
        return f"[`{os.path.basename(text)}`]"
    return m[0]


def rewrite(content, base_dir, chapter_dir, site, fragment_page=None, included=False):
    out, in_fence = [], False
    for line in content.splitlines(keepends=True):
        directive = None if in_fence else DIRECTIVE.match(line)
        if FENCE.match(line):
            in_fence = not in_fence
        elif directive:
            line = readme_section(directive[1], chapter_dir, site)
        elif not in_fence:
            if included:
                line = SELF_NAMED_LINK.sub(name_by_file, line)
            line = LINK.sub(
                lambda m: m[1]
                + rewrite_target(m[2], base_dir, chapter_dir, site, fragment_page)
                + m[3],
                line,
            )
            for token, value in site.tokens.items():
                line = line.replace(token, value)
        out.append(line)
    return "".join(out)


def walk(items):
    for item in items:
        chapter = item.get("Chapter") if isinstance(item, dict) else None
        if chapter:
            yield chapter
            yield from walk(chapter.get("sub_items", []))


def readme_is_home(book):
    """Whether the `index.md` chapter shows README.md: the variable decides, else the stub does."""
    setting = os.environ.get("VEP_RS_DOCS_README_HOME")
    if setting in ("0", "1"):
        return setting == "1"
    for chapter in walk(book["items"]):
        if chapter.get("source_path") == "index.md":
            # Lines as build_release_site.sh counts them with awk: a final newline closes
            # the last line, and nothing is stripped, so both tests read a file alike.
            lines = chapter["content"].split("\n")
            if lines[-1] == "":
                lines.pop()
            return len(lines) == 3 and lines[0] == STUB_HEADING
    return False


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "supports":
        sys.exit(0)
    context, book = json.load(sys.stdin)
    root = context["root"]
    config = context["config"]
    src = os.path.normpath(config.get("book", {}).get("src", "src")).replace(os.sep, "/")
    repo_url = config.get("output", {}).get("html", {}).get("git-repository-url")
    if not repo_url:
        sys.exit("mdbook_repo_links: output.html.git-repository-url is not set in book.toml")
    repo_file_chapters = dict(REPO_FILE_CHAPTERS)
    if readme_is_home(book):
        repo_file_chapters["index.md"] = "README.md"
    ref = os.environ.get("VEP_RS_DOCS_REF") or "main"
    blob_url = repo_url.rstrip("/")
    release = (
        os.environ.get("VEP_RS_DOCS_RELEASE")
        or latest_release(root)
        or f"[on GitHub]({blob_url}/releases/latest)"
    )
    site = Site(
        root=root,
        src=src,
        chapters=frozenset(c["source_path"] for c in walk(book["items"]) if c.get("source_path")),
        blob_url=blob_url,
        ref=ref,
        chapter_of_repo_file={v: k for k, v in repo_file_chapters.items()},
        tokens={"{{vep-rs-ref}}": ref, "{{vep-rs-release}}": release},
    )

    for chapter in walk(book["items"]):
        source = chapter.get("source_path")
        if not source:
            continue
        chapter_dir = os.path.dirname(source)
        if source in repo_file_chapters:
            chapter["content"] = unlink_headings(read_repo_file(root, repo_file_chapters[source]))
            base_dir = ""
        else:
            base_dir = os.path.join(src, chapter_dir) if chapter_dir else src
        chapter["content"] = rewrite(chapter["content"], base_dir, chapter_dir, site)
    json.dump(book, sys.stdout)


if __name__ == "__main__":
    main()
