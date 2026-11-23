#!/usr/bin/env python3
"""Finish one channel of the built documentation site and check what the build cannot.

Runs after `mdbook build` on the release channel at the site root or on the
development preview of main under `main/`. On every page it rewrites the menu
bar's `<h1 class="menu-title">` as a `<p>`, so the page's own heading is its
one `h1`, and collapses a `<title>` whose chapter name repeats the book title.
On the release channel it drops the release tag from every `<title>` and from
`og:site_name`, so a title names the project without a version; the preview
keeps `vep-rs main`.
On every page except `print.html`, `toc.html`, `404.html` and the stub mdBook
writes for a redirected path, it inserts `<link rel="canonical">` and `og:url`
for the page's address under `--base` and, on the preview,
`<meta name="vep-rs-channel" content="main">` and a `noindex` robots meta, so a
search engine lists the release and not the preview. On those pages it also
replaces the book-wide description metas with the page's first paragraph of
prose, so each search result shows its own snippet, and on the home page it
names what the project is in the title and, for a release, adds a schema.org
SoftwareApplication JSON-LD record. On the release channel it fits the title
of each of those pages to 25 to 60 characters: a title under 25 ends
` - vep-rs documentation` in place of ` - vep-rs`, a title over 60 is the
chapter name alone, and a name over 60 is cut at a word and ends in an
ellipsis. The home page's title has no such ending and is left as it is, and
`og:title` keeps every chapter name whole. In the home page's `<noscript>`
block it puts the chapter list of `toc.html` where mdBook frames that file, so
the page carries every chapter as a link without a script. It writes
`versions.json`, which the header's channel switcher reads, and for a release
`sitemap.xml`.

When `VEP_RS_DOCS_README_HOME` is `1` the home page renders the README, so its
edit link is pointed at `README.md`.

Then two checks. Every page a reader lands on has a title and a description
of its own: a title or a description two pages share, or a page left on the
book's description, fails the preview and is a warning line on the release,
whose pages come from a published tag; the release's sitemap lists every such
page. On the release a title outside 25 to 60 characters fails, and so does a
title that chapters of different names were fitted to. Every redirect
destination the built pages name must exist:
the page as a file and, when the destination carries one, the id on that page.
Every step prints one line, the page check a second one for the release's
titles; any failure exits 1 after all of them have run.

Usage::

    python3 scripts/docs/site_postbuild.py <site> --base <url> --release <tag>
        [--channel main]

Stdlib only.
"""

import argparse
import html
import html.parser
import json
import os
import posixpath
import re
import sys
import urllib.parse

# Pages mdBook generates whose address is not a reader's destination.
SKIP_PAGES = {"print.html", "toc.html", "404.html"}
PREVIEW_CHANNEL_NAME = "main"
MENU_TITLE = re.compile(r'<h1 class="menu-title">(.*?)</h1>', re.DOTALL)
TITLE = re.compile(r"<title>(.*?) - (.*?)</title>")
PAGE_TITLE = re.compile(r"(<title>)([^<]*)(</title>)")
SITE_NAME = re.compile(r'(<meta property="og:site_name" content=")([^"]*)(">)')
EDIT_HOME = re.compile(r'(href="[^"]*/edit/[^"]*/)[^"/]+/index\.md"')
REFRESH = re.compile(r'http-equiv="refresh"\s+content="0;\s*URL=([^"]+)"', re.IGNORECASE)
FRAGMENT_MAP = re.compile(r"fragmentMap\s*=\s*(\{.*?\})\s*;", re.DOTALL)
REDIRECT_URL = re.compile(r'let redirectUrl = "([^"]+)"')
SCHEME = re.compile(r"^[a-zA-Z][a-zA-Z0-9+.-]*:")
PARAGRAPH = re.compile(r"<p>(.*?)</p>", re.DOTALL)
TAG = re.compile(r"<[^>]+>")
DESCRIPTION_METAS = re.compile(
    r'(<meta (?:name="description"|property="og:description"|name="twitter:description") content=")[^"]*(">)'
)
CHAPTER_LIST = re.compile(r"<body[^>]*>\s*(.*?)\s*</body>", re.DOTALL)
CHAPTER_FRAME = re.compile(
    r'(<noscript>\s*)<iframe class="sidebar-iframe-outer" src="toc\.html"></iframe>(\s*</noscript>)'
)
HOME_TITLE = re.compile(r"<title>([^<]*)</title>")
HOME_OG_TITLES = re.compile(r'(<meta (?:property="og:title"|name="twitter:title") content=")[^"]*(">)')
OG_TITLE = re.compile(r'<meta property="og:title" content="([^"]*)">')
# A search result shows the title and description; these say what the project is in the words people search.
HOME_TAGLINE = "fast Ensembl VEP (Variant Effect Predictor) in Rust"
DESCRIPTION_MIN, DESCRIPTION_MAX = 80, 160
BOOK_NAME = "vep-rs"
# The ending of a chapter's title on the release channel, and the ending of a title that is too short with it.
TITLE_ENDING, TITLE_ENDING_LONG = f" - {BOOK_NAME}", f" - {BOOK_NAME} documentation"
# The lengths in characters a title of the release channel is fitted to.
TITLE_MIN, TITLE_MAX = 25, 60
ELLIPSIS = "\u2026"
# A title cut at a word does not end on an article, a conjunction or one of these prepositions.
SMALL_WORDS = frozenset(
    {"a", "an", "the", "and", "but", "nor", "or"}
    | {"as", "at", "before", "by", "for", "from", "in", "into", "of", "on", "over", "to", "under", "with", "without"}
)


class Ids(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name == "id" and value:
                self.ids.add(value)


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def write(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def site_pages(site):
    pages = []
    for dirpath, _, files in os.walk(site):
        for name in files:
            if name.endswith(".html"):
                rel = os.path.relpath(os.path.join(dirpath, name), site)
                pages.append(rel.replace(os.sep, "/"))
    return sorted(pages)


def page_url(base, page):
    if posixpath.basename(page) == "index.html":
        directory = posixpath.dirname(page)
        return base + (directory + "/" if directory else "")
    return base + page


def collapse_title(m):
    """`<chapter> - <book>` becomes `<book>` when the chapter is the book's name or its first word."""
    chapter, book = m[1], m[2]
    return f"<title>{book}</title>" if chapter in (book, book.split(" ")[0]) else m[0]


def drop_release(text, release):
    """The page's `<title>` and `og:site_name` without the release tag the release's book title ends with.

    The title loses the tag only after the book's name, so a chapter name that ends in the tag keeps it.
    """
    tag = " " + release

    def bare(m):
        return m[1] + m[2].removesuffix(tag) + m[3]

    def bare_title(m):
        title = m[2].removesuffix(tag)
        return m[1] + (title if title == BOOK_NAME or title.endswith(TITLE_ENDING) else m[2]) + m[3]

    return SITE_NAME.sub(bare, PAGE_TITLE.sub(bare_title, text, 1))


def fitted_title(heading):
    """The release channel's title for the chapter named `heading`.

    The name with TITLE_ENDING where that is TITLE_MIN to TITLE_MAX characters, with TITLE_ENDING_LONG where it
    is shorter, and the name alone where it is longer. A name that is itself over TITLE_MAX characters is cut
    and ends in an ellipsis: at the last space that leaves a title of those lengths, less a trailing `,`, `;`
    or `:` and trailing SMALL_WORDS while TITLE_MIN characters remain, or inside a word where it has no such
    space. A name of one character gives a title one short of TITLE_MIN.
    """
    title = heading + TITLE_ENDING
    if len(title) < TITLE_MIN:
        return heading + TITLE_ENDING_LONG
    if len(title) <= TITLE_MAX:
        return title
    if len(heading) <= TITLE_MAX:
        return heading
    room = TITLE_MAX - len(ELLIPSIS)
    space = heading.rfind(" ", TITLE_MIN - len(ELLIPSIS), room + 1)
    if space < 0:
        return heading[:room] + ELLIPSIS
    words = heading[:space]
    while True:
        last = words.rpartition(" ")[2]
        if words[-1] in ",;: ":
            shorter = words.rstrip(",;: ")
        elif last != words and last.lower() in SMALL_WORDS:
            shorter = words[: -len(last)].rstrip()
        else:
            break
        if len(shorter) + len(ELLIPSIS) < TITLE_MIN:
            break
        words = shorter
    return words + ELLIPSIS


def title_form(title):
    """Which form of `fitted_title` a title has, read from its ending: "kept", "lengthened", "cut" or "bare"."""
    if title.endswith(TITLE_ENDING_LONG):
        return "lengthened"
    if title.endswith(TITLE_ENDING):
        return "kept"
    return "cut" if title.endswith(ELLIPSIS) else "bare"


def fit_title(text):
    """The page with its `<title>` fitted by `fitted_title`, the title measured unescaped.

    A title that does not end in TITLE_ENDING is left as it is: the home page's, and a title a finished page
    already carries in another form. A title the fit does not change keeps its bytes.
    """
    m = PAGE_TITLE.search(text)
    title = html.unescape(m[2]) if m else ""
    if not title.endswith(TITLE_ENDING):
        return text
    fitted = fitted_title(title.removesuffix(TITLE_ENDING))
    if fitted == title:
        return text
    return text[: m.start(2)] + html.escape(fitted, quote=False) + text[m.end(2) :]


def page_heading(text, title):
    """The name a page's title was made from: the chapter name in its `og:title`, which the fit leaves whole,
    else the title without its ending. The home page's `og:title` is its title."""
    m = OG_TITLE.search(text)
    if m:
        return html.unescape(m[1])
    return title.removesuffix(TITLE_ENDING_LONG).removesuffix(TITLE_ENDING)


def page_description(text):
    """The page's first paragraph of prose in its main content, cut to a search snippet's length; None if none."""
    start = text.find("<main>")
    for m in PARAGRAPH.finditer(text, max(start, 0)):
        prose = " ".join(html.unescape(TAG.sub("", m[1])).split())
        if len(prose) < DESCRIPTION_MIN:
            continue
        if len(prose) > DESCRIPTION_MAX:
            prose = prose[: DESCRIPTION_MAX - 1].rsplit(" ", 1)[0].rstrip(",;:") + "\u2026"
        return prose
    return None


def describe_page(text):
    """Replaces the book-wide description metas with the page's own; a page without prose keeps the book's."""
    description = page_description(text)
    if not description:
        return text
    value = html.escape(description)
    return DESCRIPTION_METAS.sub(lambda m: m[1] + value + m[2], text)


def title_home(text, release, channel, base):
    """Names the project and what it is in the home page's title, and on the release adds its JSON-LD record.

    The title is the book's name and HOME_TAGLINE, without the chapter name a home chapter of another name
    puts before the book's.
    """
    m = HOME_TITLE.search(text)
    if m and HOME_TAGLINE not in m[1]:
        title = f"{m[1].rpartition(' - ')[2]}: {HOME_TAGLINE}"
        text = text.replace(m[0], f"<title>{html.escape(title, quote=False)}</title>", 1)
        text = HOME_OG_TITLES.sub(lambda t: t[1] + html.escape(title) + t[2], text)
    if not channel and "application/ld+json" not in text and "</head>" in text:
        record = {
            "@context": "https://schema.org",
            "@type": "SoftwareApplication",
            "name": "vep-rs",
            "alternateName": "vep-rs: Ensembl VEP in Rust",
            "description": page_description(text) or HOME_TAGLINE,
            "applicationCategory": "Bioinformatics",
            "operatingSystem": "Linux, macOS",
            "programmingLanguage": "Rust",
            "softwareVersion": release.lstrip("v"),
            "license": "https://www.apache.org/licenses/LICENSE-2.0",
            "url": base,
            "codeRepository": "https://github.com/natera-open-source/vep-rs",
            "downloadUrl": "https://github.com/natera-open-source/vep-rs/releases",
            "offers": {"@type": "Offer", "price": "0", "priceCurrency": "USD"},
        }
        script = '<script type="application/ld+json">' + json.dumps(record).replace("</", "<\\/") + "</script>"
        text = text.replace("</head>", f"    {script}\n    </head>", 1)
    return text


def chapter_list(site):
    """The list of every chapter that mdBook writes to `toc.html` at the site root; None if the site has none."""
    path = os.path.join(site, "toc.html")
    m = CHAPTER_LIST.search(read(path)) if os.path.isfile(path) else None
    return m[1] if m else None


def list_chapters(text, chapters):
    """Puts the chapter list where a `<noscript>` block frames `toc.html`; the page and the count of links listed.

    The sidebar does not scroll, so the list sits in a box of the class that does. Its links are relative to the
    site root, so the page is the home page.
    """
    box = f'<div class="sidebar-scrollbox">{chapters}</div>'
    text, framed = CHAPTER_FRAME.subn(lambda m: m[1] + box + m[2], text, 1)
    return text, chapters.count(' href="') if framed else 0


def finish_pages(site, pages, reader_pages, base, channel, readme_home, release):
    problems, inserted, present = [], 0, 0
    chapters, listed = chapter_list(site), 0
    for page in pages:
        path = os.path.join(site, page)
        text = MENU_TITLE.sub(r'<p class="menu-title">\1</p>', read(path))
        text = TITLE.sub(collapse_title, text, 1)
        if not channel:
            text = drop_release(text, release)
        if readme_home and page == "index.html":
            text = EDIT_HOME.sub(r'\1README.md"', text)
        if page in reader_pages:
            text = describe_page(text)
            if page == "index.html":
                text = title_home(text, release, channel, base)
                if chapters:
                    text, listed = list_chapters(text, chapters)
            if not channel:
                text = fit_title(text)
            if 'rel="canonical"' in text:
                present += 1
            elif "</head>" not in text:
                problems.append(f"{page}: no </head> to insert the canonical link before")
            else:
                url = html.escape(page_url(base, page))
                tags = [f'<link rel="canonical" href="{url}">', f'<meta property="og:url" content="{url}">']
                if channel:
                    tags.append(f'<meta name="vep-rs-channel" content="{channel}">')
                    tags.append('<meta name="robots" content="noindex">')
                text = text.replace("</head>", "\n        ".join(tags) + "\n    </head>", 1)
                inserted += 1
        write(path, text)
    what = "canonical, og:url" + (", channel, noindex" if channel else "")
    print(
        f"{len(pages)} pages finished: menu title demoted, head tags ({what}) inserted on {inserted},"
        f" {present} already carried a canonical link"
    )
    if listed:
        print(f"home page: {listed} chapter links listed in its <noscript> block")
    return problems


def write_versions(site, base, release, channel):
    root = urllib.parse.urlsplit(base).path
    if channel:
        suffix = channel + "/"
        if not root.endswith(suffix):
            return [f"--base path {root} does not end in {suffix}, so the release path cannot be derived"]
        root = root[: -len(suffix)]
    preview_path = root + PREVIEW_CHANNEL_NAME + "/"
    versions = {
        "release": release,
        "channels": [
            {"name": release, "path": root},
            {"name": PREVIEW_CHANNEL_NAME, "path": preview_path},
        ],
    }
    write(os.path.join(site, "versions.json"), json.dumps(versions, indent=2) + "\n")
    print(f"versions.json: release {release}, channels at {root} and {preview_path}")
    return []


def write_sitemap(site, pages, base):
    lines = ['<?xml version="1.0" encoding="UTF-8"?>']
    lines.append('<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">')
    for page in pages:
        lines.append(f"  <url><loc>{html.escape(page_url(base, page))}</loc></url>")
    lines.append("</urlset>")
    write(os.path.join(site, "sitemap.xml"), "\n".join(lines) + "\n")
    print(f"sitemap.xml: {len(pages)} pages")


def shared_heads(what, heads):
    """One finding for each title or description that more than one page has."""
    return [
        f"{what} {head!r} is on {len(pages)} pages: {', '.join(pages)}"
        for head, pages in heads.items()
        if len(pages) > 1
    ]


def check_page_heads(site, reader_pages, channel):
    """Every page a reader lands on has a title and a description of its own, and the release's sitemap lists it.

    A title or a description two pages share, or a page left on the book's description, fails the preview,
    whose pages can still be changed; the release's pages come from a published tag, so each is a warning there.
    The release's titles are fitted, and two outcomes of the fit fail there: one title for chapters of
    different names, and a title outside TITLE_MIN to TITLE_MAX characters.
    """
    titles, headings, descriptions, on_book = {}, {}, {}, []
    for page in sorted(reader_pages):
        text = read(os.path.join(site, page))
        m = PAGE_TITLE.search(text)
        title = html.unescape(m[2]) if m else ""
        titles.setdefault(title, []).append(page)
        headings.setdefault(title, set()).add(page_heading(text, title))
        description = page_description(text)
        if description:
            descriptions.setdefault(description, []).append(page)
        else:
            on_book.append(page)
    fitted_alike = {} if channel else {title: pages for title, pages in titles.items() if len(headings[title]) > 1}
    shared_titles = shared_heads("title", {title: titles[title] for title in titles if title not in fitted_alike})
    shared_descriptions = shared_heads("description", descriptions)
    print(
        f"page heads: {len(reader_pages)} pages, {len(shared_titles) + len(fitted_alike)} titles shared,"
        f" {len(shared_descriptions)} descriptions shared, {len(on_book)} on the book's description"
    )
    findings = shared_titles + shared_descriptions
    findings += [
        f"{page}: no paragraph of {DESCRIPTION_MIN} characters or more, so its description is the book's"
        for page in on_book
    ]
    if channel:
        return findings
    for finding in findings:
        print(f"WARNING: [site_postbuild] {finding}")
    failures = [
        f"title {title!r} is the fit of {len(headings[title])} headings: {', '.join(pages)}"
        for title, pages in fitted_alike.items()
    ]
    outside = []
    for title, pages in titles.items():
        if not TITLE_MIN <= len(title) <= TITLE_MAX:
            limit = f"under {TITLE_MIN}" if len(title) < TITLE_MIN else f"over {TITLE_MAX}"
            outside += [f"{page}: title {title!r} is {len(title)} characters, {limit}" for page in pages]
    forms = [title_form(title) for title, pages in titles.items() for _ in pages]
    kept, lengthened, bare, cut = (forms.count(form) for form in ("kept", "lengthened", "bare", "cut"))
    print(
        f'page titles: {len(reader_pages)} pages, {kept} end "{TITLE_ENDING}",'
        f' {lengthened} end "{TITLE_ENDING_LONG}", {bare} without an ending, {cut} cut,'
        f" {len(outside)} outside {TITLE_MIN} to {TITLE_MAX} characters"
    )
    failures += outside
    sitemap = os.path.join(site, "sitemap.xml")
    listed = read(sitemap).count("<loc>") if os.path.isfile(sitemap) else 0
    if listed != len(reader_pages):
        failures.append(f"sitemap.xml lists {listed} pages, a reader lands on {len(reader_pages)}")
    return failures


def redirect_destinations(text):
    """Every destination a page's redirect script or refresh tag names, in page order."""
    destinations = []
    for m in FRAGMENT_MAP.finditer(text):
        destinations.extend(json.loads(m[1]).values())
    destinations.extend(m[1] for m in REDIRECT_URL.finditer(text))
    destinations.extend(m[1] for m in REFRESH.finditer(text))
    return list(dict.fromkeys(destinations))


def check_redirects(site, pages, base):
    base_path = urllib.parse.urlsplit(base).path
    ids = {}
    broken, checked, external = [], 0, 0
    for page in pages:
        for dest in redirect_destinations(read(os.path.join(site, page))):
            if SCHEME.match(dest):
                external += 1
                continue
            checked += 1
            path, _, fragment = dest.partition("#")
            if path.endswith("/"):
                path += "index.html"
            if not path:
                target = page
            elif path.startswith("/"):
                if not path.startswith(base_path) and path != base_path.rstrip("/"):
                    broken.append(f"{page}: redirect to {dest} leaves the site")
                    continue
                target = path[len(base_path):]
            else:
                target = posixpath.normpath(posixpath.join(posixpath.dirname(page), path))
            if not target:
                target = "index.html"
            if not os.path.isfile(os.path.join(site, target)):
                broken.append(f"{page}: redirect to {dest}: no such file in the site")
                continue
            if fragment and target.endswith(".html"):
                if target not in ids:
                    parser = Ids()
                    parser.feed(read(os.path.join(site, target)))
                    ids[target] = parser.ids
                if fragment not in ids[target]:
                    broken.append(f"{page}: redirect to {dest}: no id {fragment!r} on {target}")
    print(f"redirect destinations: {checked} checked, {external} external not fetched, {len(broken)} broken")
    return broken


def main():
    ap = argparse.ArgumentParser(description="Finish one channel of the built documentation site and check it.")
    ap.add_argument("site", help="the built site, mdBook's output directory")
    ap.add_argument("--base", required=True, help="the URL the site is served at, ending in a slash")
    ap.add_argument("--release", required=True, help="the release tag the site root serves (vX.Y.Z)")
    ap.add_argument("--channel", help="the channel built under <base>; omitted for the release at the root")
    args = ap.parse_args()
    base = args.base if args.base.endswith("/") else args.base + "/"
    readme_home = os.environ.get("VEP_RS_DOCS_README_HOME") == "1"

    pages = site_pages(args.site)
    if not pages:
        sys.exit(f"ERROR: [site_postbuild] no .html page under {args.site}")
    reader_pages = {
        p for p in pages if p not in SKIP_PAGES and not REFRESH.search(read(os.path.join(args.site, p)))
    }

    failures = []
    failures += finish_pages(args.site, pages, reader_pages, base, args.channel, readme_home, args.release)
    failures += write_versions(args.site, base, args.release, args.channel)
    if not args.channel:
        write_sitemap(args.site, sorted(reader_pages), base)
    failures += check_page_heads(args.site, reader_pages, args.channel)
    failures += check_redirects(args.site, pages, base)

    for failure in failures:
        print(f"ERROR: [site_postbuild] {failure}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
