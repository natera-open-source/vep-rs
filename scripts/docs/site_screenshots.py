#!/usr/bin/env python3
"""Screenshot and probe the built documentation site from phone to desktop widths.

This is the one script under `scripts/` with a third-party runtime: it needs
Playwright (`pip install playwright`) and the Chromium it downloads, and it
never runs in CI. It serves the built site over HTTP at the path the site is
published under (`/vep-rs/` for a release, `/vep-rs/main/` for the preview), so
the 404 page's assets, relative links, the fetch of `versions.json` and the
version switcher behave as on the published site, opens every page at
every width, writes viewport-sized screenshots (never a full-page image) and
records a layout probe per page and width in `probe.json`: horizontal overflow
of the document and of the content column, every element whose right edge
leaves the viewport, inline code broken across lines (split into spans that
fit their block on one line and long spans wider than it), the smallest font size
on visible text, the characters per line of the text measure (by a ruler of
100 `n` characters and by the first paragraph's own text), and every code
block or table wrapper that scrolls, with or without a visible scrollbar.
`--axe` adds an axe-core run (WCAG 2.0 to 2.2, levels A and AA) through
`axe-playwright-python` when that package is importable. `contact-sheet.html`
beside the screenshots shows every page's tiles with the probe summary under
each and opens from a file URL. `summary.md` beside it holds one Markdown table
row per page and width: overflow, wrapped inline code, long-span wraps,
scrollers without a visible scrollbar, minimum font size and characters per
line, then the totals. `--summary` prints that table after a run, or, given
`probe.json` files, prints theirs without shooting and without Playwright.

Usage::

    python3 scripts/docs/site_screenshots.py --site target/book --out /tmp/shots \\
        [--pages index,cli-reference] [--widths 375,1280] [--themes light,navy] [--axe]
    python3 scripts/docs/site_screenshots.py --summary /tmp/shots/probe.json

In the Playwright container image, with the site and this directory mounted
(the Python package version must match the image's browser build)::

    docker run --rm -v "$PWD/target/book:/site:ro" -v "$PWD/scripts/docs:/scripts:ro" \\
        -v /tmp/shots:/out mcr.microsoft.com/playwright/python:v1.63.0-noble \\
        bash -c 'pip install -q playwright==1.63.0 axe-playwright-python && \\
                 python3 -I /scripts/site_screenshots.py --site /site --out /out --axe'

Screenshots are named `<page>-<width>-<theme>-<tile>.png`: `<page>` is the page
path with `/` replaced by `-`, `<tile>` is `top`, then `s1`, `s2` and on, each
starting where the screen before it ended. At the `--to-end-widths` (default
375) the tiles run to the foot of the page, at the other full widths to `s2`.
A last tile that would repeat rows already shot shows only the rows below them.
Scrollbars are classic, as on Linux and Windows, so a scrolling block's
scrollbar height is measurable.
"""

import argparse
import functools
import html
import http.server
import json
import os
import re
import sys
import threading

try:
    from axe_playwright_python.sync_playwright import Axe
except ImportError:
    Axe = None

WIDTHS = [320, 360, 375, 390, 414, 768, 820, 1024, 1280, 1366, 1440, 1536, 1920, 2560]
FULL_WIDTHS = [375, 768, 1024, 1280, 1440, 1920]
AXE_WIDTHS = [375, 1280]
# The usual device or window height at each width.
HEIGHTS = {320: 568, 360: 800, 375: 667, 390: 844, 414: 896, 768: 1024, 820: 1180, 1024: 768,
           1280: 800, 1366: 768, 1440: 900, 1536: 864, 1920: 1080, 2560: 1440}
AXE_TAGS = ["wcag2a", "wcag2aa", "wcag21a", "wcag21aa", "wcag22aa"]
SKIP_PAGES = {"print.html", "toc.html", "404.html"}
TO_END_WIDTHS = [375]
# A remainder shorter than a line of text holds only the page's bottom margin.
MIN_NEW_ROWS = 24

# mdBook's menu bar can stay pinned over the top of a scrolled screen. Each
# tile after the first starts that many rows higher, so every row shows in some tile.
COVERED_JS = r"""() => {
  const bar = document.getElementById('mdbook-menu-bar');
  const r = bar ? bar.getBoundingClientRect() : null;
  return r && r.top <= 1 && r.bottom > 0 ? Math.ceil(r.bottom) : 0;
}"""

# The probe runs in the page. `tol` is the tolerance in pixels. Overflow is
# checked on the document and on the content column, because the theme makes
# the content column its own scroll container and a wide element inside it
# never widens the document.
PROBE_JS = r"""(tol) => {
  const doc = document.documentElement, vw = doc.clientWidth;
  const content = document.getElementById('mdbook-content') || document.body;
  const skip = new Set(['SCRIPT', 'STYLE', 'NOSCRIPT', 'TEMPLATE', 'LINK', 'META']);
  const shown = (el) => !skip.has(el.tagName) &&
    el.checkVisibility({visibilityProperty: true, opacityProperty: true});
  const text = (el, n) => (el.textContent || '').replace(/\s+/g, ' ').trim().slice(0, n);
  const describe = (el) => ({tag: el.tagName.toLowerCase(),
    class: typeof el.className === 'string' ? el.className : '', text: text(el, 60)});

  const over = [];
  for (const el of document.body.querySelectorAll('*')) {
    const box = el.closest('pre, .table-wrapper');
    if ((box && box !== el) || !shown(el)) continue;
    const r = el.getBoundingClientRect();
    if (r.width + r.height > 0 && r.right > vw + tol) over.push(el);
  }
  const overflowing = over.filter((el) => !over.some((o) => o !== el && el.contains(o)))
    .map((el) => Object.assign(describe(el), {right: Math.round(el.getBoundingClientRect().right)}));

  // A span on two or more lines is `wrapped` when its one-line width fits the
  // content box of its block, and a long-span wrap when it is wider than that box.
  const blockOf = (el) => {
    let p = el.parentElement;
    while (p && ['inline', 'contents'].includes(getComputedStyle(p).display)) p = p.parentElement;
    return p || document.body;
  };
  const inner = (el) => {
    const cs = getComputedStyle(el);
    return el.clientWidth - parseFloat(cs.paddingLeft) - parseFloat(cs.paddingRight);
  };
  const oneLine = (el) => {
    const clone = el.cloneNode(true);
    clone.style.cssText = 'position:absolute;left:0;top:0;visibility:hidden;white-space:nowrap;width:max-content';
    el.parentElement.appendChild(clone);
    const width = clone.getBoundingClientRect().width;
    clone.remove();
    return width;
  };
  const wrapped = [], longWraps = [];
  for (const el of document.querySelectorAll('code')) {
    if (el.closest('pre') || !shown(el)) continue;
    const rects = [...el.getClientRects()].filter((r) => r.width > 0);
    const lines = new Set(rects.map((r) => Math.round(r.top))).size;
    if (lines < 2) continue;
    const block = blockOf(el);
    const entry = {text: text(el, 60), rects: rects.length, lines, width: Math.round(oneLine(el)),
      block: block.tagName.toLowerCase(), block_width: Math.round(inner(block))};
    (entry.width <= entry.block_width + tol ? wrapped : longWraps).push(entry);
  }

  const minFont = (root) => {
    let size = null, where = null;
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    for (let node = walker.nextNode(); node; node = walker.nextNode()) {
      const el = node.parentElement;
      if (!node.nodeValue.trim() || !el || !shown(el)) continue;
      const range = document.createRange();
      range.selectNodeContents(node);
      const r = range.getBoundingClientRect();
      if (r.width === 0 || r.height === 0 || r.right <= 0 || r.left >= vw) continue;
      const px = parseFloat(getComputedStyle(el).fontSize);
      if (size === null || px < size) { size = px; where = el; }
    }
    return {px: size, element: where ? describe(where) : null};
  };
  const contentFont = minFont(content), pageFont = minFont(document.body);

  const column = [...content.querySelectorAll('main p')].find((p) => shown(p) && text(p, 80).length > 40)
    || content.querySelector('main') || content;
  const glyph = (sample) => {
    const ruler = document.createElement('span');
    ruler.textContent = sample;
    ruler.style.cssText = 'position:absolute;left:0;top:0;visibility:hidden;white-space:nowrap;font:inherit';
    column.appendChild(ruler);
    const width = ruler.getBoundingClientRect().width / sample.length;
    ruler.remove();
    return width > 0 ? Math.round(column.clientWidth / width * 10) / 10 : null;
  };
  const sample = text(column, 200);

  const scrolls = (el) => el && el.clientWidth > 0 && el.scrollWidth > el.clientWidth + tol;
  const hasBar = (el) => {
    const cs = getComputedStyle(el);
    return el.offsetHeight - el.clientHeight - parseFloat(cs.borderTopWidth) -
      parseFloat(cs.borderBottomWidth) > 0;
  };
  const scroller = (kind, label, el) => ({kind, text: label, scroll_width: el.scrollWidth,
    client_width: el.clientWidth, scrollbar_visible: hasBar(el)});
  const scrollers = [];
  for (const pre of document.querySelectorAll('pre')) {
    const code = pre.querySelector('code');
    const el = scrolls(pre) ? pre : scrolls(code) ? code : null;
    if (el && shown(pre)) {
      scrollers.push(scroller('pre', (pre.textContent || '').trim().split('\n')[0].slice(0, 60), el));
    }
  }
  for (const w of document.querySelectorAll('.table-wrapper')) {
    if (shown(w) && scrolls(w)) scrollers.push(scroller('table', text(w.querySelector('tr') || w, 60), w));
  }

  return {
    title: document.title, document_height: doc.scrollHeight,
    client_width: vw, scroll_width: doc.scrollWidth,
    content_client_width: content.clientWidth, content_scroll_width: content.scrollWidth,
    overflow: doc.scrollWidth > vw + tol || content.scrollWidth > content.clientWidth + tol,
    overflowing_total: overflowing.length, overflowing: overflowing.slice(0, 40),
    wrapped_inline_code: wrapped, long_span_wraps: longWraps,
    min_font_px: contentFont.px, min_font_element: contentFont.element,
    min_font_page_px: pageFont.px, min_font_page_element: pageFont.element,
    measure_px: column.clientWidth,
    chars_per_line: glyph('n'.repeat(100)),
    chars_per_line_text: sample.length >= 40 ? glyph(sample) : null,
    scrollers,
  };
}"""


def site_base(site):
    # mdBook writes the book's site-url into 404.html alone, as its <base href>.
    try:
        with open(os.path.join(site, "404.html"), encoding="utf-8") as f:
            found = re.search(r'<base href="([^"]*)"', f.read())
    except OSError:
        found = None
    return found.group(1) if found else "/"


def serve(site, base):
    class Handler(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def translate_path(self, path):
            if not path.startswith(base):
                return ""
            return super().translate_path("/" + path[len(base):])

    handler = functools.partial(Handler, directory=site)
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}{base}"


def find_pages(site, wanted):
    if wanted and wanted != "all":
        names = [p.strip().removesuffix(".html") for p in wanted.split(",") if p.strip()]
    else:
        names = []
        for dirpath, _, files in os.walk(site):
            for name in files:
                rel = os.path.relpath(os.path.join(dirpath, name), site).replace(os.sep, "/")
                if name.endswith(".html") and rel not in SKIP_PAGES:
                    names.append(rel[:-5])
        names.sort()
    slugs = {}
    for name in names:
        if not os.path.isfile(os.path.join(site, name + ".html")):
            sys.exit(f"site_screenshots: no page {name}.html under {site}")
        slug = name.replace("/", "-")
        if slug in slugs:
            sys.exit(f"site_screenshots: pages {slugs[slug]} and {name} share the file name {slug}")
        slugs[slug] = name
    return names


def ints(text):
    return [int(v) for v in text.split(",") if v.strip()]


def axe_summary(violations):
    return {"count": len(violations), "violations": [
        {"id": v["id"], "impact": v.get("impact"), "help": v["help"], "nodes": len(v["nodes"]),
         "targets": [", ".join(map(str, n["target"])) for n in v["nodes"][:10]]} for v in violations]}


def tile_name(stem, tile):
    return f"{stem}-{tile}.png"


def tile_path(out, stem, tile):
    return os.path.join(out, tile_name(stem, tile))


def shoot_tiles(page, width, limit, out, stem):
    """Screenshot the page from the top down in at most `limit` tiles.

    Returns the tile names and whether rows remain below the last tile."""
    def snap(name, clip=None):
        page.wait_for_timeout(150)
        page.screenshot(path=tile_path(out, stem, name), clip=clip, animations="disabled")
        names.append(name)

    names = []
    page.evaluate("window.scrollTo(0, 0)")
    snap("top")
    total, view = page.evaluate("[document.documentElement.scrollHeight, document.documentElement.clientHeight]")
    shown = view
    while shown < total and len(names) < limit:
        page.evaluate(f"window.scrollTo(0, {int(shown)})")
        page.wait_for_timeout(150)
        want = shown - page.evaluate(COVERED_JS)
        page.evaluate(f"window.scrollTo(0, {int(want)})")
        y = page.evaluate("window.scrollY")
        repeat = shown - y
        if view - repeat < MIN_NEW_ROWS:
            break
        # Below `want` the page has reached its end: the rows above `repeat` were shot already.
        clip = {"x": 0, "y": repeat, "width": width, "height": view - repeat} if y < want else None
        snap(f"s{len(names)}", clip)
        shown = y + view
    return names, shown + MIN_NEW_ROWS <= total


def shoot(args, names, base_url, probe):
    # Imported here so that --summary on probe files runs without Playwright.
    from playwright.sync_api import sync_playwright

    out, themes, widths = args.out, probe["themes"], probe["widths"]
    full_widths, axe_widths = set(probe["full_widths"]), set(probe["axe_widths"])
    to_end = set(probe["to_end_widths"])
    run_axe = args.axe and Axe is not None
    with sync_playwright() as pw:
        # Headless Chromium hides every scrollbar unless this default is dropped.
        browser = pw.chromium.launch(ignore_default_args=["--hide-scrollbars"])
        for theme in themes:
            context = browser.new_context(device_scale_factor=args.scale)
            context.add_init_script(f"localStorage.setItem('mdbook-theme', {json.dumps(theme)})")
            page = context.new_page()
            for name in names:
                slug = name.replace("/", "-")
                entry = probe["pages"].setdefault(name, {"slug": slug, "widths": {}})
                url = base_url + name + ".html"
                for width in widths:
                    full = width in full_widths
                    if not full and theme != themes[0]:
                        continue
                    height = HEIGHTS.get(width, 900)
                    page.set_viewport_size({"width": width, "height": height})
                    page.goto(url, wait_until="load")
                    page.evaluate("document.fonts.ready")
                    # A repeat load lays the page out before the cached web font applies and can
                    # leave a block overflowing by a few pixels with no scrollbar; a resize relays it.
                    page.set_viewport_size({"width": width + 1, "height": height})
                    page.set_viewport_size({"width": width, "height": height})
                    rec = entry["widths"].setdefault(str(width), {"viewport": [width, height], "tiles": []})
                    if theme == themes[0]:
                        rec.update(page.evaluate(PROBE_JS, 1))
                    limit = (args.max_tiles if width in to_end else 3) if full else 1
                    stem = f"{slug}-{width}-{theme}"
                    tiles, more = shoot_tiles(page, width, limit, out, stem)
                    rec["tiles"].extend(tile_name(stem, tile) for tile in tiles)
                    if more and width in to_end:
                        print(f"{name} at {width} px: the page continues past {limit} tiles")
                    if run_axe and full and width in axe_widths:
                        page.evaluate("window.scrollTo(0, 0)")
                        result = Axe().run(page, options={"runOnly": {"type": "tag", "values": AXE_TAGS}})
                        rec.setdefault("axe", {})[theme] = axe_summary(result.response["violations"])
                if theme == themes[-1]:
                    recs = entry["widths"].values()
                    over = [w for w, r in entry["widths"].items() if r["overflow"]]
                    tiles = sum(len(r["tiles"]) for r in recs)
                    print(f"{name}: {len(recs)} widths, {tiles} tiles, overflow at {', '.join(over) or 'none'}")
            context.close()
        browser.close()


def summary_html(rec):
    def bad(flag, text):
        return f'<b class="bad">{text}</b>' if flag else text

    scroll = {kind: [s for s in rec["scrollers"] if s["kind"] == kind] for kind in ("pre", "table")}
    over = f"yes ({rec['overflowing_total']} elements)" if rec["overflow"] else "no"
    page_font = "" if rec["min_font_page_px"] == rec["min_font_px"] else f" (page {rec['min_font_page_px']} px)"
    parts = [
        bad(rec["overflow"], f"overflow {over}"),
        f"{rec['chars_per_line']} n/line, {rec['chars_per_line_text']} chars/line ({rec['measure_px']} px)",
        f"min font {rec['min_font_px']} px{page_font}",
    ]
    for kind, items in scroll.items():
        hidden = sum(not s["scrollbar_visible"] for s in items)
        parts.append(bad(hidden, f"{len(items)} {kind} scroll ({hidden} without scrollbar)"))
    parts.append(bad(rec["wrapped_inline_code"], f"{len(rec['wrapped_inline_code'])} inline code wrapped"))
    parts.append(f"{len(rec['long_span_wraps'])} long-span wraps")
    for theme, axe in rec.get("axe", {}).items():
        ids = ", ".join(f"{v['id']} {v['nodes']}" for v in axe["violations"])
        detail = f" ({ids})" if ids else ""
        parts.append(bad(axe["count"], f"axe {theme}: {axe['count']} violations{detail}"))
    return " | ".join(parts)


def write_sheet(probe, out):
    widths = ", ".join(map(str, probe["widths"]))
    full_widths = ", ".join(map(str, probe["full_widths"]))
    doc = [
        "<!doctype html><meta charset=utf-8><title>Documentation site contact sheet</title>",
        "<style>body{font:14px/1.4 system-ui,sans-serif;margin:1em 2em}h2{margin-top:2em}",
        ".tiles img{max-height:240px;max-width:360px;border:1px solid #999;margin:2px;vertical-align:top}",
        ".probe{font:12px/1.5 monospace;margin:.3em 0 1em}.bad{color:#b00000}nav a{margin-right:1em}</style>",
        (f"<h1>Contact sheet</h1><p>site {html.escape(probe['site'])} at {html.escape(probe['base'])}, "
         f"themes {', '.join(probe['themes'])}, widths {widths}, full widths {full_widths}, "
         f"axe {probe['axe']}</p>"),
        "<nav>" + " ".join(f'<a href="#{e["slug"]}">{html.escape(n)}</a>'
                           for n, e in probe["pages"].items()) + "</nav>",
    ]
    for name, entry in probe["pages"].items():
        first = next(iter(entry["widths"].values()))
        title = html.escape(first.get("title", ""))
        doc.append(f'<section><h2 id="{entry["slug"]}">{html.escape(name)}</h2><p>{title}</p>')
        for width, rec in entry["widths"].items():
            imgs = "".join(f'<a href="{f}"><img src="{f}" loading="lazy" alt="{f}"></a>' for f in rec["tiles"])
            doc.append(f'<h3>{width} px</h3><div class="tiles">{imgs}</div>'
                       f'<p class="probe">{summary_html(rec)}</p>')
        doc.append("</section>")
    with open(os.path.join(out, "contact-sheet.html"), "w", encoding="utf-8") as f:
        f.write("\n".join(doc) + "\n")


def summary_md(probe):
    """One Markdown table row per page and width, then the totals.

    A probe without `long_span_wraps` shows them as n/a."""
    def num(value):
        return "n/a" if value is None else f"{value:g}"

    cells = []
    for name, entry in probe["pages"].items():
        for width, rec in entry["widths"].items():
            spans = rec.get("long_span_wraps")
            cells.append({"page": name, "width": int(width), "tiles": len(rec["tiles"]),
                          "overflow": rec["overflow"], "wrapped": len(rec["wrapped_inline_code"]),
                          "long": None if spans is None else len(spans),
                          "bare": sum(not s["scrollbar_visible"] for s in rec["scrollers"]),
                          "font": rec["min_font_px"], "cpl": rec["chars_per_line_text"]})
    lines = [f"Probe of {probe['site']} at {probe.get('base', '/')}, themes {', '.join(probe['themes'])}.",
             "", ("| Page | Width | Tiles | Overflow | Wrapped code | Long-span wraps | Scrollers without a bar "
                  "| Min font px | Chars per line |"), "| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |"]
    lines += [f"| {c['page']} | {c['width']} | {c['tiles']} | {'yes' if c['overflow'] else 'no'} | {c['wrapped']} "
              f"| {num(c['long'])} | {c['bare']} | {num(c['font'])} | {num(c['cpl'])} |" for c in cells]
    def of(n):
        return f"{n} of {len(cells)} cells"

    def count(key):
        return f"{sum(c[key] or 0 for c in cells)} in {of(sum(bool(c[key]) for c in cells))}"

    split = all(c["long"] is not None for c in cells)
    fonts = sorted((c["font"], c["page"], c["width"]) for c in cells if c["font"] is not None)
    off = [c for c in cells if c["width"] >= 768 and c["cpl"] is not None and not 60 <= c["cpl"] <= 80]
    wrapped = ("wrapped inline code that fits its block" if split
               else "inline code on two or more lines, fit not recorded")
    smallest = "{:g} px ({} at {:d} px)".format(*fonts[0]) if fonts else "n/a"
    lines += ["", f"- overflow: {of(sum(c['overflow'] for c in cells))}",
              f"- {wrapped}: {count('wrapped')}",
              f"- long-span wraps: {count('long') if split else 'not recorded in this probe'}",
              f"- scrollers without a visible bar: {count('bare')}",
              f"- minimum font: {smallest}",
              f"- characters per line outside 60 to 80 at 768 px and wider: {of(len(off))}"]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--site", help="the built site directory")
    ap.add_argument("--base", help="the URL path the site is served at; default the site-url of its 404.html")
    ap.add_argument("--out", help="where the screenshots, probe.json, summary.md and the sheet go")
    ap.add_argument("--pages", default="all", help="comma-separated page paths without .html, or all")
    ap.add_argument("--widths", default=",".join(map(str, WIDTHS)))
    ap.add_argument("--full-widths", default=",".join(map(str, FULL_WIDTHS)),
                    help="widths shot in every theme with scrolled tiles")
    ap.add_argument("--to-end-widths", default=",".join(map(str, TO_END_WIDTHS)),
                    help="full widths whose tiles run to the foot of the page")
    ap.add_argument("--max-tiles", type=int, default=100, help="tiles per page at a --to-end-widths width")
    ap.add_argument("--themes", default="light,navy")
    ap.add_argument("--axe", action="store_true", help="run axe-core at the --axe-widths")
    ap.add_argument("--axe-widths", default=",".join(map(str, AXE_WIDTHS)))
    ap.add_argument("--scale", type=float, default=1.0, help="device scale factor of the screenshots")
    ap.add_argument("--summary", nargs="*", metavar="PROBE",
                    help="print summary.md after the run; given probe.json files, print theirs and stop")
    args = ap.parse_args()

    if args.summary:
        for path in args.summary:
            with open(path, encoding="utf-8") as f:
                print(f"## {path}\n\n{summary_md(json.load(f))}")
        return
    if not (args.site and args.out):
        ap.error("--site and --out are required unless --summary names probe files")
    site = os.path.abspath(args.site)
    path = (args.base or site_base(site)).strip("/")
    base = "/" + path + "/" if path else "/"
    names = find_pages(site, args.pages)
    os.makedirs(args.out, exist_ok=True)
    if args.axe and Axe is None:
        print("axe skipped")
    probe = {
        "site": site, "base": base, "themes": args.themes.split(","),
        "widths": ints(args.widths), "full_widths": ints(args.full_widths),
        "to_end_widths": ints(args.to_end_widths), "axe_widths": ints(args.axe_widths),
        "axe": "run" if args.axe and Axe is not None else "skipped", "pages": {},
    }
    server, base_url = serve(site, base)
    try:
        shoot(args, names, base_url, probe)
    finally:
        server.shutdown()
    with open(os.path.join(args.out, "probe.json"), "w", encoding="utf-8") as f:
        json.dump(probe, f, indent=1)
    write_sheet(probe, args.out)
    summary = summary_md(probe)
    with open(os.path.join(args.out, "summary.md"), "w", encoding="utf-8") as f:
        f.write(summary)
    if args.summary is not None:
        print(summary)
    tiles = sum(len(r["tiles"]) for e in probe["pages"].values() for r in e["widths"].values())
    over = sum(r["overflow"] for e in probe["pages"].values() for r in e["widths"].values())
    print(f"{tiles} screenshots of {len(names)} pages at {len(probe['widths'])} widths, "
          f"{over} page-width cells overflow; {os.path.join(args.out, 'probe.json')} and "
          f"{os.path.join(args.out, 'contact-sheet.html')}")


if __name__ == "__main__":
    main()
