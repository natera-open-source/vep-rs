"""Render the manuscript figures in every output format.

Usage::

    python3 manuscript/figures/generate_figures.py
    python3 manuscript/figures/generate_figures.py --only fig1
    python3 manuscript/figures/generate_figures.py --only figS1

`fig1_concordance_throughput` is the main-text figure: (a) concordance per dataset for
both Rust engines, (b) three-engine wall time on both architectures, (c) vep-rs speedup
over VEP on both architectures. `figS1_per_class_concordance` is supplementary Figure
S1, one panel carrying every Sequence Ontology term of Table S3.

Inputs are read from ../data/: f1_observations.csv, f1_observations_fastvep.csv and
wall_times.csv for Figure 1, f1_by_consequence_class.csv for Figure S1. Output lands
beside this script as .png, .pdf, .eps and .tif, plus a MANIFEST.txt listing the sha256
of every .png, .pdf and .eps written. Delivered widths and Figure S1's row pitch are
measured after layout (`_assert_delivered_width`, `_assert_row_pitch`), because
`savefig.bbox` is tight and the crop depends on the installed fonts. The make_*.py
scripts beside this one build the other released CSVs; no figure reads them.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import numpy as np
import pandas as pd

try:
    import scienceplots  # noqa: F401  -- registers styles via import side effect

    _HAS_SCIENCEPLOTS = True
except ImportError:
    _HAS_SCIENCEPLOTS = False

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FIG_DIR = Path(__file__).resolve().parent


# wall_times.csv can hold `median_of_20` rows from more than one sweep on the same host
# class, so every selector below is a sweep id and no figure filters on a host or
# instance type. An empty selection would still draw a figure, so every selection asserts
# its row count and an empty or partial one is a hard error.
VEP_RS_SWEEP_ID = "sweep-20260920T012411Z"

# fastVEP and Perl carry their own pins because vep-rs and its comparators are measured
# in separate sweeps; a shared constant would select zero rows for one side.
FASTVEP_WALLTIME_SWEEP_ID = "sweep-20260918T224520Z"
PERL_WALLTIME_SWEEP_ID = "sweep-20260918T224520Z"

# vep-rs SV wall-time pin (the vep-rs column of supplementary Table S7): the sweep the
# structural-variant wall-time rows carry, kept as its own constant beside
# VEP_RS_SWEEP_ID. Unused by either figure: panel (b) plots the six SNP/indel datasets
# and no SV cell.
VEP_RS_SV_SWEEP_ID = "sweep-20260920T012411Z"

# fastVEP concordance pin (the fastVEP rows of Table S2 and Figure 1 panel (a)). Pinned by
# sweep id, never by a notes keyword: a keyword can match rows from more than one sweep,
# while a sweep id names exactly the rows that sweep measured.
FASTVEP_F1_SWEEP_ID = "sweep-20260918T224520Z"

# A comparator's two structural-variant cells sit on their own sweep, separate from its
# SNP and indel cells; each carries its own pin so a selector names the sweep the CSV
# rows carry. Panel (a) reads FASTVEP_SV_F1_SWEEP_ID for fastVEP's two SV rows; the two
# wall-time pins name the sweep of the comparators' Table S7 cells and feed no figure.
FASTVEP_SV_WALLTIME_SWEEP_ID = "sweep-20260920T012411Z"
PERL_SV_WALLTIME_SWEEP_ID = "sweep-20260920T012411Z"
FASTVEP_SV_F1_SWEEP_ID = "sweep-20260920T012411Z"


def hu(value: float, places: int = 2) -> str:
    """Render a float at `places` decimals under ROUND_HALF_UP.

    The manuscript tables round half away from zero; f-string formatting rounds half to
    even, so a value whose next digit is exactly 5 would render one unit lower here than
    in the tables.
    """
    return str(Decimal(str(value)).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))


def sig3(value: float) -> str:
    """Render a ratio of medians at three significant figures under ROUND_HALF_UP, the
    precision the manuscript prints every speedup at (S1.2): 176.29 -> 176, 78.26 -> 78.3,
    64.955 -> 65.0. A carry that adds a fourth figure (99.95 -> 100.0) is requantized."""
    d = Decimal(str(value))
    exponent = d.adjusted() - 2
    r = d.quantize(Decimal(1).scaleb(exponent), rounding=ROUND_HALF_UP)
    if r.adjusted() > d.adjusted():
        exponent += 1
        r = r.quantize(Decimal(1).scaleb(exponent), rounding=ROUND_HALF_UP)
    if exponent > 0:
        r = r.quantize(Decimal(1))
    return f"{r:f}"


def _sweep_pin(wt, sweep_id: str):
    """Boolean mask selecting rows tagged with `sweep_id` in the notes field."""
    return wt["notes"].astype(str).str.contains(f"sweep_id={sweep_id}", na=False)


def _fastvep_f1_pin(fastvep_f1, suites):
    """Select the fastVEP concordance rows `FASTVEP_F1_SWEEP_ID` pins, or exit.

    Returns those rows for `suites`. Exits when the count is not exactly `len(suites)`:
    too few means the pin names rows the CSV does not carry (an empty selection would
    still draw a figure), too many means two sweeps matched and the panel would plot
    duplicates.
    """
    sel = fastvep_f1[
        _sweep_pin(fastvep_f1, FASTVEP_F1_SWEEP_ID) & (fastvep_f1["suite"].isin(suites))
    ].copy()
    if len(sel) != len(suites):
        present = sorted(
            fastvep_f1["notes"]
            .astype(str)
            .str.extract(r"sweep_id=([^;]+)")[0]
            .dropna()
            .unique()
        )
        raise SystemExit(
            f"fastVEP F1 pin sweep_id={FASTVEP_F1_SWEEP_ID} selected {len(sel)} row(s) "
            f"for {len(suites)} suites in f1_observations_fastvep.csv "
            f"(sweeps present: {present}). "
            "Point FASTVEP_F1_SWEEP_ID at a sweep the CSV carries."
        )
    return sel


def _assert_one_row_per_cell(df, engine_suite_cols, label: str) -> None:
    """Exit when a filter leaves more than one row per `engine_suite_cols` cell.

    Each figure plots one measurement per cell; a groupby over duplicate rows would
    average distinct sweeps or builds. `label` names the panel in the message.
    """
    dupes = df.groupby(engine_suite_cols).size()
    dupes = dupes[dupes > 1]
    if not dupes.empty:
        raise SystemExit(
            f"{label}: filter left >1 row for {len(dupes)} (engine, suite) cell(s): "
            f"{dupes.to_dict()}. Taking a median across them would mix sweeps/builds. "
            "Pin the sweep id (and condition=/arch= tags) so exactly one row survives."
        )


# Okabe-Ito 8-color colorblind-safe palette (Wong 2011, Nature Methods).
OK = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
}

# One hue per engine and nothing else encoded in hue: architecture and eligibility are
# carried by marker fill. A metric encoded as the engine hue tinted toward white puts two
# different engines within OKLab dE 7 of each other under normal vision, closer under
# protanopia, so the group cannot be read by engine. Blue, orange and green stay within
# one lightness band, above a chroma floor and apart under the common colour-vision
# deficiencies (vermillion sits too close to orange), and blue against orange also
# separates in grayscale, L* about 45 against 70, for a black-and-white print.
ENGINES = {
    "vep_rs": OK["blue"],
    "fastvep": OK["orange"],
    "perl": OK["green"],
}

SUITE_LABELS = {
    "clinvar_full_grch37": "ClinVar GRCh37",
    "clinvar_full_grch38": "ClinVar GRCh38",
    "gnomad_v2.1.1_chr21": "gnomAD v2.1.1 chr21",
    "gnomad_v4.1_chr21": "gnomAD v4.1 chr21",
    "1kg_integrated_chr21": "1000G Phase 3 chr21",
    "1kg_highcov_chr21": "1000G high-cov chr21",
    # wall_times.csv names the two ClinVar suites without the "_full" infix that
    # f1_observations_fastvep.csv uses; both spellings are aliased rather than a column
    # renamed, which would break the F1 lookups.
    "clinvar_grch37": "ClinVar GRCh37",
    "clinvar_grch38": "ClinVar GRCh38",
}


def setup_style() -> None:
    """Apply publication-quality rcParams once at startup."""
    if _HAS_SCIENCEPLOTS:
        plt.style.use(["science", "no-latex"])
    plt.rcParams.update(
        {
            "font.family": "serif",
            # STIX Two Text first, a Times-family face; DejaVu Serif is a wide screen
            # face that reads as mismatched beside a Times-set body, and a bare
            # "STIX" names a family matplotlib does not ship.
            "font.serif": ["STIX Two Text", "Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "figure.dpi": 100,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05,
            "pdf.fonttype": 42,  # Type 42 (TrueType), so PDF text stays searchable
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            # The `science` style sets xtick.top/ytick.right True and both
            # minor.visible True, which draws ticks on the two edges the spine
            # settings above remove, so these four keys are turned off after
            # plt.style.use.
            "xtick.top": False,
            "ytick.right": False,
            "xtick.minor.visible": False,
            "ytick.minor.visible": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "axes.linewidth": 0.6,
            "axes.grid": True,
            # A solid pale grid, not an alpha'd dark one: the PostScript backend
            # renders partial transparency opaque, so a `grid.alpha` below 1 gives
            # the .eps a heavier grid than the .png. This grey is black at 0.3 on
            # white, so all four output formats agree.
            "grid.alpha": 1.0,
            "grid.color": "#D1D5DB",
            "grid.linewidth": 0.4,
            "legend.frameon": False,
        }
    )


def _save(fig, name: str) -> list[Path]:
    """Save `fig` as PNG, PDF, EPS and TIFF under `name` and return the paths.

    PNG takes the global `savefig.dpi`, because its digest is the byte contract
    MANIFEST.txt records. TIFF is the raster fallback: LZW is lossless and shrinks
    line art by orders of magnitude at its dpi.
    """
    per_format_kwargs = {
        "tif": {"dpi": 1200, "pil_kwargs": {"compression": "tiff_lzw"}}
    }
    paths = []
    for ext in ("png", "pdf", "eps", "tif"):
        path = FIG_DIR / f"{name}.{ext}"
        fig.savefig(path, **per_format_kwargs.get(ext, {}))
        paths.append(path)
    plt.close(fig)
    return paths


# 178 mm is a double-column width. (a) spans the top row and (b) and (c) share the row
# beneath it, so (a) gets the full width its bars' decades of log axis need. Type is 7 pt
# minimum at final size.
FIG1_WIDTH_IN = 178 / 25.4
FIG1_HEIGHT_IN = 130 / 25.4

# Panel (a) and Figure S1 both plot F1 on an axis logarithmic in 1 - F1 and reversed, so
# further right is higher F1: on a linear 0-1 axis every vep-rs bar is full height and
# identical, while in log(1 - F1) the two engines separate by orders of magnitude. A logit
# axis is infinite at both ends, and terms sit at F1 exactly 0 and exactly 1; log(1 - F1)
# needs a clip at one end only, which is this floor.
DISCORD_FLOOR = 1e-7

# The worse end of both axes: 1 - F1 = 1.0 is F1 = 0. Panel (a)'s bars are anchored here
# and grow toward better, so a longer bar is a better one.
CONCORDANCE_ANCHOR = 1.0


# Shared by both figures so the two cannot label the same transform differently. Positions are
# 1 - F1; labels are the F1 they correspond to. `1 - 1e-7` formatted naively is "1.0", which
# would print the perfect-agreement tick as if it were exact agreement, so the label is built
# from the exponent instead: 1e-k -> "0." + k nines.
def _f1_tick_labels(positions: list[float]) -> list[str]:
    labels = []
    for pos in positions:
        if pos >= 1.0:
            labels.append("0")
        else:
            nines = int(round(-math.log10(pos)))
            labels.append("0." + "9" * nines)
    return labels


def _fig1_panel_a(ax, df, fv_df, sv_rows) -> None:
    """(a) Concordance per dataset, both Rust engines, horizontal grouped bars."""
    labels, vr, fv = [], [], []
    for _, r in df.iterrows():
        labels.append(SUITE_LABELS.get(r["suite"], r["suite"]))
        vr.append(max(1.0 - float(r["raw_f1"]), DISCORD_FLOOR))
    for _, r in fv_df.iterrows():
        fv.append(max(1.0 - float(r["raw_f1"]), DISCORD_FLOOR))
    for lbl, v, f in sv_rows:
        labels.append(lbl)
        vr.append(max(1.0 - v, DISCORD_FLOOR))
        fv.append(max(1.0 - f, DISCORD_FLOOR))

    y = np.arange(len(labels))
    h = 0.38
    # Each bar spans [datum, CONCORDANCE_ANCHOR] in 1 - F1 and renders from the left spine
    # rightward once the axis is reversed. A negative width with `left` at the anchor keeps
    # both endpoints explicit, which a log scale requires.
    ax.barh(y - h / 2, [v - CONCORDANCE_ANCHOR for v in vr], height=h,
            left=CONCORDANCE_ANCHOR, color=ENGINES["vep_rs"], label="vep-rs", linewidth=0)
    ax.barh(y + h / 2, [v - CONCORDANCE_ANCHOR for v in fv], height=h,
            left=CONCORDANCE_ANCHOR, color=ENGINES["fastvep"], label="fastVEP", linewidth=0)
    ax.set_xscale("log")
    # Descending limits reverse the axis: the better end is the right-hand one.
    ax.set_xlim(CONCORDANCE_ANCHOR, DISCORD_FLOOR)
    # One tick per decade, as Figure S1 draws, so a bar ending near 0.9 or 0.999 reads to
    # within its decade; the 1e-7 limit itself is unlabelled because a "0.9999999" label
    # would be centred on the right spine and hang past it, widening the tight bbox.
    _ticks = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    ax.set_xticks(_ticks)
    ax.set_xticklabels(_f1_tick_labels(_ticks), fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    # A semicolon joins the two clauses and the minus is U+2212, matching the captions'
    # "1 − F1".
    ax.set_xlabel("Concordance (F1), log(1 − F1) scale; further right is better", fontsize=8)
    ax.grid(axis="y", visible=False)
    # No per-bar value labels: the caption points at the tables, so the figure states no
    # number of its own. The legend sits above the axes; inside, it collides with bar ends
    # in every corner, since the longest bars reach the right spine.
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=7,
              handlelength=1.2, columnspacing=1.0)


# The eligibility threshold the reported floor is defined over, in pooled VEP tuples.
# Mirrors the same-named constant in scripts/concordance/compare_vep_outputs.py; a term
# below it is drawn open in Figure S1 and daggered in Supplementary Table S3, so the
# renderer, the comparator and the table caption all have to agree on it.
PER_CLASS_MIN_TUPLES = 1000

# The declared figsize is the input to the tight-bbox crop, not the delivered size, so both
# numbers are held by `_assert_delivered_width` and `_assert_row_pitch` rather than reasoned
# to. The legend above and the term labels to the left sit inside the crop, and their
# metrics depend on the installed fonts (a host without STIX Two Text crops several mm
# wider), so a margin measured on one host is not a margin on another. What overshoots is
# fixed-width text, so narrowing the declared width narrows the delivered width about one
# for one.
FIGS1_WIDTH_IN = 142 / 25.4
FIGS1_HEIGHT_IN = 200 / 25.4

# 1.5x the 7 pt minimum label size, in mm. Below this, 7 pt tick labels on adjacent rows
# touch.
FIGS1_MIN_ROW_PITCH_MM = 1.5 * 7 / 72 * 25.4


# The supplement's text width: letter paper (215.9 mm) at 1 in (25.4 mm) margins leaves
# 165.1 mm. A figure delivered wider than this is scaled down to fit, which takes the 7 pt
# labels below 7 pt.
SUPP_TEXT_WIDTH_MM = 165.0


def _assert_legend_within_axes(fig, ax, name: str) -> None:
    """Exit when `ax`'s legend extends past either side of `ax`.

    A legend drawn above its axes and wider than them overhangs the neighbouring panel or
    the page edge; its width is set by the font, which differs between hosts, so the fit
    is measured on each render.
    """
    fig.canvas.draw()
    _r = fig.canvas.get_renderer()
    _leg = ax.get_legend().get_window_extent(_r)
    _axes = ax.get_window_extent(_r)
    if _leg.x0 < _axes.x0 or _leg.x1 > _axes.x1:
        raise SystemExit(
            f"ERROR: [generate_figures] {name} legend spans {(_leg.width / fig.dpi) * 25.4:.1f} mm "
            f"against axes of {(_axes.width / fig.dpi) * 25.4:.1f} mm, so it overhangs a side. "
            "Shorten its entries or use fewer columns."
        )


def _assert_delivered_width(fig, name: str, limit_mm: float) -> None:
    """Exit when `fig`'s tight bounding box is wider than `limit_mm`.

    `savefig.bbox` is "tight", so the delivered width is the bounding box of every drawn
    artist, legend and tick labels included, not the declared `figsize`; a page layout
    absorbs an overshoot by downscaling the figure and its type. `pdfinfo` on the saved
    file measures the same quantity, but only after the artifact is written.
    """
    _w_mm = fig.get_tightbbox().width * 25.4
    if _w_mm > limit_mm:
        raise SystemExit(
            f"ERROR: [generate_figures] {name} is delivered {_w_mm:.1f} mm wide, past the "
            f"{limit_mm:.0f} mm text width, so it would be downscaled and its 7 pt type would "
            "fall below 7 pt. Narrow the declared figsize; the tight bbox exceeds it."
        )
    # Printed on a pass as well: the overshoot is font-dependent, so the margin is read off
    # each host's own render rather than inferred from the assertion not firing.
    print(f"{name}: delivered width {_w_mm:.2f} mm (limit {limit_mm:.1f} mm)")


def _assert_row_pitch(fig, ax, n_rows: int) -> None:
    """Exit when adjacent rows of `ax` are closer than `FIGS1_MIN_ROW_PITCH_MM`.

    The pitch is measured after layout through the axis's own data transform, because the
    height an axis receives depends on what the tight bbox crops, on the legend above it
    and on the x label below it, none of which `figsize` predicts. Rows sit one data unit
    apart and the axis spans more than `n_rows` units, since matplotlib pads each end of
    the y axis by 5% of the data range, so axis height over row count overstates the
    pitch. `n_rows` also asserts that every row lies inside the drawn y limits.
    """
    fig.canvas.draw()
    # x must be a valid coordinate on a log axis, so 1.0; only the y pitch is read.
    (_, _y0), (_, _y1) = ax.transData.transform([(1.0, 0.0), (1.0, 1.0)])
    _pitch = abs(_y1 - _y0) / fig.dpi * 25.4
    _lo, _hi = sorted(ax.get_ylim())
    if not (_lo <= 0 and _hi >= n_rows - 1):
        raise SystemExit(
            f"ERROR: [generate_figures] figS1 draws rows 0..{n_rows - 1} but its y limits are "
            f"{ax.get_ylim()}, so some rows lie outside the axis"
        )
    if _pitch < FIGS1_MIN_ROW_PITCH_MM:
        raise SystemExit(
            f"ERROR: [generate_figures] figS1 gives each of {n_rows} rows "
            f"{_pitch:.2f} mm of axis, under the {FIGS1_MIN_ROW_PITCH_MM:.2f} mm that 7 pt "
            "labels need to clear each other. Raise FIGS1_HEIGHT_IN."
        )
    print(f"figS1: row pitch {_pitch:.2f} mm over {n_rows} rows "
          f"(minimum {FIGS1_MIN_ROW_PITCH_MM:.2f} mm)")


def figS1_per_class_concordance() -> list[Path]:
    """Supplementary Figure S1: per-Sequence-Ontology-term concordance, every term.

    Every row of Table S3 is drawn at 7 pt, ranked worst first on vep-rs rather than by
    tuple count: the floor term carries orders of magnitude fewer pooled tuples than
    `intron_variant`, so a volume order buries the value the caption names, and Table S3
    is already the volume order.

    A rule separates the terms VEP emits, one mark per engine, from the terms only fastVEP
    emits, which have no vep-rs mark and a fastVEP F1 of zero by construction. A term below
    `PER_CLASS_MIN_TUPLES` is drawn open and daggered, both, because open against filled
    survives a greyscale print where a hue would not. A right-pointing mark at the right
    limit is exact agreement, both one-sided counts zero: 1 - F1 = 0 has no place on a log
    axis, and clipped to the floor it would read as 0.9999999, so the mark is chosen on the
    integer columns rather than on the six-decimal `f1` (which prints 1.000000 for a term
    with a few one-sided tuples in ten million) and points toward the better end of the
    reversed axis; such a near-exact term is drawn at its true 1 - F1, just left of the
    limit.

    A vep-rs row whose `adj_f1` differs from its `f1` carries a second, smaller diamond at
    the adjusted value (the S2 mask's exclusions removed from both sides), joined to the raw
    mark by a hairline so the pair reads as one term's two views; exact agreement after the
    mask takes the same right-pointing shape, smaller. Rows where the two agree draw
    one mark, so the figure adds ink only where the mask changes the reading.
    """
    per_class_path = DATA_DIR / "f1_by_consequence_class.csv"
    if not per_class_path.is_file():
        raise SystemExit(
            f"ERROR: [generate_figures] {per_class_path} is missing; figS1 has no data. "
            "This CSV is committed, so its absence is a broken tree rather than a pending "
            "measurement."
        )

    raw = pd.read_csv(per_class_path)
    # 1 - F1 from the row's integers, not from the six-decimal `f1` string: a term with a
    # handful of one-sided tuples in ten million prints as 1.000000 yet sits at 1 - F1 of a
    # few 1e-7, which this axis can show (its floor is DISCORD_FLOOR), and a term whose two
    # tuple sets are identical has 1 - F1 = 0 exactly, which it cannot. Only the latter is
    # drawn at the limit with the right-pointing mark.
    _one_sided = raw["vep_only"].astype(int) + raw["engine_only"].astype(int)
    _both = raw["vep_tuples"].astype(int) + raw["engine_tuples"].astype(int)
    raw["exact"] = _one_sided == 0
    raw["discordance"] = (_one_sided / _both.where(_both > 0, 1)).clip(lower=DISCORD_FLOOR)
    if "adj_f1" not in raw.columns:
        raise SystemExit(
            f"ERROR: [generate_figures] {per_class_path} carries no adj_f1 column; the "
            "pooled CSV is written by compare_vep_outputs.py --pool-consequence-class "
            "--per-term, and Figure S1 draws the adjusted series from it."
        )
    per_engine = {e: g.set_index("term") for e, g in raw.groupby("engine")}
    vep_rs, fastvep = per_engine["vep-rs"], per_engine.get("fastvep")

    # fastVEP shares vep-rs's row order, so the two are read against each other.
    ranked = vep_rs.sort_values("discordance", ascending=False)
    # Terms only fastVEP emits, kept in Table S3's own order so the figure and the table agree.
    fv_only = (
        [t for t in fastvep.index if t not in vep_rs.index] if fastvep is not None else []
    )
    labels = list(ranked.index) + fv_only

    fig, ax = plt.subplots(figsize=(FIGS1_WIDTH_IN, FIGS1_HEIGHT_IN))

    for engine, key in (("vep_rs", "vep-rs"), ("fastvep", "fastvep")):
        table = per_engine.get(key)
        if table is None:
            continue
        for row_i, term in enumerate(labels):
            if term not in table.index:
                continue
            _row = table.loc[term]
            # Shape encodes exact agreement and fill encodes eligibility, independently: a
            # term can be both, and a solid exact-agreement marker would falsify the caption's
            # "an open mark is a term below the threshold". Eligibility is a property of the
            # term, not the engine: `vep_tuples` is the VEP-side total, so both engines' marks
            # on a shared row take the same fill; it is read per row because a fastVEP-only
            # term has no vep-rs row, and there `vep_tuples` is 0.
            _exact = bool(_row["exact"])
            _eligible = int(_row["vep_tuples"]) >= PER_CLASS_MIN_TUPLES
            ax.scatter(
                [float(_row["discordance"])], [row_i], s=15,
                marker=">" if _exact else "o",
                facecolors=ENGINES[engine] if _eligible else "white",
                edgecolors=ENGINES[engine], linewidths=0.7, zorder=3,
            )
            # The adjusted view, vep-rs only, where the mask touches the term: a smaller
            # diamond at the adjusted 1 - F1, computed from the integers with the excluded
            # pairs removed from both sides, joined to the raw mark by a hairline. The
            # exclusion columns are empty for fastVEP (pandas reads them as NaN) and zero
            # where no excluded pair touches the term, and neither case draws anything.
            if key == "vep-rs" and pd.notna(_row.get("vep_only_excluded")):
                _ex_v, _ex_e = int(_row["vep_only_excluded"]), int(_row["engine_only_excluded"])
                if _ex_v or _ex_e:
                    _adj_one_sided = (int(_row["vep_only"]) - _ex_v) + (int(_row["engine_only"]) - _ex_e)
                    _adj_both = int(_row["vep_tuples"]) + int(_row["engine_tuples"]) - _ex_v - _ex_e
                    _adj_exact = _adj_one_sided == 0
                    _adj_disc = max(_adj_one_sided / _adj_both, DISCORD_FLOOR)
                    ax.plot(
                        [float(_row["discordance"]), _adj_disc], [row_i, row_i],
                        color=ENGINES[engine], linewidth=0.6, alpha=0.6, zorder=2,
                    )
                    ax.scatter(
                        [_adj_disc], [row_i], s=9,
                        marker=">" if _adj_exact else "D",
                        facecolors=ENGINES[engine], edgecolors="white",
                        linewidths=0.4, zorder=4,
                    )

    # The rule sits between rows and consumes no row slot, so the pitch assertion measures
    # what the labels get.
    if fv_only:
        ax.axhline(len(ranked) - 0.5, color="#9CA3AF", linewidth=0.6, zorder=1)
        # Right-aligned: the fastVEP-only marks sit at 1 - F1 = 1.0, the left spine of the
        # reversed axis, so the empty band below the rule is on the right. Centred on the
        # first fastVEP-only row rather than just under the rule: 7 pt type is about 2.4 mm
        # tall, and a label a fraction of a row under the rule is struck through by it. The
        # opaque white box masks the gridlines the text crosses.
        ax.annotate(
            "terms only fastVEP emits (no VEP tuples)",
            xy=(0.985, len(ranked)), xycoords=("axes fraction", "data"),
            ha="right", va="center", fontsize=7, color="#6B7280",
            bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none"),
        )

    ax.set_xscale("log")
    # Descending limits reverse the axis. The right limit sits beyond the floor: the
    # exact-agreement marks draw at DISCORD_FLOOR, and with the limit there the triangle is
    # centred on the spine and renders as a half glyph.
    ax.set_xlim(1.6, DISCORD_FLOOR * 0.8)
    ax.set_yticks(np.arange(len(labels)))
    # The dagger is decided by the term's `vep_tuples`, never by parsing its label. A
    # fastVEP-only term reads its own row, since it has no vep-rs row.
    _ticks = []
    for t in labels:
        _src = vep_rs if t in vep_rs.index else fastvep
        _n = int(_src.loc[t, "vep_tuples"])
        _ticks.append(t + ("" if _n >= PER_CLASS_MIN_TUPLES else " †"))
    ax.set_yticklabels(_ticks, fontsize=7)
    ax.invert_yaxis()
    # Ticks stop at 1e-6: a "0.9999999" label at the 1e-7 limit is centred on the right
    # spine, so half of it hangs outside the axes and widens the tight bbox. The
    # exact-agreement marks past 1e-6 are explained by the legend rather than read off a
    # tick.
    _xticks = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    ax.set_xticks(_xticks)
    ax.set_xticklabels(_f1_tick_labels(_xticks), fontsize=7)
    # Same wording as panel (a): semicolon between the clauses, U+2212 minus.
    ax.set_xlabel("Concordance (F1), log(1 − F1) scale; further right is better", fontsize=8)
    ax.grid(axis="y", visible=False)

    # Two engines, the eligibility fill and the exact-agreement mark, so the key is complete
    # without the caption. matplotlib fills a legend column-major, so at ncol=2 the engines
    # share one column and the two mark conventions the other.
    _key = [
        Line2D([], [], marker="o", linestyle="none", markersize=3.4,
               color=ENGINES["vep_rs"], markerfacecolor=ENGINES["vep_rs"],
               markeredgecolor=ENGINES["vep_rs"], label="vep-rs"),
        Line2D([], [], marker="o", linestyle="none", markersize=3.4,
               color=ENGINES["fastvep"], markerfacecolor=ENGINES["fastvep"],
               markeredgecolor=ENGINES["fastvep"], label="fastVEP"),
        Line2D([], [], marker="o", linestyle="none", markersize=3.4, color="0.35",
               markerfacecolor="white", markeredgecolor="0.35",
               label=f"< {PER_CLASS_MIN_TUPLES:,} VEP tuples (†)"),
        Line2D([], [], marker=">", linestyle="none", markersize=3.4, color="0.35",
               markerfacecolor="white", markeredgecolor="0.35",
               label="identical tuple sets (F1 = 1 exactly)"),
        Line2D([], [], marker="D", linestyle="none", markersize=2.6,
               color=ENGINES["vep_rs"], markerfacecolor=ENGINES["vep_rs"],
               markeredgecolor="white", label="vep-rs adjusted F1 (S2 mask), where it differs"),
    ]
    # Above the axes and far narrower than them, so the legend costs height in the top margin
    # and no delivered width.
    ax.legend(
        _key, [h.get_label() for h in _key],
        loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=2, fontsize=7,
        columnspacing=1.0, handlelength=1.2, handletextpad=0.4,
    )

    _assert_row_pitch(fig, ax, len(labels))
    _assert_delivered_width(fig, "figS1", SUPP_TEXT_WIDTH_MM)
    return _save(fig, "figS1_per_class_concordance")


def _fig1_panel_b(ax, wall) -> None:
    """(b) Wall time per dataset per engine, both architectures, log y.

    Architecture is encoded by marker fill, filled ARM and open x86, which keeps the hue
    budget at one per engine. No connecting line: the x axis is six named datasets in
    source-family order, not monotone in any quantity, so a polyline over it draws a
    profile that is an artefact of the category order; sorting by tuple count to rescue
    the line would split the ClinVar GRCh37/GRCh38 pair the text argues on. The two
    architectures are dodged horizontally instead, because within one engine they differ
    by a few per cent, which on a log axis overlaps the filled and open marks.
    """
    _DODGE = 0.11
    for engine, key in (("perl", "perl"), ("fastvep", "fastvep"), ("vep_rs", "vep-rs")):
        for arch, filled in (("arm64", True), ("x86_64", False)):
            sub = wall[(wall["engine"] == key) & (wall["arch"] == arch)]
            if sub.empty:
                continue
            sub = sub.sort_values("order")
            ax.plot(
                sub["order"] + (-_DODGE if filled else _DODGE),
                sub["wall_time_sec"].astype(float),
                marker="o", markersize=3.2, linestyle="none",
                color=ENGINES[engine],
                markerfacecolor=ENGINES[engine] if filled else "white",
                markeredgecolor=ENGINES[engine], markeredgewidth=0.7,
                label={"perl": "VEP", "fastvep": "fastVEP", "vep_rs": "vep-rs"}[engine]
                if filled else None,
            )
    ax.set_yscale("log")
    ax.set_ylabel("Wall time (s), log scale", fontsize=8)
    ax.set_xticks(range(len(FIG1_DATASET_ORDER)))
    # Rotated labels extend below and left of the axis, which the tight bbox absorbs in
    # height, not width: the figure's leftmost extent is panel (a)'s dataset names.
    ax.set_xticklabels(
        [SUITE_LABELS.get(s, s) for s in FIG1_DATASET_ORDER],
        rotation=35, ha="right", fontsize=7,
    )
    ax.grid(axis="x", visible=False)
    # Proxy handles for the architecture key: fill is orthogonal to the three hues, so no
    # plotted artist means "ARM" independent of an engine. Neutral grey keeps them from
    # reading as a fourth engine.
    _arch_key = [
        Line2D([], [], marker="o", linestyle="none", markersize=3.2, color="0.35",
               markerfacecolor="0.35", markeredgecolor="0.35", label="ARM Graviton4"),
        Line2D([], [], marker="o", linestyle="none", markersize=3.2, color="0.35",
               markerfacecolor="white", markeredgecolor="0.35", label="x86 Intel"),
    ]
    _engine_handles, _engine_labels = ax.get_legend_handles_labels()
    # Above the axes: inside, the legend lands on the VEP marks at the gnomAD cells, the
    # panel's tallest points. matplotlib fills a legend column-major, so at ncol=2 the five
    # entries split 3 and 2, engines in one column and architectures in the other.
    ax.legend(
        _engine_handles + _arch_key,
        _engine_labels + [h.get_label() for h in _arch_key],
        loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=7,
        columnspacing=0.9, handlelength=1.2, handletextpad=0.4,
    )


def _fig1_panel_c(ax, ratios, arm_geomean: float, x86_geomean: float) -> None:
    """(c) Speedup of vep-rs over VEP per dataset, both architectures, log x.

    `ratios` holds the per-cell values from the pinned rows of wall_times.csv, never a
    literal here. The x ticks are the rungs of a doubling ladder (25×, 50×, 100×, 200×, ...)
    that the data reach, plus the rung below the lowest speedup, labelled "100×" and so on:
    the axis spans less than a decade, where the log locator places no major tick at all
    and labels the minor ticks in scientific notation, and evenly stepped multiples would
    read as a linear axis. No parity line at 1×: it would stretch the axis to span 1× to
    the maximum ratio and compress every mark into the rightmost sixth of the panel, and no
    dataset is near parity. The two geometric-mean lines, dashed ARM and dotted x86, stay
    because they lie inside the data range; they are named in the legend rather than beside
    the lines, because a label beside a line in the right half of the axis runs past the
    right spine at a wider font (the delivered width is font-dependent, and a build host's
    font need not be the author's), and the legend is one column so its width stays inside the axes.
    """
    labels = [SUITE_LABELS.get(s, s) for s in FIG1_DATASET_ORDER]
    y = np.arange(len(labels))
    for arch, filled, dy in (("arm64", True, -0.16), ("x86_64", False, 0.16)):
        vals = [ratios[(s, arch)] for s in FIG1_DATASET_ORDER]
        ax.scatter(
            vals, y + dy, s=14,
            facecolors=ENGINES["vep_rs"] if filled else "white",
            edgecolors=ENGINES["vep_rs"], linewidths=0.7, zorder=3,
            label="ARM Graviton4" if filled else "x86 Intel",
        )
    # Named for the architecture each summarises, since marks from both architectures sit
    # around both lines. "×" rather than "x", as the manuscript writes every ratio.
    ax.axvline(arm_geomean, color="#6B7280", linewidth=0.7, linestyle="--", zorder=1,
               label=f"ARM geo. mean {sig3(arm_geomean)}×")
    ax.axvline(x86_geomean, color="#6B7280", linewidth=0.7, linestyle=":", zorder=1,
               label=f"x86 geo. mean {sig3(x86_geomean)}×")
    ax.set_xscale("log")
    # Reading the limits runs the autoscale over the marks drawn above, so the tick set
    # follows the data; the docstring says why the log locator's own ticks are unusable.
    _lo, _hi = ax.get_xlim()
    # A doubling ladder (25x, 50x, 100x, 200x, ...) reads as a log axis, where evenly
    # stepped multiples (100x, 200x, 300x) read as a linear one; the ladder is cut to the
    # rungs the data reach, plus the rung below the lowest speedup so it is read against a
    # labelled line.
    _ladder = [25.0 * 2 ** k for k in range(0, 8)]
    _xticks = [t for t in _ladder if t <= _hi]
    _below = [t for t in _xticks if t <= _lo]
    _xticks = ([_below[-1]] if _below else []) + [t for t in _xticks if t > _lo]
    # Pulling the left limit down to the first tick draws that tick on the spine, so the
    # lowest speedup is read against a labelled line.
    ax.set_xlim(left=min(_lo, _xticks[0]))
    ax.set_xticks(_xticks)
    ax.set_xticklabels([f"{t:g}×" for t in _xticks], fontsize=7)
    # The log scale labels its minor ticks in scientific notation beside the majors; only
    # the labels go, the marks stay.
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("vep-rs speedup over VEP, log scale", fontsize=8)
    ax.grid(axis="y", visible=False)
    # Above the axes: inside, the legend lands on the low-speedup marks in the lower right.
    # One column: two columns of the geometric-mean entries are wider than the axes.
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1, fontsize=7,
              handlelength=1.6)


FIG1_DATASET_ORDER = [
    "clinvar_full_grch37",
    "clinvar_full_grch38",
    "gnomad_v2.1.1_chr21",
    "gnomad_v4.1_chr21",
    "1kg_integrated_chr21",
    "1kg_highcov_chr21",
]



def fig1_concordance_throughput() -> list[Path]:
    """Figure 1: concordance and throughput, three panels.

    Data sources, one per panel:
      (a) f1_observations.csv + f1_observations_fastvep.csv, raw F1 only
      (b) wall_times.csv
      (c) wall_times.csv, VEP-over-vep-rs ratios derived per dataset per architecture

    Every filter below pins its sweep id or binary digest and asserts its row count: an
    empty selection would still draw a figure rather than failing.
    """
    f1 = pd.read_csv(DATA_DIR / "f1_observations.csv")
    fastvep_f1 = pd.read_csv(DATA_DIR / "f1_observations_fastvep.csv")
    wt = pd.read_csv(DATA_DIR / "wall_times.csv")

    PUBLICATION_BINARY_MD5 = "23a294cf"
    S_TO_SUITE = {
        "s01": "clinvar_full_grch37",
        "s02": "clinvar_full_grch38",
        "s03": "gnomad_v2.1.1_chr21",
        "s04": "gnomad_v4.1_chr21",
        "s05": "1kg_integrated_chr21",
        "s06": "1kg_highcov_chr21",
    }
    pub = f1[
        (f1["binary_md5"] == PUBLICATION_BINARY_MD5)
        & (f1["suite"].isin(S_TO_SUITE.keys()))
    ].copy()
    pub["suite"] = pub["suite"].map(S_TO_SUITE)
    pub["order"] = pub["suite"].map({v: i for i, v in enumerate(FIG1_DATASET_ORDER)})
    pub = pub.sort_values("order").reset_index(drop=True)
    if len(pub) != 6:
        raise SystemExit(
            f"ERROR: [generate_figures] panel (a) selected {len(pub)} of 6 vep-rs F1 rows "
            f"at binary {PUBLICATION_BINARY_MD5}; an empty or partial selection still draws "
            "a figure, so this fails instead"
        )

    FASTVEP_SUITE_TO_DISPLAY = {
        "clinvar_full_grch37": "clinvar_full_grch37",
        "clinvar_full_grch38": "clinvar_full_grch38",
        "gnomad_v2.1.1_chr21": "gnomad_v2.1.1_chr21",
        "gnomad_v4.1_chr21": "gnomad_v4.1_chr21",
        "1kg_phase3_chr21": "1kg_integrated_chr21",
        "1kg_highcov_chr21": "1kg_highcov_chr21",
    }
    fv = _fastvep_f1_pin(fastvep_f1, list(FASTVEP_SUITE_TO_DISPLAY.keys()))
    fv["suite"] = fv["suite"].map(FASTVEP_SUITE_TO_DISPLAY)
    fv["order"] = fv["suite"].map({v: i for i, v in enumerate(FIG1_DATASET_ORDER)})
    fv = fv.sort_values("order").reset_index(drop=True)

    # The two structural-variant cells, appended to panel (a) after the six SNP/indel rows
    # and separated from them by their labels. Read from the SV rows of the same two CSVs.
    def _sv(df, sid, col="raw_f1"):
        sub = df[df["suite"] == sid]
        if len(sub) != 1:
            raise SystemExit(
                f"ERROR: [generate_figures] expected exactly 1 {sid} row for panel (a), "
                f"found {len(sub)}"
            )
        return float(sub.iloc[0][col])

    sv_pub = f1[f1["binary_md5"] == PUBLICATION_BINARY_MD5]
    # fastVEP's SV rows are selected by their own sweep pin, like every other selector.
    fv_sv = fastvep_f1[_sweep_pin(fastvep_f1, FASTVEP_SV_F1_SWEEP_ID)]
    sv_rows = [
        # The two CSVs name the SV cells differently, so both spellings are written out
        # rather than derived.
        ("SV per-VCF GRCh37", _sv(sv_pub, "s07"), _sv(fv_sv, "sv_grch37")),
        ("SV per-VCF GRCh38", _sv(sv_pub, "s08"), _sv(fv_sv, "sv_grch38")),
    ]

    # Panels (b) and (c). arch, condition and sweep_id live inside the notes field rather
    # than in columns of their own, so each is matched as a substring there; `arch` is then
    # lifted into a column because both panels group on it.
    _notes = wt["notes"].astype(str)
    base = (
        (wt["cache_state"] == "warm")
        & (wt["run_index"] == "median_of_20")
        & _notes.str.contains("methodology=cross_instance_parallel_clone", na=False)
        & _notes.str.contains("condition=sites_only", na=False)
        & wt["suite"].isin(FIG1_DATASET_ORDER + ["clinvar_grch37", "clinvar_grch38"])
    )
    # wall_times.csv spells the two ClinVar suites without the "_full" infix.
    ALIAS = {"clinvar_grch37": "clinvar_full_grch37", "clinvar_grch38": "clinvar_full_grch38"}
    rows = []
    for engine, sweep in (
        ("vep-rs", VEP_RS_SWEEP_ID),
        ("fastvep", FASTVEP_WALLTIME_SWEEP_ID),
        ("perl", PERL_WALLTIME_SWEEP_ID),
    ):
        sub = wt[base & (wt["engine"] == engine) & _sweep_pin(wt, sweep)].copy()
        sub["arch"] = sub["notes"].astype(str).str.extract(r"arch=([^;]+)")[0]
        sub["suite"] = sub["suite"].replace(ALIAS)
        sub = sub[sub["suite"].isin(FIG1_DATASET_ORDER)]
        _assert_one_row_per_cell(sub, ["engine", "suite", "arch"], f"panel (b) {engine}")
        sub["order"] = sub["suite"].map({v: i for i, v in enumerate(FIG1_DATASET_ORDER)})
        rows.append(sub)
    wall = pd.concat(rows, ignore_index=True)
    for arch in ("arm64", "x86_64"):
        for engine in ("vep-rs", "fastvep", "perl"):
            n = len(wall[(wall["arch"] == arch) & (wall["engine"] == engine)])
            if n != 6:
                raise SystemExit(
                    f"ERROR: [generate_figures] panel (b) has {n} of 6 {engine}/{arch} "
                    "cells; a partial selection still draws a figure"
                )

    ratios, arm, x86 = {}, [], []
    for arch in ("arm64", "x86_64"):
        for suite in FIG1_DATASET_ORDER:
            def _one(engine):
                m = wall[(wall["engine"] == engine) & (wall["arch"] == arch)
                         & (wall["suite"] == suite)]
                return float(m.iloc[0]["wall_time_sec"])
            r = _one("perl") / _one("vep-rs")
            ratios[(suite, arch)] = r
            (arm if arch == "arm64" else x86).append(r)
    arm_geomean = math.exp(sum(math.log(v) for v in arm) / len(arm))
    x86_geomean = math.exp(sum(math.log(v) for v in x86) / len(x86))

    # (a) spans the top row and (b) and (c) sit beneath it; 1x3 would give each panel 59 mm.
    # wspace holds (c)'s y-tick labels, the six dataset names, which are drawn outside its
    # left spine into this gutter: at a narrower gutter they sit on top of (b)'s marks,
    # legible as text and invisible to every structural check. The room comes out of the
    # panels, not the page: a figure authored past the column width is delivered wider and
    # scaled down with its 7 pt type, which `_assert_delivered_width` below catches.
    fig = plt.figure(figsize=(FIG1_WIDTH_IN, FIG1_HEIGHT_IN))
    gs = fig.add_gridspec(2, 2, hspace=0.75, wspace=0.62)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    _fig1_panel_a(ax_a, pub, fv, sv_rows)
    _fig1_panel_b(ax_b, wall)
    _fig1_panel_c(ax_c, ratios, arm_geomean, x86_geomean)

    # Panel labels outside the axes, at the same offset in all three.
    for ax, label in zip((ax_a, ax_b, ax_c), ("(a)", "(b)", "(c)")):
        ax.annotate(
            label, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(-30, 22),
            textcoords="offset points", fontweight="bold", fontsize=9,
        )
    # The gutter settings above move the delivered width, so the limit is asserted rather
    # than assumed; (c)'s legend is measured against its axes for the same reason.
    _assert_delivered_width(fig, "fig1", FIG1_WIDTH_IN * 25.4)
    _assert_legend_within_axes(fig, ax_c, "fig1 (c)")
    return _save(fig, "fig1_concordance_throughput")


FIGURES = {
    "fig1": fig1_concordance_throughput,
    "figS1": figS1_per_class_concordance,
}


def write_manifest(paths: list[Path], partial: bool) -> None:
    """Write MANIFEST.txt, preserving entries for figures this run did not render.

    `--only` renders one figure. Overwriting the manifest with just that figure's
    paths would drop the other figure's digests, so a partial render merges into
    the existing manifest rather than replacing it.
    """
    # TIFF is absent from the manifest: it is gitignored, so a fresh clone has no file
    # for its entry to describe.
    manifest = FIG_DIR / "MANIFEST.txt"
    entries: dict[str, str] = {}
    if partial and manifest.is_file():
        for ln in manifest.read_text().splitlines():
            parts = ln.split("  ")
            if len(parts) == 2:
                entries[parts[1]] = parts[0]
    for p in paths:
        if p.suffix == ".tif":
            continue
        entries[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
    with manifest.open("w") as fh:
        fh.write("# Manuscript figures sha256 manifest\n")
        fh.write("# Generated by manuscript/figures/generate_figures.py\n\n")
        for name in sorted(entries):
            fh.write(f"{entries[name]}  {name}\n")
    print(f"manifest -> {manifest} ({len(entries)} entries)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all", action="store_true", help="Render every figure (default)."
    )
    parser.add_argument(
        "--only",
        choices=list(FIGURES.keys()),
        default=None,
        help="Render only the named figure (for fast iteration).",
    )
    args = parser.parse_args()

    setup_style()
    targets = [args.only] if args.only else list(FIGURES.keys())

    all_paths: list[Path] = []
    for name in targets:
        print(f"rendering {name}...")
        all_paths.extend(FIGURES[name]())

    write_manifest(all_paths, partial=args.only is not None)
    print(f"\nwrote {len(all_paths)} files in {FIG_DIR}")


if __name__ == "__main__":
    main()
