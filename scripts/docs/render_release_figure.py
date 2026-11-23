#!/usr/bin/env python3
"""Render the release's comparison figure, and the values it draws, from the provenance records.

Records. Without --records the README's `## Current release: vX.Y.Z` heading names the release and
the records are the files of `docs/concordance-provenance/` that state it; with --records they are
the files named (a directory stands for its `*.json` files). A record's kind is read from the
record: the release record (`record: release-addendum`), the population record
(`population-measurements`), a comparator's dataset record (`comparator-suites`) and wall-time
record (`comparator-walltime`), and the plugin record (`record_type: plugins`). A record states its
release under `release`, a comparator's dataset record under `scored_with_release.version`. One
release record and one population record are required; the plugin record is optional and without
it panels (c) and (d) are not drawn. Panel (e) is drawn when the release record holds the block
`genome_checks`.

Tools. `comparators.tsv` beside this script lists the tools, one row each: `id`, `name`, `role`
(`released`, `comparator` or `reference`), `record_keys` (the names the records file the tool
under), `colour_slot`, `concordance_panel` and `speed_panel` (`yes` or `no`) and
`headline_extra_args` (the extra arguments of the run a panel reads when the records hold more than
one run of the tool; `-` for none). A tool appears where its row allows it and a record holds its
figure. Each tool's version is the one the records state, and one tool at two versions is refused.

Panels. (a) Raw F1 per dataset on an axis logarithmic in 1 - F1 and reversed, one bar per tool from
the plot's left edge, with a tool's adjusted F1 as a marker: the genome-wide ClinVar datasets (the
release record's `concordance` cells), the whole-genome population datasets
(`population_concordance`, and the population record's `comparator_concordance`), the pool of those
whole-genome datasets from their summed integer counts, as the README pools them, and the
structural-variant sets, each row the counts of its files pooled, as its label says. A comparator's
adjusted counts are the `adjusted` block of its dataset record's cell and, on a population dataset,
its counts less `excluded_perl` and `excluded_rust`; where a record holds none the marker is left
out, with a printed note and one line in the caption. (b) Wall time per dataset and architecture on
a logarithmic axis that starts at the decade under its smallest value: the population record's
`wall_time` sums and, for ClinVar, the `per_cell_aggregates` medians of the release record and of
the comparators' wall-time record. Of the release record's entries of a cell and architecture the
one that names no `target_cpu` is read, and the value's key says so, since the record may hold
others that name one; the axis label and the caption say which lines are medians and
which are sums. The caption calls a cell of the record's `sum_of_cell_medians_sec` a separately
timed part of its dataset, and adds the clause of the population record's `methodology_note` that
says what a dataset's sum is of, where it has one. The caption and the alternative text say in how
many runs a dataset, or a separately timed part of one, is timed, tool by tool where the tools
differ, from the number each record states beside a value: `n_clones` of a ClinVar entry,
`n_clones_per_cell` of a population sum's entry or, where the entry states none, of its engine's
block. Each number has its keys among the values file's `sources`; where a record states none
beside a value, the clause is left out, with a printed note.
No tool is named in the panel: a circle wears its tool's colour, and in the legend a circle in that
colour stands beside the swatch of each tool that has bars too. (c) Raw F1 per gate row of the
plugin record's `rows`, grouped by the kind of data a row is scored on, the counts summed over a
row's entries: a bar per tool at the feature projection with the tool's adjusted marker, as on
panel (a), and a plain cross where the record lists the plugin as `absent` for a comparator, in the
room left of the tick that reads 0, the reason named once in the legend. A comparator's adjusted
counts are the adjusted counts of its entries; where they hold none the marker is left out, with a
printed note and the caption's one line. Where an entry states under `adjusted_undetermined` that
one or more of the comparator's values meet the condition of a declared class that is not applied
(`values_meeting_condition`), its adjusted F1 is not determined on that row: the marker is left
out, with a printed note, and one caption line names the rows, its record keys among the values
file's `sources`; a count of 0 there sets nothing aside, and the marker is drawn from the entry's
adjusted counts. A per-allele figure is drawn on no panel: it stays in the record. Every gate row
so shows each tool the record scores: a row that neither scores a comparator per feature nor lists
the plugin as absent for it is refused, as are two reasons for one comparator's absent plugins. A
row scored on real data sums, per tool, the genome-wide ClinVar entries of its key and the entries
that the record's `population_data` holds for the plugin on the whole-genome population datasets of
the row's assembly, which the population record's `datasets` name: the bar from the summed raw
counts and the marker from the summed adjusted counts, where every summed entry holds them. Such an
entry has the form of a `real_data` entry and names its dataset as the population record does; a
tool is drawn on no dataset the record does not hold for it, and a value's key names every entry
it sums. The record's `population_absent` states, with a reason, each such dataset on which a
plugin has no count. The group of those rows is labelled with the plain names of the datasets the
released tool's rows sum: ClinVar's, then each population dataset's, one that the record labels in
several releases going by the leading words its labels share. (d)
The plugins each tool reads, from the record's `support`: drawn when the record holds a count for
the released tool and every comparator of the registry, and otherwise left out with one line in the
caption; the reference has no row. The panel's title names the tools it draws, and the note of what
the counts are of, with the number of plugins the reference lists, stands under the panel. A panel
carries its letter once: where the plugin panel's rows take both columns, its second part is titled
as the first's continued and carries no letter. (e) The ready-made Ensembl VEP caches that read
exact in the genome check, of those checked, from the release record's `genome_checks`. The released
tool's row is one bar, the block's total of caches that are exact, labelled `<exact> of <checked>`
with thousands separators. Each comparator of the registry that the block's `absent` lists has a
row with no bar: a cross in the tool's colour at the left of the plot, where a bar would start, and
the entry's reason beside it, so the legend's one reason for a tool's crosses stays the plugin
record's. Every value names its key of the block. A comparator with a mark in another panel that
`absent` does not list is refused by name, as is an entry without a statement; an entry stated at
another version than the records state for its tool is refused with the other versions, as one tool
at two versions; an entry whose engine the registry lists as no comparator is noted and left out;
and a block without `absent` draws the released tool's row alone. The panel stands in the second
column under (d) and its note.

Marks. On an F1 axis a value is a place, not a length. The plot's left edge stands a fixed room
left of the tick that reads 0, a marker's half width and four pixels more, and every bar starts at
that edge: a marker so stands whole inside the plot at its counts' place, never over the axis line,
and a bar is seen left of its own marker on every row, at an F1 near 0 too. A bar ends at its
counts' place, however near 0 that is. An F1 of exactly 1 is drawn at the axis floor, whose tick
reads 1 behind a break mark, since 1 has no place on an axis of 1 - F1 in decades; the floor is the
first decade under the smallest 1 - F1 the panels hold, the decade below for one that lies on a
decade, and no nearer 1 than the seventh, so an F1 below 1 stands short of it. The F1 panels label
one set of decades, the ticks that read 0 and 1 always among them. Every row carries its adjusted
F1 as a marker, so that none reads as unscored: a triangle at the floor for an F1 of exactly 1, at
the end of a full bar as beside a shorter one, and a diamond for one below 1. A row at the maximum
carries no number: its full bar and the triangle at its end say it. A marker is drawn whole, over
every bar, with a thin outline in the ground's colour, and no gap is cut into it where it stands
over a bar's tip; a thin dashed line joins the end of its tool's bar to it, at exactly 1 and below
1 alike, for every tool, with a clause in the caption for each F1 panel and no entry in the legend.
The line is drawn wherever a marker stands right of its bar's end and the two do not touch, for a
comparator as for the released tool. A marker that stands left of its bar's end, an adjusted F1
under the raw one, has no line: it stands alone on its bar. A marker touches its bar where its body, inside its outline, stands half a pixel or
more over the bar's end, and there no line is drawn: a marker is set on whole pixels, and a body
that reaches less far over the bar's end may leave ground before it. The line is drawn at any
stretch between the two, however short,
and nothing lengthens it: it starts with a dash at the bar's end and runs through the marker's
outline to the body, where it ends with a dash, so that no ground stands before the body. A
diamond's body narrows towards its bar, so on a pixel row off the lane's line its edge stands
further from the bar's end: past the last dash the line is drawn on, solid and in the body's
colour, far enough into a diamond that it ends on the body on every pixel row it inks. Its
dashes and the gaps between them are fitted alike to its length, no dash shorter than half again
the line's width. A line too short for two such dashes and the gap between them takes a finer
pattern, a shorter first dash before a gap and a last dash of their least lengths, as long as the
first dash is a pixel and a half; a shorter line cannot show two dashes at the image's size and is
a solid stub, and so the caption and the alternative text call the join a line. Every F1 is 2I/(A+B) from the record's integer counts, rounded half
up to nine places, as the README's tables print it. A mark's place on the axis is taken from the
counts and never from the printed figure, and an F1 is exactly 1 only when its counts are equal:
one that prints as 1.000000000 from unequal counts stands short of the floor, as a bar and as a
diamond, and carries no drawn label. Every wall time is rounded as the README prints it.

Texts. The caption says once that on an F1 axis a value is a place, not a length. The caption and the
alternative text say for each F1 panel that a marker is joined to the end of its bar by a line
where the two differ (a line, since a join too short to show two dashes is drawn solid), and once each that a marker standing at the end of its bar means that the raw
and the adjusted F1 are equal as drawn.
Where a comparator has an adjusted marker, the caption quotes what its adjusted figure applies from
the population record's `comparator_concordance.note`, the sentence of it that speaks of the
adjusted figure, and says nothing on it where the record holds no such sentence, nor in a draft,
whose figures the note of a final record may not hold for. It states that the
released tool matches the reference exactly except for the documented divergences only over records
whose every adjusted F1 of the released tool is exactly 1; over any other record it says in that
sentence's place that the render is a draft of a state that is not shipped. The alternative text
says, panel by panel, what the panel holds as it is drawn: its tools and the mark each has, for the
wall-time panel the order of the tools' circles on a row from the left, since that panel tells its
tools apart by colour alone, and for
the plugin panel what its groups of rows are scored on, in the caption's clauses, on how many rows
a comparator has a figure and on how many a cross, and on how
many the released tool's adjusted F1 is exactly 1. Both texts say in one sentence which datasets
each tool's counts cover on the rows scored on real data, from the entries each row sums: the
population datasets a tool's rows sum beside ClinVar, assembly by assembly and each as the
population record labels it, every row of the tool that covers fewer with what it covers, and a
tool whose rows cover ClinVar alone; the sentence has its key among the values file's `sources`.
No panel names a version: the legend alone names
each tool with the version its records state. An F1 axis is labelled F1 against the reference, and
the legend says that a bar is the raw F1. Both texts also say, from the record's own counts and
each with its key among the values file's `sources`, what the image leaves unsaid: the row names of
the plugin panel where they are more than its plugins, the rows on which the tools' figures are
summed over different numbers of the row's sets, each such number kept again by itself with the
entries it counts, and how many of the released tool's adjusted
markers stand where the axis reads 1, past its break; and they state the count each tool reads in
the panel of counts. The reason a comparator's crosses stand for, which both texts give, has one
entry among the `sources` too, and so has each number of rows the alternative text gives of the
plugin panel and, where the order of the tools' circles differs between rows, of the wall-time
panel, with what it counts of the panel's values. A sentence whose record holds nothing to say is
left out. For panel (e) the caption says what the bar counts and against which reference, the
counts of the block's sets under the block's labels, that the caches that are not exact are named
on the supported-genomes page where the block counts one, and for each comparator with a row one
sentence built from the tool's name, the version its entry states and the clauses of the entry's
statements: what the comparator's own documentation says, and that the tool therefore has no
counterpart in the panel. The alternative text says what the panel holds as it is drawn, the bar
with its label and each cross with its reason, then the sets' counts and the same sentence. The
reference, each set's counts and each comparator's sentence have their keys among the `sources`.

Two layers. The values layer uses the standard library alone and writes the values file: `image`,
`canvas`, `alt_text`, `caption`, `notes` (the lines drawn inside the image that are no panel's),
`sources` (each figure or sentence of the two texts and of a note that no panel's value holds, as
`text` with the `record` and `key` it is drawn from; one whose record holds no such key is left out
of the texts), `tools`, `records` (each file with its digest) and `panels` (each `id`, `title`,
`groups` of `label` and `rows` and, where a line is drawn under the panel, its `note`; a row has
`label` and `values`, each value `tool`, `kind` (`bar`, `marker` or `text`), `value`, `record` and
`key`, and an F1 its `counts`: matched, reference and tool; a value that draws on two records names
both, joined by `and`). The same records give the same bytes. The drawing layer needs the packages
of `requirements-figure.txt`; it draws the values file and nothing else on a canvas of 1080 by 1350
pixels and writes the PNG with the values file's digest in a text chunk (`values_sha256`).
--values-only stops after the values file.

Refused, each as one `ERROR: [render_release_figure] ...` line: a missing or unreadable record, two
records of one kind where one is read, a record lacking a field the released tool's marks read, two
wall times of one tool, dataset and architecture, a registry that is malformed or lists more tools
than the figure has fills, a structural-variant cell whose files do not sum to its counts, adjusted
counts of the released tool that are unequal and round to 1 at nine places (the released version's
adjusted F1 is exactly 1 or the record is not a shipped state), a plugin record whose `gate` block
does not say it passed for the released tool, and records that do not all state one release.
--draft draws the last two as a marked draft: a line inside the image, which also opens the caption
and the alternative text, and never under `docs/assets/`. A draft also draws what panel (e) refuses,
each with a printed note: it leaves out the row of a comparator the block does not list or lists
without a statement, and draws an entry stated at another version, naming each version. Of the plugin record's entries and
statements on the population datasets, one on a dataset the population record does not list or of
another assembly than its dataset's is refused, as is a statement beside an entry; a row scored on
real data on which the released tool lacks a population dataset of its assembly without the
record's statement is also refused, and in a draft it is drawn from the datasets the record holds,
with a printed note. With a plugin record, a population dataset whose label holds no word is refused.
A comparator's or the reference's missing
figure loses its mark with a printed note, and so does a tool or a dataset the figure has no place
for: an engine the registry does not list, a dataset on which the comparator alone is scored. The
drawing layer refuses the same way a drawing whose text is small, faint or set over other text,
whose plots have no room on the canvas, whose rows would be set tighter than 18 pixels to fit it
(a figure with more rows takes a taller canvas), whose wall-time circles cannot be set apart on
their line, and one that does not set each panel's letter once. It also refuses a values file whose
tools hold no reference for an F1 axis to name, as where no record holds a wall time of the
reference, and one whose wall-time panel holds no group or no value.

usage: render_release_figure.py [--repo DIR] [--records PATH ...] [--registry FILE]
                                [--values FILE] [--png FILE] [--values-only] [--draft]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import struct
import sys
import zlib
from decimal import ROUND_HALF_UP, Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path

COMPONENT = "render_release_figure"
IMAGE = "docs/assets/comparison.png"
VALUES = "docs/assets/comparison.values.json"
CANVAS = (1080, 1350)
STAMP_KEY = "values_sha256"
HEADING = re.compile(r"^## Current release: v(\d+\.\d+\.\d+)\s*$", re.MULTILINE)

REGISTRY_COLUMNS = (
    "id", "name", "role", "record_keys", "colour_slot", "concordance_panel", "speed_panel", "headline_extra_args",
)
ROLES = ("released", "comparator", "reference")
MAX_TOOLS = 8
# One hue per tool and nothing else encoded in hue: blue, orange and green of the Okabe-Ito palette stay
# apart under the common colour-vision deficiencies, and blue against orange also separates in grayscale.
FILLS = {"1": "#0072B2", "2": "#E69F00", "3": "#009E73"}

KINDS = {
    "release-addendum": "release",
    "population-measurements": "population",
    "comparator-suites": "datasets",
    "comparator-walltime": "walltime",
}
KIND_NAMES = {
    "release": "release record",
    "population": "population record",
    "datasets": "comparator dataset record",
    "walltime": "comparator wall-time record",
    "plugins": "plugin record",
}
CLINVAR_CELLS = (("s01", "ClinVar", "GRCh37"), ("s02", "ClinVar", "GRCh38"))
SV_CELLS = (("s07", "GRCh37"), ("s08", "GRCh38"))
ARCHES = (("arm64", "ARM"), ("x86_64", "x86"))
WHOLE_GENOME_SUFFIX = ", whole genome"
NINE_PLACES = Decimal("0.000000001")

# Per kind of data a gate row is scored on: the plugin record's entry list, the group label and the
# caption's clause. `{plugins}` and `{their}` are filled from the rows of the group. The label of the rows
# scored on real data goes on to name the population datasets the released tool's rows sum.
PLUGIN_GROUPS = (
    ("clinvar", "real_data", "Genome-wide ClinVar", "plugins verified on real data"),
    ("committed_table", "committed_table", "Committed table", "{plugins} on {their} committed table"),
    (
        "synthetic_fixture", "synthetic_fixtures", "Synthetic fixtures (data that cannot be redistributed)",
        "plugins whose data cannot be redistributed on synthetic fixtures, marked synthetic",
    ),
)
CUSTOM_GROUP = (
    "custom_annotation", "Custom annotation (--custom sources)",
    "custom-annotation sources on the genome-wide ClinVar datasets",
)
# A gate row scored on real data sums, per tool, the genome-wide ClinVar entries of its key and the entries the
# record holds for the plugin on the whole-genome population datasets of the row's assembly. Those entries have
# a list of their own, each in the form of a `real_data` entry under the population record's name of its
# dataset; a second list states, with a reason, each such dataset on which a plugin has no count.
REAL_KIND = "clinvar"
POPULATION_LIST = "population_data"
POPULATION_ABSENT = "population_absent"
SYNTHETIC_KIND = "synthetic_fixture"
SYNTHETIC_TAG = "synthetic"
ENTRY_TARGETS = ("dataset", "corpus")
RAW_COUNTS = ("intersection", "perl", "rust")
COMPARATOR_COUNTS = ("intersection", "perl", "engine_tuples")
ADJUSTED_COUNTS = ("adj_intersection", "adj_perl", "adj_rust")
# A comparator entry's statement that its adjusted figure is not determined, and the count inside it: the
# tool's values on the row that meet the condition of a declared class that is not applied. At 1 or more the
# class would set values aside and cannot, so no adjusted F1 stands; at 0 it sets nothing aside.
UNDETERMINED, UNDETERMINED_COUNT = "adjusted_undetermined", "values_meeting_condition"
# Per plugin whose reference ran its module with stand-in files, the caption's clause.
STAND_IN_CLAUSES = {
    "LoFTEE": "LoFTEE's reference ran with stand-ins for the MaxEntScan files its module bundles",
}

PANEL_TITLES = {
    "a": "Concordance per dataset",
    "b": "Wall time per dataset and architecture",
    "c": "Concordance per plugin",
    "d": "Plugins read by {tools}",
    "e": "Ready-made Ensembl VEP caches, exact of checked",
}
PANEL_AXES = {"a": "f1", "b": "seconds", "c": "f1", "d": "count", "e": "count"}
# What a row of each F1 panel is, in the alternative text.
F1_ROWS = {"a": "dataset", "c": "plugin"}
# The groups of panel (a); `{count}` is the number of whole-genome datasets the pooled row pools, in words.
DATASET_GROUPS = (
    "Genome-wide ClinVar", "Whole-genome population datasets", "The {count} whole-genome datasets, pooled", "Structural variants",
)
# The marks of the legend's entries by kind, a row of the legend each: tools, the marks of an F1,
# architectures.
LEGEND_KINDS = (("fill", "fill_dot", "dot"), ("raw", "exact", "adjusted", "absent"), ("filled", "open"))
NUMBER_WORDS = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
COUNT_WORDS = {1: "One panel", 2: "Two panels", 3: "Three panels", 4: "Four panels", 5: "Five panels"}

# What follows a panel's title over its second part, which carries no letter.
CONTINUED = ", continued"
DRAFT_LINE = "Draft: not drawn from one release's final records."
COUNTS_NOTE = (
    "Counts are of the {supported} plugin{plural} {released} supports, of {reference}'s {listed}. Plugins {released} does not "
    "port are not counted; the README lists them by reason."
)
SUPPORT_LEFT_OUT = "The panel of the plugins each tool reads is left out: the records hold no count for {tools}."
# The caption's sentence on what a comparator's adjusted figure applies, in the words of the population
# record's note on it; a record without such a note gives no sentence, and a draft carries none.
ADJUSTED_NOTE = "On a comparator's adjusted F1 the record's note reads: \"{sentence}\""
ADJUSTED_LEFT_OUT = "{tools}'s adjusted F1 is left out where no record holds its adjusted counts."
ADJUSTED_UNDETERMINED = "{tool}'s adjusted F1 is not determined on {rows}."
ADJUSTED_EMPTY = "{tool}'s adjusted F1 is left out on {rows}, where every tuple of the row is set aside and none is left to score."
# The caption's sentence on concordance, which stands only over records whose every adjusted F1 of the released
# tool is exactly 1, and the sentence that stands in its place over any other record.
CONCORDANCE_SENTENCE = "vep-rs matches Ensembl VEP exactly except for the documented divergences."
UNSHIPPED_SENTENCE = "This is a draft of a state that is not shipped: not every adjusted F1 of {released} is exactly 1."
RECORDS_SENTENCE = "Every value is in the records listed under Provenance records."
# The caption's one sentence on how an F1 axis is read: a bar starts in the room left of the tick that reads 0.
# It names no tick: the caption's only marks are bars, markers, lines and crosses.
PLACE_SENTENCE = (
    "On an F1 axis a value is a place, not a length: every bar starts at the left edge of its plot, left of where an F1 of 0 stands, "
    "and ends at the place of its F1."
)
# The wall-time texts' words for a cell of the records, which is what one median is taken over: a ClinVar cell is
# its whole dataset on one architecture, and a population dataset's cells are the parts its sum adds up.
SUMMED_PARTS = "its dataset's separately timed parts"
RUNS_SUBJECT = "each dataset, or each separately timed part of one,"
# The one sentence of the caption and of the alternative text on a marker without a line from its bar:
# the line's length is the distance between the places of the two figures.
AGREE_SENTENCE = "A marker standing at the end of its bar means that the tool's raw and adjusted F1 are equal as drawn."
# Sentences of the caption and the alternative text on what the image leaves unsaid, each from the record's own
# counts: the plugin panel's row names where they are more than its plugins, the rows whose tools are summed
# over different numbers of sets, and the place of an F1 of exactly 1. The last two strings are the keys the
# values file names for the sentences before them.
NAMES_SENTENCE = "The plugin panel's {rows} rows carry {names} row names for {plugins} plugin{plural}: {doubled}."
UNEVEN_SENTENCE = (
    "A row's figure sums the per-feature entries that the record holds for the tool on the row's own dataset or corpus, each scored by "
    "itself, and the tools differ in the number of those entries on {count} row{plural}: {rows}."
)
# The sentence of the caption and of the alternative text on the datasets a tool's counts cover on the plugin
# panel's rows scored on real data, which names the bar by no second word: its clause for a tool whose rows sum
# population datasets, with the builds' datasets and the rows that cover fewer, its clause for a tool whose
# every such row is of ClinVar alone, and the key the values file names for the sentence.
COVERAGE_SENTENCE = "On a row of a plugin verified on real data, {clauses}."
COVERAGE_SUMMED = (
    "{tool}'s counts are summed over the genome-wide ClinVar dataset of the row's genome build and whole-genome population datasets of "
    "that build ({builds}){fewer}"
)
COVERAGE_FEWER = ", except on {rows}"
COVERAGE_ALONE = "{tool}'s counts are of the genome-wide ClinVar dataset of the row's genome build alone"
COVERAGE_KEY = (
    "the per-feature entries of each engine on a row of the plugin panel scored on real data, in real_data and population_data: their datasets, each "
    "named by its label among the population record's datasets"
)
EXACT_SENTENCE = (
    "An F1 of exactly 1 has no place on a scale of 1 - F1: it is drawn where the axis reads 1, past a break in the axis line, and {exact} "
    "of the {markers} adjusted marker{plural} of {released} stand{one} there."
)
# The alternative text's sentence, and its key, on a plugin panel whose groups the layout splits over the two
# columns: read without the image, that text alone says so.
CONTINUES_SENTENCE = "The plugin panel runs over both columns: its groups from {group} on continue in the second column, under the same title."
CONTINUES_KEY = "rows: scored_on (the plugin panel's groups in the order the figure gives them, which the layout splits over the two columns)"
UNEVEN_KEY = "the per-feature entries of each engine on a row of the plugin panel (plugin, assembly, dataset or corpus): their number"
# The key of the reason the caption quotes for a tool's crosses: the record states it once per plugin it lists as
# absent for the tool's engine, and one reason for all of them.
ABSENT_KEY = "absent[engine={engine}].reason (one reason for every plugin the record lists as absent for the engine)"
# What follows the name of a tool's per-feature entries on a row in the key of the number of sets the texts print.
SETS_KEY = "the number of entries"
EXACT_KEY = "the adjusted counts of {released} on every row of the F1 panels: equal counts are an F1 of exactly 1"
# The keys of the numbers of rows the alternative text gives of a panel, each a count of the values file's own rows
# of that panel: the rows that hold a value of a tool and kind, the rows on which a tool's marker has equal counts,
# every row, the rows that hold a marker, and the rows whose markers stand in one order of their tools.
PANEL_ROWS_KEY = "panels[id={panel}]: the rows that hold a value of tool={tool}, kind={kind}: their number"
PANEL_EXACT_KEY = "panels[id={panel}]: the rows whose value of tool={tool}, kind=marker has equal counts: their number"
PANEL_TOTAL_KEY = "panels[id={panel}]: its rows: their number"
PANEL_TIMED_KEY = "panels[id={panel}]: the rows that hold a value of kind=marker: their number"
PANEL_ORDER_KEY = "panels[id={panel}]: the rows whose values of kind=marker, from the smallest, are of {tools}: their number"
RAW_LEGEND = "Bar: raw F1"
JOIN_CLAUSE = "joined to the end of its bar by a line where the two differ"
F1_AXIS_SCALE = "log(1 − F1) scale; further right is better"
SECONDS_AXIS = "Wall time (s), log scale; further left is better"
SECONDS_STATISTICS = "{medians}: medians; other datasets: sums of medians"
COUNT_AXIS = "Plugins read"
# Panel (e), the genomes checked from their ready-made caches: the release record's block, the axis label, the
# clause that opens the panel's part of the caption and of the alternative text, and their sentences. Every
# number, label and name of them is read from the block. A comparator's sentence gives what its own
# documentation says, in the clauses of the block's entry for it, and nothing of a run.
GENOMES_BLOCK = "genome_checks"
GENOMES_AXIS = "Caches exact"
GENOMES_CLAUSE = "the ready-made Ensembl VEP caches that are exact, of those checked"
GENOMES_COUNTED = (
    "The bar counts the caches on which {released}'s output reads exact on the genome check's generated input, of those checked "
    "against the stored output of {reference}."
)
GENOMES_SETS = "By set, exact of checked: {sets}."
GENOMES_NOT_EXACT = "The caches that are not exact are named on the supported-genomes page."
GENOMES_ABSENT = "{tool}'s documentation says {says}, so it has no counterpart in this panel, where a cross stands in its row."
GENOMES_DRAWN = "a bar for {released}, labelled {label}"
GENOMES_CROSSED = "for {tool} a cross ({reason}) and no bar"
GENOMES_SET_KEY = GENOMES_BLOCK + ".sets[label={label}].label, exact, checked"
GENOMES_ABSENT_KEY = GENOMES_BLOCK + ".absent[engine={engine}].stated_by.version, statements[].says"
DRAFT_HELP = "draw a marked draft from a plugin record whose gate block does not read passed for the released tool, or from records of several releases"


def fail(message: str) -> None:
    sys.exit(f"ERROR: [{COMPONENT}] {message}")


def read_registry(path: Path) -> list[dict]:
    """The registry's tools in file order, each {id, name, role, keys, fill, concordance, speed,
    headline}. Refused by name: a missing file, another header, a row of another
    width, more than MAX_TOOLS rows, more rows than the figure has fills, a repeated id or colour slot,
    a cell outside its column's values, and a registry without exactly one released tool and one
    reference."""
    if not path.is_file():
        fail(f"no registry {path}")
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    header = tuple(lines[0].split("\t")) if lines else ()
    if header != REGISTRY_COLUMNS:
        fail(f"{path.name} opens with the columns {', '.join(header) or 'none'}; the registry's are {', '.join(REGISTRY_COLUMNS)}")
    rows = lines[1:]
    if len(rows) > MAX_TOOLS:
        fail(
            f"{path.name} lists {len(rows)} tools; past {MAX_TOOLS} a comparison is drawn as small multiples, one panel per "
            "tool, which this figure is not"
        )
    if len(rows) > len(FILLS):
        fail(f"{path.name} lists {len(rows)} tools and the figure holds {len(FILLS)} fills, one per tool; a tool past them has no colour of its own")
    tools = []
    for number, line in enumerate(rows, start=2):
        cells = line.split("\t")
        where = f"{path.name} line {number}"
        if len(cells) != len(REGISTRY_COLUMNS):
            fail(f"{where} has {len(cells)} cells; a row has {len(REGISTRY_COLUMNS)}")
        row = dict(zip(REGISTRY_COLUMNS, (cell.strip() for cell in cells)))
        for column, allowed in (
            ("role", ROLES), ("colour_slot", tuple(FILLS)), ("concordance_panel", ("yes", "no")),
            ("speed_panel", ("yes", "no")),
        ):
            if row[column] not in allowed:
                fail(f"{where} carries {column} {row[column]!r}; it is one of {', '.join(allowed)}")
        keys = tuple(row["record_keys"].split())
        if not row["id"] or not row["name"] or not keys:
            fail(f"{where} names no id, name or record key")
        tools.append({
            "id": row["id"], "name": row["name"], "role": row["role"], "keys": keys, "fill": FILLS[row["colour_slot"]],
            "concordance": row["concordance_panel"] == "yes", "speed": row["speed_panel"] == "yes",
            "headline": "" if row["headline_extra_args"] == "-" else " ".join(row["headline_extra_args"].split()),
        })
    for field in ("id", "fill"):
        seen = [tool[field] for tool in tools]
        if len(set(seen)) != len(seen):
            fail(f"{path.name} repeats {'an id' if field == 'id' else 'a colour slot'}")
    for role in ("released", "reference"):
        count = sum(tool["role"] == role for tool in tools)
        if count != 1:
            fail(f"{path.name} lists {count} tools of the role {role}; the figure reads one")
    return tools


def tool_of(tools: list[dict], key) -> dict | None:
    """The registry tool the records file under `key`, or None."""
    return next((tool for tool in tools if key in tool["keys"]), None)


def repo_path(path: Path, repo: Path) -> str:
    """A record's path under the repository, or its file name when it lies outside."""
    try:
        return path.resolve().relative_to(repo.resolve()).as_posix()
    except ValueError:
        return path.name


def load_record(path: Path, repo: Path) -> dict:
    """A record as {path, file, sha256, data, kind}; `kind` is None for a file the figure does not read, one
    whose JSON is no object among them."""
    try:
        raw = path.read_bytes()
    except OSError:
        fail(f"cannot read {path}")
    try:
        data = json.loads(raw)
    except ValueError:
        fail(f"{path.name} is not JSON")
    kind = None
    if isinstance(data, dict):
        named = data.get("record")
        kind = "plugins" if data.get("record_type") == "plugins" else KINDS.get(named) if isinstance(named, str) else None
    return {"path": path, "file": repo_path(path, repo), "sha256": hashlib.sha256(raw).hexdigest(), "data": data, "kind": kind}


def stated_release(data: dict):
    """The release a record states: its `release`, else the version under `scored_with_release`."""
    if data.get("release") is not None:
        return data["release"]
    scored = data.get("scored_with_release")
    return scored.get("version") if isinstance(scored, dict) else None


def current_release(readme: Path) -> str:
    """The version the README's current-release heading names."""
    if not readme.is_file():
        fail(f"no {readme}")
    found = HEADING.findall(readme.read_text(encoding="utf-8"))
    if len(found) != 1:
        fail(f"{readme.name} has {len(found)} '## Current release: vX.Y.Z' headings; one is read")
    return found[0]


def discover(repo: Path, paths: list[Path] | None) -> list[dict]:
    """The records the figure reads. Named files must each be a record the figure reads; the members of a
    named directory, and of the provenance directory when nothing is named, are kept when they are one, and
    the provenance directory's also only when they state the README's release."""
    release = None
    if paths is None:
        release = f"v{current_release(repo / 'README.md')}"
        paths = [repo / "docs" / "concordance-provenance"]
    records = []
    for path in paths:
        if path.is_dir():
            members = [load_record(member, repo) for member in sorted(path.glob("*.json"))]
            records += [record for record in members if record["kind"] is not None]
        elif path.is_file():
            record = load_record(path, repo)
            if record["kind"] is None:
                fail(f"{path.name} is no record the figure reads: a release, population, comparator or plugin record")
            records.append(record)
        else:
            fail(f"no {path}")
    if release is not None:
        records = [record for record in records if stated_release(record["data"]) == release]
    return records


def extra_args(record: dict, tool: dict) -> str | None:
    """The extra arguments a comparator record states for `tool` (`""` for a run given none), or None when the
    record holds no run of the tool."""
    data = record["data"]
    if record["kind"] == "datasets":
        block = data if data.get("engine") in tool["keys"] else None
    else:
        engines = data.get("engines") if isinstance(data.get("engines"), dict) else {}
        block = next((engines[key] for key in tool["keys"] if isinstance(engines.get(key), dict)), None)
    if block is None:
        return None
    return " ".join(str(block.get("engine_extra_args") or "").split())


def one_run_per_tool(records: list[dict], kind: str, tools: list[dict]) -> dict[str, dict]:
    """Per tool id, the one record of `kind` a panel reads. Where the records hold more than one run of a tool
    the registry names a headline run for, the run with those extra arguments is read and the others are set
    aside for every tool; a tool left with two records is refused by name."""
    candidates = [record for record in records if record["kind"] == kind]
    set_aside: list[int] = []
    for tool in tools:
        runs = [record for record in candidates if extra_args(record, tool) is not None]
        if len(runs) > 1 and tool["headline"]:
            headline = [record for record in runs if extra_args(record, tool) == tool["headline"]]
            if len(headline) != 1:
                fail(
                    f"{len(runs)} {KIND_NAMES[kind]}s hold a run of {tool['name']} ({', '.join(r['file'] for r in runs)}) and "
                    f"{len(headline)} of them ran it with {tool['headline']!r}; a panel reads one run"
                )
            set_aside += [id(record) for record in runs if record is not headline[0]]
    chosen = {}
    for tool in tools:
        runs = [r for r in candidates if id(r) not in set_aside and extra_args(r, tool) is not None]
        if len(runs) > 1:
            fail(f"{len(runs)} {KIND_NAMES[kind]}s hold a run of {tool['name']} ({', '.join(r['file'] for r in runs)}); a panel reads one run")
        if runs:
            chosen[tool["id"]] = runs[0]
    return chosen


def select(records: list[dict], tools: list[dict], draft: bool) -> dict:
    """The records by the part they play: {release, population, plugins (or None), datasets and walltime (per
    tool id), read (every record a panel reads, in a fixed order)}. Refused by name: no record or several of a
    kind read once, records that do not all state one release, and a plugin record that did not pass its gate
    for the released tool; --draft admits the last two."""
    chosen: dict = {}
    for kind, required in (("release", True), ("population", True), ("plugins", False)):
        found = [record for record in records if record["kind"] == kind]
        if len(found) > 1:
            fail(f"{len(found)} {KIND_NAMES[kind]}s are given ({', '.join(r['file'] for r in found)}); the figure reads one")
        if required and not found:
            fail(f"no {KIND_NAMES[kind]} is among the records")
        chosen[kind] = found[0] if found else None
    chosen["datasets"] = one_run_per_tool(records, "datasets", tools)
    chosen["walltime"] = one_run_per_tool(records, "walltime", tools)
    read = [chosen["release"], chosen["population"]]
    for kind in ("datasets", "walltime"):
        for record in chosen[kind].values():
            if all(record is not other for other in read):
                read.append(record)
    if chosen["plugins"] is not None:
        read.append(chosen["plugins"])
    chosen["read"] = read

    releases = {record["file"]: stated_release(record["data"]) for record in read}
    for name, release in releases.items():
        if release is not None and not isinstance(release, str):
            fail(f"{name} states a release that is no text: {release!r}")
    if (len(set(releases.values())) != 1 or None in releases.values()) and not draft:
        stated = "; ".join(f"{name} states {release or 'none'}" for name, release in releases.items())
        fail(f"the records do not all state one release ({stated}); --draft draws a marked draft from them")
    released = next(tool for tool in tools if tool["role"] == "released")
    plugins = chosen["plugins"]
    if plugins is not None and not draft:
        gate = plugins["data"].get("gate") if isinstance(plugins["data"].get("gate"), dict) else {}
        if gate.get("passed") is not True or gate.get("engine") not in released["keys"]:
            fail(
                f"{plugins['file']} did not pass its gate for {released['name']} (gate engine {gate.get('engine')!r}, passed "
                f"{gate.get('passed')!r}); --draft draws a marked draft from it"
            )
    return chosen


def block_of(parent, key: str, where: str, kind: type = dict):
    """`parent[key]` as an object or a list, refused by name when absent or of another shape."""
    value = parent.get(key) if isinstance(parent, dict) else None
    if not isinstance(value, kind):
        fail(f"{where} carries no {key}")
    return value


def count(block, key: str, where: str) -> int:
    """A record's tuple count: a non-negative integer, refused by name otherwise."""
    value = block.get(key) if isinstance(block, dict) else None
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        fail(f"{where} carries no count under {key}: {value!r}")
    return value


def f1(intersection: int, a: int, b: int, where: str) -> str:
    """2I/(A+B) at nine decimal places, rounded half up; an exact 1 reads 1.000000000."""
    if a + b == 0:
        fail(f"{where} counts no tuple on either side; an F1 needs one")
    return format((Decimal(2 * intersection) / Decimal(a + b)).quantize(NINE_PLACES, rounding=ROUND_HALF_UP), "f")


def released_adjusted_f1(matched: int, reference: int, engine: int, where: str, released: str) -> str:
    """The released tool's adjusted F1 at nine places, as the README prints it. Counts that are unequal and
    round to 1 are refused by name: the figure would print an exact 1 for a row that is not exact."""
    value = f1(matched, reference, engine, where)
    if Decimal(value) == 1 and not matched == reference == engine:
        fail(
            f"{where} adjusted counts are unequal and their F1 rounds to 1 at nine places (reference {reference}, {released} "
            f"{engine}, matched {matched}); the released version's adjusted F1 is exactly 1 or the record is not a shipped state"
        )
    return value


def seconds(value, places: int, where: str) -> str:
    """A wall time rounded half up to `places` decimals, from the number's decimal text; refused by name
    unless it is finite and positive."""
    try:
        number = Decimal(str(value))
    except InvalidOperation:
        fail(f"{where} is not a number: {value!r}")
    if isinstance(value, bool) or not number.is_finite() or number <= 0:
        fail(f"{where} is not a positive number: {value!r}")
    return format(number.quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP), "f")


def short_label(label: str) -> str:
    return str(label).replace(WHOLE_GENOME_SUFFIX, "")


def and_list(items: list[str]) -> str:
    if len(items) <= 1:
        return "".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


def image_release(block, where: str) -> str:
    """The Ensembl VEP release a reference block's image tag names: the tag without its `release_` prefix."""
    image = block.get("image") if isinstance(block, dict) else None
    if not isinstance(image, str) or ":" not in image:
        fail(f"{where} names no image tag")
    return image.rsplit(":", 1)[1].removeprefix("release_")


def version_text(tool: dict, version: str) -> str:
    """A version as the figure prints it: a reference's without its `release_` prefix, any other tool's behind
    one `v`."""
    return version.removeprefix("release_") if tool["role"] == "reference" else f"v{version.removeprefix('v')}"


class Reader:
    """What the panels share while the records are read: the tools, the chosen records, the versions each
    tool's records state and the notes printed for marks a record does not hold."""

    def __init__(self, tools: list[dict], chosen: dict, draft: bool):
        self.tools = tools
        self.chosen = chosen
        self.draft = draft
        self.released = next(tool for tool in tools if tool["role"] == "released")
        self.reference = next(tool for tool in tools if tool["role"] == "reference")
        self.versions: dict[str, list[str]] = {}
        self.drawn: set[str] = set()
        self.printed: list[str] = []
        self.unadjusted: list[str] = []
        self.undetermined: dict[str, list[tuple[str, str]]] = {}
        self.emptied: dict[str, list[tuple[str, str]]] = {}
        self.uneven: list[tuple[str, dict[str, int], dict[str, str]]] = []
        self.population: dict[str, tuple[str, str]] = {}
        self.coverage: list[tuple[dict, str, str, list[str]]] = []
        self.real_label = ""
        self.row_names: list[str] = []
        self.sources: list[dict] = []
        self.runs: dict[str, list[tuple[int | None, str, str]]] = {}
        self.genomes: tuple[dict, dict[str, dict] | None] | None = None
        self.genome_sentences: dict | None = None

    def note(self, text: str) -> None:
        if text not in self.printed:
            self.printed.append(text)

    def run_count(self, tool: dict, runs, record: dict, key: str) -> None:
        """Keep the number of runs a record states for the cell, or for each cell, of a wall-time value of
        `tool`, with its record and key; where the record states none there, keep that and print a note."""
        stated = runs if isinstance(runs, int) and not isinstance(runs, bool) and runs > 0 else None
        if stated is None:
            self.note(f"note: {record['file']} states no number of runs beside a wall time of {tool['name']}; the texts leave the number of runs out")
        if (stated, record["file"], key) not in self.runs.setdefault(tool["id"], []):
            self.runs[tool["id"]].append((stated, record["file"], key))

    def version(self, tool: dict, version, where: str) -> None:
        """Record a version a record states for `tool`: a reference's without its `release_` prefix, any
        other tool's behind one `v`."""
        if not isinstance(version, str) or not version:
            fail(f"{where} states no version of {tool['name']}")
        text = version_text(tool, version)
        seen = self.versions.setdefault(tool["id"], [])
        if text not in seen:
            seen.append(text)

    def value(self, tool: dict, kind: str, text: str, record: dict, key: str, counts: list[int] | None = None) -> dict:
        """One value of a row; an F1 carries the integer counts it is computed from as (matched, reference,
        tool), which give its mark's place on the axis."""
        self.drawn.add(tool["id"])
        value = {"tool": tool["id"], "kind": kind, "value": text, "record": record["file"], "key": key}
        if counts is not None:
            value["counts"] = list(counts)
        return value

    def adjusted(self, tool: dict, counts: list[int], where: str) -> str:
        """A tool's adjusted F1 from its adjusted counts: the released tool's is held to the rule of a shipped
        state."""
        if tool["role"] == "released":
            return released_adjusted_f1(*counts, where, tool["name"])
        return f1(*counts, where)

    def comparator_marker(self, values: list[dict], tool: dict, adjusted: list[int] | None, where: str, record: dict, key: str,
                          label: str) -> None:
        """Add a comparator's adjusted marker to a row's values where its record holds the adjusted counts;
        where it holds none, print a note and keep the tool's name for the caption's line."""
        if adjusted is None:
            self.note(f"note: no record holds {tool['name']}'s adjusted counts on {label}; the figure omits that marker")
            if tool["name"] not in self.unadjusted:
                self.unadjusted.append(tool["name"])
            return
        values.append(self.value(tool, "marker", self.adjusted(tool, adjusted, where), record, key, adjusted))

    def undetermined_marker(self, tool: dict, label: str, record: dict, key: str) -> None:
        """Leave out a comparator's adjusted marker on a row whose record states that its adjusted figure is
        not determined: print a note and keep the row and its key for the caption's line."""
        self.note(f"note: {record['file']} states that {tool['name']}'s adjusted F1 is not determined on {label}; the figure omits that marker")
        self.undetermined.setdefault(tool["name"], []).append((label, key))

    def emptied_marker(self, tool: dict, label: str, record: dict, key: str) -> None:
        """Leave out a tool's adjusted marker on a row whose adjusted counts hold no tuple on either side: every
        tuple is set aside and no F1 is left. Print a note and keep the row and its key for the caption's line."""
        self.note(f"note: {record['file']} sets aside every tuple of {tool['name']} on {label}; the figure omits that marker")
        self.emptied.setdefault(tool["name"], []).append((label, key))

    def source(self, text: str, record: dict, key: str) -> str:
        """Keep, for the values file, the record and key a figure or sentence of the caption, of the alternative
        text or of a note is drawn from, where no panel's value holds it; returns the text."""
        self.sources.append({"text": text, "record": record["file"], "key": key})
        return text

    def label(self, tool: dict) -> str:
        """A tool's name and the version, or versions, its records state."""
        return f"{tool['name']} {' and '.join(self.versions[tool['id']])}"


def concordance_tools(reader: Reader) -> list[dict]:
    return [tool for tool in reader.tools if tool["concordance"] and tool["role"] != "reference"]


def read_versions(reader: Reader) -> None:
    """Every version the records state for a tool, in panel order; a block naming an engine the registry does
    not list is noted and left out."""
    chosen = reader.chosen
    release, population = chosen["release"], chosen["population"]
    reader.version(reader.released, release["data"].get("release"), release["file"])
    for family, block in block_of(release["data"], "perl_reference", release["file"]).items():
        reader.version(reader.reference, image_release(block, f"{release['file']} perl_reference.{family}"), release["file"])
    for tool_id, record in chosen["datasets"].items():
        tool = next(tool for tool in reader.tools if tool["id"] == tool_id)
        reader.version(tool, record["data"].get("engine_version"), record["file"])
    comparator = population["data"].get("comparator_concordance")
    if isinstance(comparator, dict):
        tool = tool_of(reader.tools, comparator.get("engine"))
        if tool is None:
            reader.note(f"note: {population['file']} scores {comparator.get('engine')}, which the registry does not list; the figure omits it")
        else:
            reader.version(tool, comparator.get("engine_version"), f"{population['file']} comparator_concordance")
    for key, block in block_of(population["data"], "wall_time", population["file"]).items():
        tool = tool_of(reader.tools, key)
        if tool is None:
            reader.note(f"note: {population['file']} times {key}, which the registry does not list; the figure omits it")
        else:
            reader.version(tool, block.get("engine_version") if isinstance(block, dict) else None, f"{population['file']} wall_time.{key}")
    for record in {id(r): r for r in chosen["walltime"].values()}.values():
        for key, block in block_of(record["data"], "engines", record["file"]).items():
            tool = tool_of(reader.tools, key)
            if tool is None:
                reader.note(f"note: {record['file']} times {key}, which the registry does not list; the figure omits it")
            else:
                reader.version(tool, block.get("version") if isinstance(block, dict) else None, f"{record['file']} engines.{key}")
    plugins = chosen["plugins"]
    if plugins is not None:
        for index, block in enumerate(block_of(plugins["data"], "engines", plugins["file"], list)):
            tool = tool_of(reader.tools, block.get("id") if isinstance(block, dict) else None)
            if tool is None:
                reader.note(f"note: {plugins['file']} scores {block.get('id') if isinstance(block, dict) else None}, which the registry does not list; the figure omits it")
            else:
                reader.version(tool, block.get("engine_version"), f"{plugins['file']} engines[{index}]")
        for index, block in enumerate(block_of(plugins["data"], "perl_reference", plugins["file"], list)):
            reader.version(reader.reference, block.get("ensembl_vep") if isinstance(block, dict) else None, f"{plugins['file']} perl_reference[{index}]")
    read_genomes(reader)
    if not reader.draft:
        for tool in reader.tools:
            versions = reader.versions.get(tool["id"], [])
            if len(versions) > 1:
                fail(f"the records name {tool['name']} at {and_list(versions)}; a tool is drawn at one version, and --draft draws a marked draft")


def read_genomes(reader: Reader) -> None:
    """The release record's block of the genomes checked, where it holds one, and the versions it states: the
    reference's, from its image tag, and each comparator's, from the entry that says why the comparator has no
    counterpart. The block's `absent` entries are kept by tool id, or None where the block holds no such list.
    An entry names its engine once; one whose engine the registry lists as no comparator is noted and left
    out. A version stated here is a version the records state: one tool at two versions is refused with the
    others."""
    release = reader.chosen["release"]
    block = release["data"].get(GENOMES_BLOCK)
    if block is None:
        return
    where = f"{release['file']} {GENOMES_BLOCK}"
    if not isinstance(block, dict):
        fail(f"{where} is no object")
    reader.version(reader.reference, image_release(block_of(block, "reference", where), f"{where}.reference"), f"{where}.reference")
    entries: dict[str, dict] | None = None
    if "absent" in block:
        entries = {}
        for position, item in enumerate(block_of(block, "absent", where, list)):
            at = f"{where}.absent[{position}]"
            engine = item.get("engine") if isinstance(item, dict) else None
            if not isinstance(engine, str) or not engine:
                fail(f"{at} names no engine")
            tool = tool_of(reader.tools, engine)
            if tool is None or tool["role"] != "comparator":
                reader.note(f"note: {at} names {engine}, which the registry lists as no comparator; the figure omits it")
                continue
            if tool["id"] in entries:
                fail(f"{where}.absent names {tool['name']} twice; a comparator has one entry")
            reason = item.get("reason")
            if not isinstance(reason, str) or not reason:
                fail(f"{at} states no reason")
            stated = block_of(item, "stated_by", at)
            reader.version(tool, stated.get("version"), f"{at}.stated_by")
            entries[tool["id"]] = {"engine": engine, "reason": reason, "stated_by": stated, "where": at}
    reader.genomes = (block, entries)


def check_one_configuration(reader: Reader) -> None:
    """A comparator's inputs name one set of extra arguments, so its marks of every panel come from one
    configuration; refused by name otherwise."""
    population = reader.chosen["population"]
    for tool in reader.tools:
        named: dict[str, str] = {}
        wall = block_of(population["data"], "wall_time", population["file"])
        for key in tool["keys"]:
            if isinstance(wall.get(key), dict):
                named[f"{population['file']} wall_time.{key}"] = " ".join(str(wall[key].get("engine_extra_args") or "").split())
        comparator = population["data"].get("comparator_concordance")
        if isinstance(comparator, dict) and comparator.get("engine") in tool["keys"]:
            named[f"{population['file']} comparator_concordance"] = " ".join(str(comparator.get("engine_extra_args") or "").split())
        for kind in ("datasets", "walltime"):
            record = reader.chosen[kind].get(tool["id"])
            if record is not None:
                named[record["file"]] = extra_args(record, tool)
        if len(set(named.values())) > 1:
            stated = "; ".join(f"{where}: {words or 'none'}" for where, words in named.items())
            fail(f"the records ran {tool['name']} with different extra arguments ({stated}); its marks come from one configuration")


def pool_counts(pool: dict | None, tool: dict, raw: list[int], adjusted: list[int] | None = None) -> None:
    """Add one row's counts (matched, reference, engine) of a tool to the pool of the rows it is drawn on,
    and count the rows that gave adjusted counts."""
    if pool is None:
        return
    sums = pool.setdefault(tool["id"], {"rows": 0, "raw": [0, 0, 0], "adjusted_rows": 0, "adjusted": [0, 0, 0]})
    sums["rows"] += 1
    sums["raw"] = [a + b for a, b in zip(sums["raw"], raw)]
    if adjusted is not None:
        sums["adjusted_rows"] += 1
        sums["adjusted"] = [a + b for a, b in zip(sums["adjusted"], adjusted)]


def dataset_cell(reader: Reader, cell: str, label: str, pool: dict | None = None) -> dict:
    """One row of panel (a) from a cell of the release record's `concordance`: the released tool's raw bar and
    adjusted marker, and a comparator's raw bar from its dataset record, with its adjusted marker where the
    cell holds an `adjusted` block. The counts go into `pool` when given."""
    release = reader.chosen["release"]
    values = []
    for tool in concordance_tools(reader):
        if tool["role"] == "released":
            where = f"{release['file']} concordance.{cell}"
            counts = block_of(block_of(release["data"], "concordance", release["file"]), cell, f"{release['file']} concordance")
            raw = [count(counts, key, where) for key in RAW_COUNTS]
            adjusted = [count(counts, key, where) for key in ADJUSTED_COUNTS]
            pool_counts(pool, tool, raw, adjusted)
            values.append(reader.value(tool, "bar", f1(*raw, where), release, f"concordance.{cell}: {', '.join(RAW_COUNTS)}", raw))
            values.append(reader.value(
                tool, "marker", reader.adjusted(tool, adjusted, where), release, f"concordance.{cell}: {', '.join(ADJUSTED_COUNTS)}", adjusted,
            ))
            continue
        record = reader.chosen["datasets"].get(tool["id"])
        counts = record["data"].get("concordance", {}).get(cell) if record is not None else None
        if not isinstance(counts, dict):
            reader.note(f"note: no record holds {tool['name']}'s concordance on {label}; the figure omits that bar")
            continue
        where = f"{record['file']} concordance.{cell}"
        raw = [count(counts, key, where) for key in COMPARATOR_COUNTS]
        held = counts.get("adjusted")
        adjusted = [count(held, key, f"{where}.adjusted") for key in COMPARATOR_COUNTS] if isinstance(held, dict) else None
        pool_counts(pool, tool, raw, adjusted)
        values.append(reader.value(tool, "bar", f1(*raw, where), record, f"concordance.{cell}: {', '.join(COMPARATOR_COUNTS)}", raw))
        reader.comparator_marker(
            values, tool, adjusted, f"{where}.adjusted", record, f"concordance.{cell}.adjusted: {', '.join(COMPARATOR_COUNTS)}", label,
        )
    return {"label": label, "values": values}


def population_row(reader: Reader, released_row: dict, comparator_rows: dict, key: str, label: str, pool: dict) -> dict:
    """One row of panel (a) from a population dataset: the released tool from the release record's
    `population_concordance`, a comparator from the population record's `comparator_concordance`, with its
    adjusted marker where its row states `excluded_perl` and `excluded_rust`. The counts go into `pool`."""
    release, population = reader.chosen["release"], reader.chosen["population"]
    values = []
    for tool in concordance_tools(reader):
        if tool["role"] == "released":
            where = f"{release['file']} population_concordance.{key}"
            matched, perl, engine = (count(released_row, name, where) for name in ("intersection", "perl_tuples", "vep_rs_tuples"))
            set_aside = [count(released_row, name, where) if name in released_row else 0 for name in ("excluded_perl", "excluded_rust")]
            raw, adjusted = [matched, perl, engine], [matched, perl - set_aside[0], engine - set_aside[1]]
            pool_counts(pool, tool, raw, adjusted)
            values.append(reader.value(
                tool, "bar", f1(*raw, where), release, f"population_concordance.{key}: intersection, perl_tuples, vep_rs_tuples", raw,
            ))
            values.append(reader.value(
                tool, "marker", reader.adjusted(tool, adjusted, where), release,
                f"population_concordance.{key}: intersection, perl_tuples - excluded_perl, vep_rs_tuples - excluded_rust", adjusted,
            ))
            continue
        row = comparator_rows.get(tool["id"])
        if row is None:
            reader.note(f"note: no record holds {tool['name']}'s concordance on {label}; the figure omits that bar")
            continue
        where = f"{population['file']} comparator_concordance.{key}"
        raw = [count(row, name, where) for name in ("intersection", "perl_tuples", "engine_tuples")]
        adjusted = None
        if "excluded_perl" in row and "excluded_rust" in row:
            adjusted = [raw[0], raw[1] - count(row, "excluded_perl", where), raw[2] - count(row, "excluded_rust", where)]
        pool_counts(pool, tool, raw, adjusted)
        values.append(reader.value(
            tool, "bar", f1(*raw, where), population, f"comparator_concordance.{key}: intersection, perl_tuples, engine_tuples", raw,
        ))
        reader.comparator_marker(
            values, tool, adjusted, where, population,
            f"comparator_concordance.{key}: intersection, perl_tuples - excluded_perl, engine_tuples - excluded_rust", label,
        )
    return {"label": label, "values": values}


def pooled_row(reader: Reader, pool: dict, rows: int) -> dict:
    """The pool of panel (a)'s whole-genome rows, as the README pools them: each tool's F1 from its counts summed
    over the ClinVar cells and the population datasets, its adjusted F1 from its summed adjusted counts. A
    comparator without a figure on every one of the rows loses its bar with a printed note, and one without
    adjusted counts on every one of them its marker."""
    release, population = reader.chosen["release"], reader.chosen["population"]
    label = f"All {NUMBER_WORDS.get(rows, rows)}"
    cells = ", ".join(f"concordance.{cell}" for cell, *_ in CLINVAR_CELLS)
    values = []
    for tool in concordance_tools(reader):
        sums = pool.get(tool["id"])
        if sums is None or sums["rows"] != rows:
            reader.note(f"note: no record holds {tool['name']}'s concordance on every dataset of {label}; the figure omits that bar")
            continue
        if tool["role"] == "released":
            where = f"{release['file']} population_concordance.per_dataset and {cells}"
            values.append(reader.value(
                tool, "bar", f1(*sums["raw"], where), release,
                f"summed over population_concordance.per_dataset (intersection, perl_tuples, vep_rs_tuples) and {cells} ({', '.join(RAW_COUNTS)})",
                sums["raw"],
            ))
            values.append(reader.value(
                tool, "marker", reader.adjusted(tool, sums["adjusted"], where), release,
                "summed over population_concordance.per_dataset (intersection, perl_tuples - excluded_perl, vep_rs_tuples - excluded_rust) "
                f"and {cells} ({', '.join(ADJUSTED_COUNTS)})", sums["adjusted"],
            ))
            continue
        scored = reader.chosen["datasets"][tool["id"]]
        where = f"{population['file']} comparator_concordance.per_dataset and {scored['file']} {cells}"
        # The comparator's pool draws on two records, and its source names both.
        both = {"file": f"{population['file']} and {scored['file']}"}
        values.append(reader.value(
            tool, "bar", f1(*sums["raw"], where), both,
            "summed over comparator_concordance.per_dataset (intersection, perl_tuples, engine_tuples) and, in "
            f"{scored['file']}, {cells} (intersection, perl, engine_tuples)", sums["raw"],
        ))
        reader.comparator_marker(
            values, tool, sums["adjusted"] if sums["adjusted_rows"] == rows else None, where, both,
            "summed over comparator_concordance.per_dataset (intersection, perl_tuples - excluded_perl, engine_tuples - excluded_rust) "
            f"and, in {scored['file']}, the adjusted blocks of {cells} (intersection, perl, engine_tuples)", label,
        )
    return {"label": label, "values": values}


def pooled_files(release: dict, concordance: dict, cell: str) -> str:
    """The name of a structural-variant row, which says what its figure is: the cell's counts pool its files,
    and the name says so with their number once the files' own counts are seen to sum to the cell's. A cell
    whose files sum to other counts is refused by name; one without a list of files is named without one."""
    where = f"{release['file']} concordance.{cell}"
    block = block_of(concordance, cell, f"{release['file']} concordance")
    files = block.get("per_file")
    if not isinstance(files, dict) or not files:
        return "SV per-VCF"
    pooled = [count(block, key, where) for key in RAW_COUNTS]
    summed = [sum(count(entry, key, f"{where}.per_file") for entry in files.values()) for key in ("intersection", "perl_tuples", "engine_tuples")]
    if summed != pooled:
        fail(
            f"{where} counts {pooled[0]} matched of {pooled[1]} and {pooled[2]}, and its {len(files)} files sum to {summed[0]} of "
            f"{summed[1]} and {summed[2]}; the row's figure is the files' counts pooled"
        )
    return f"SV per-VCF, {len(files)} files pooled"


def panel_datasets(reader: Reader) -> list[dict]:
    """Panel (a): the genome-wide ClinVar datasets, the whole-genome population datasets, the pool of both,
    whose group says in full how many whole-genome datasets it pools, and the structural-variant sets."""
    release, population = reader.chosen["release"], reader.chosen["population"]
    concordance = block_of(release["data"], "concordance", release["file"])
    pool: dict = {}
    clinvar = [dataset_cell(reader, cell, f"{name} {assembly}", pool) for cell, name, assembly in CLINVAR_CELLS]

    scored = block_of(release["data"], "population_concordance", release["file"])
    comparator = population["data"].get("comparator_concordance")
    comparator_tool = tool_of(reader.tools, comparator.get("engine")) if isinstance(comparator, dict) else None
    if comparator_tool is not None and not (comparator_tool["concordance"] and comparator_tool["role"] == "comparator"):
        comparator_tool = None
    rows = []
    for row in block_of(scored, "per_dataset", f"{release['file']} population_concordance", list):
        dataset = row.get("dataset") if isinstance(row, dict) else None
        if not isinstance(dataset, str) or not row.get("label") or not row.get("assembly"):
            fail(f"{release['file']} population_concordance.per_dataset holds an entry without dataset, label and assembly")
        theirs = {}
        if comparator_tool is not None:
            match = [r for r in block_of(comparator, "per_dataset", f"{population['file']} comparator_concordance", list) if isinstance(r, dict) and r.get("dataset") == dataset]
            if len(match) == 1:
                theirs[comparator_tool["id"]] = match[0]
        rows.append(population_row(reader, row, theirs, f"per_dataset[dataset={dataset}]", f"{short_label(row['label'])} {row['assembly']}", pool))
    if comparator_tool is not None:
        # The rows are the release record's: a dataset its block does not list has no row, in the pool either.
        listed = {row["dataset"] for row in scored["per_dataset"]}
        for row in block_of(comparator, "per_dataset", f"{population['file']} comparator_concordance", list):
            dataset = row.get("dataset") if isinstance(row, dict) else None
            if not isinstance(dataset, str) or not dataset:
                fail(f"{population['file']} comparator_concordance.per_dataset holds an entry that names no dataset")
            if dataset not in listed:
                reader.note(
                    f"note: {population['file']} scores {comparator_tool['name']} on {dataset}, which the release record's "
                    "population_concordance does not list; the figure omits it"
                )
    pooled = [pooled_row(reader, pool, len(clinvar) + len(rows))]

    structural = []
    for cell, assembly in SV_CELLS:
        structural.append(dataset_cell(reader, cell, f"{pooled_files(release, concordance, cell)} {assembly}"))
    count_word = NUMBER_WORDS.get(len(clinvar) + len(rows), len(clinvar) + len(rows))
    labels = [label.format(count=count_word) for label in DATASET_GROUPS]
    return [dict(zip(("label", "rows"), group)) for group in zip(labels, (clinvar, rows, pooled, structural))]


def panel_wall_time(reader: Reader) -> list[dict]:
    """Panel (b): one group per dataset and one row per architecture, each tool the registry admits to a
    speed panel as a marker. ClinVar reads whole-file medians, a population dataset the sum of its
    per-chromosome medians; two medians, or two sums, of one tool, dataset and architecture are refused by
    name."""
    release, population = reader.chosen["release"], reader.chosen["population"]
    tools = [tool for tool in reader.tools if tool["speed"]]
    wall = block_of(population["data"], "wall_time", population["file"])
    groups = []
    for cell, name, assembly in CLINVAR_CELLS:
        rows = []
        for arch, arch_label in ARCHES:
            values = []
            for tool in tools:
                if tool["role"] == "released":
                    record, engine = release, None
                    matches = [
                        a for a in block_of(release["data"], "per_cell_aggregates", release["file"], list)
                        if isinstance(a, dict) and a.get("suite_id") == cell and a.get("arch") == arch and "target_cpu" not in a
                    ]
                    if not matches:
                        fail(f"{release['file']} has no per_cell_aggregates entry for {cell} on {arch}")
                else:
                    record = reader.chosen["walltime"].get(tool["id"])
                    aggregates = record["data"].get("per_cell_aggregates") if record is not None else None
                    matches = [
                        a for a in aggregates if isinstance(a, dict) and a.get("engine") in tool["keys"]
                        and a.get("suite_id") == cell and a.get("arch") == arch
                    ] if isinstance(aggregates, list) else []
                    if not matches:
                        reader.note(f"note: no record holds {tool['name']}'s wall time on {name} {assembly} ({arch_label}); the figure omits that mark")
                        continue
                    engine = matches[0]["engine"]
                if len(matches) > 1:
                    whose = "" if engine is None else f" of {engine}"
                    fail(f"{record['file']} has {len(matches)} per_cell_aggregates entries{whose} for {cell} on {arch}; the figure reads one")
                # The released tool's median is read from the entry that names no CPU target, and its key names that
                # condition: a record may hold further entries of the cell and architecture that name one.
                selector = f"suite_id={cell},arch={arch},no target_cpu" if engine is None else f"engine={engine},suite_id={cell},arch={arch}"
                key = f"per_cell_aggregates[{selector}].median_sec"
                values.append(reader.value(tool, "marker", seconds(matches[0].get("median_sec"), 2, f"{record['file']} {key}"), record, key))
                reader.run_count(tool, matches[0].get("n_clones"), record, f"per_cell_aggregates[{selector}].n_clones")
            rows.append({"label": arch_label, "values": values})
        groups.append({"label": f"{name} {assembly}", "rows": rows})
    for dataset in block_of(population["data"], "datasets", population["file"], list):
        if not isinstance(dataset, dict) or not all(dataset.get(k) for k in ("dataset", "label", "assembly")):
            fail(f"{population['file']} datasets holds an entry without dataset, label and assembly")
        label = f"{short_label(dataset['label'])} {dataset['assembly']}"
        rows = []
        for arch, arch_label in ARCHES:
            values = []
            for tool in tools:
                engine = next((key for key in tool["keys"] if isinstance(wall.get(key), dict)), None)
                sums = wall[engine].get("per_dataset") if engine is not None else None
                matches = [
                    s for s in sums if isinstance(s, dict) and s.get("dataset") == dataset["dataset"] and s.get("arch") == arch
                ] if isinstance(sums, list) else []
                key = f"wall_time.{engine}.per_dataset[dataset={dataset['dataset']},arch={arch}].sum_of_cell_medians_sec"
                if len(matches) > 1:
                    fail(f"{population['file']} has {len(matches)} wall_time sums of {tool['name']} for {dataset['dataset']} on {arch}; the figure reads one")
                if not matches:
                    if tool["role"] == "released":
                        fail(f"{population['file']} has no wall_time sum of {tool['name']} for {dataset['dataset']} on {arch}")
                    reader.note(f"note: no record holds {tool['name']}'s wall time on {label} ({arch_label}); the figure omits that mark")
                    continue
                values.append(reader.value(
                    tool, "marker", seconds(matches[0].get("sum_of_cell_medians_sec"), 0, f"{population['file']} {key}"), population, key,
                ))
                # A sum's entry states the runs of each of its cells, or its engine's block states them for every cell.
                held, where = (matches[0], key.rsplit(".", 1)[0]) if "n_clones_per_cell" in matches[0] else (wall[engine], f"wall_time.{engine}")
                reader.run_count(tool, held.get("n_clones_per_cell"), population, f"{where}.n_clones_per_cell")
            rows.append({"label": arch_label, "values": values})
        groups.append({"label": label, "rows": rows})
    return groups


def plugin_entries(plugins: dict) -> dict:
    """Every entry of the plugin record's entry lists as {(list, plugin, assembly, dataset or corpus, engine,
    projection): [(entry, where)]}, in record order."""
    index: dict = {}
    for name in (*(group[1] for group in PLUGIN_GROUPS), POPULATION_LIST, CUSTOM_GROUP[0]):
        entries = plugins["data"].get(name)
        if entries is None and name in (POPULATION_LIST, CUSTOM_GROUP[0]):
            continue
        for position, entry in enumerate(block_of(plugins["data"], name, plugins["file"], list)):
            where = f"{plugins['file']} {name}[{position}]"
            if not isinstance(entry, dict):
                fail(f"{where} is not an object")
            target = next((entry[k] for k in ENTRY_TARGETS if entry.get(k)), None)
            key = (name, entry.get("plugin"), entry.get("assembly"), target, entry.get("engine"), entry.get("projection"))
            if not all(isinstance(part, str) and part for part in key):
                fail(f"{where} carries no plugin, assembly, dataset or corpus, engine and projection")
            index.setdefault(key, []).append((entry, where))
    return index


def summed_counts(sets: list, keys: tuple[str, ...]) -> list[int] | None:
    """The sums of `keys` (matched, reference, engine) over an entry key's sets, or None when the sums hold
    no tuple on either side."""
    matched, perl, engine = (sum(count(entry, key, where) for entry, where in sets) for key in keys)
    return [matched, perl, engine] if perl + engine else None


def entry_selector(row_key: tuple, engine: str, projection: str, sets: list) -> str:
    """The values file's name for an entry key's sets: the entries of one engine and projection on a row."""
    name, plugin, assembly, target = row_key
    target_name = next(k for k in ENTRY_TARGETS if sets[0][0].get(k))
    return f"{name}[plugin={plugin},assembly={assembly},{target_name}={target},engine={engine},projection={projection}]"


def entry_key(parts: list[tuple[tuple, list]], engine: str, projection: str, counts: tuple[str, ...]) -> str:
    """The values file's key for a figure summed over the sets of one or more entry keys, each given with its
    sets: it names every key whose entries the figure sums."""
    number = sum(len(sets) for _, sets in parts)
    summed = f", summed over {number} entries" if number > 1 else ""
    return f"{' and '.join(entry_selector(key, engine, projection, sets) for key, sets in parts)}: {', '.join(counts)}{summed}"


def dataset_families(labels: list[str]) -> list[str]:
    """The plain names of the datasets `labels` name, in their order: labels that open with one word name one
    dataset in several releases, which goes by the leading words they share. Each label holds a word, which
    `population_datasets` holds the record to."""
    families: dict[str, list[str]] = {}
    for label in labels:
        words = label.split()
        shared = families.setdefault(words[0], words)
        families[words[0]] = shared[: next((at for at, pair in enumerate(zip(shared, words)) if pair[0] != pair[1]), min(len(shared), len(words)))]
    return [" ".join(words) for words in families.values()]


def population_datasets(reader: Reader) -> dict[str, tuple[str, str]]:
    """The whole-genome population datasets as {dataset: (label, assembly)}, in the population record's order,
    each label as panels (a) and (b) spell it. A label that holds no word so spelled is refused by name: the
    plugin panel's group and texts name a dataset by its label's words."""
    population = reader.chosen["population"]
    held = {}
    for position, dataset in enumerate(block_of(population["data"], "datasets", population["file"], list)):
        if not isinstance(dataset, dict) or not all(isinstance(dataset.get(k), str) and dataset[k] for k in ("dataset", "label", "assembly")):
            fail(f"{population['file']} datasets holds an entry without dataset, label and assembly")
        label = short_label(dataset["label"])
        if not label.split():
            fail(
                f"{population['file']} datasets[{position}].label holds no word: {dataset['label']!r}; the plugin panel names a population "
                "dataset by the words of its label"
            )
        held[dataset["dataset"]] = (label, dataset["assembly"])
    return held


def population_statements(reader: Reader, index: dict, datasets: dict) -> dict[tuple, str]:
    """What the plugin record holds on the population datasets, held to the population record: returns its
    statements that a plugin has no count on a dataset as {(plugin, assembly, dataset): reason}. Refused by
    name: an entry or a statement on a dataset the population record does not list or of another assembly
    than its dataset's; a statement without a reason; and a statement beside an entry."""
    plugins, population = reader.chosen["plugins"], reader.chosen["population"]

    def held_to_the_records(where: str, assembly, dataset) -> None:
        if dataset not in datasets:
            fail(f"{where} is on {dataset}, which {population['file']} does not list among its datasets")
        if datasets[dataset][1] != assembly:
            fail(f"{where} states {assembly} for {dataset}, which {population['file']} lists on {datasets[dataset][1]}")

    scored = set()
    for key, sets in index.items():
        if key[0] == POPULATION_LIST:
            held_to_the_records(sets[0][1], *key[2:4])
            scored.add(key[1:4])
    stated: dict[tuple, str] = {}
    if plugins["data"].get(POPULATION_ABSENT) is None:
        return stated
    for position, item in enumerate(block_of(plugins["data"], POPULATION_ABSENT, plugins["file"], list)):
        where = f"{plugins['file']} {POPULATION_ABSENT}[{position}]"
        named = tuple(item.get(key) for key in ("plugin", "assembly", "dataset")) if isinstance(item, dict) else ()
        if len(named) != 3 or not all(isinstance(part, str) and part for part in named):
            fail(f"{where} names no plugin, assembly and dataset")
        reason = item.get("reason")
        if not isinstance(reason, str) or not reason:
            fail(f"{where} states no reason")
        held_to_the_records(where, *named[1:])
        if named in scored:
            fail(f"{where} states that {named[0]} has no count on {named[2]}, and {POPULATION_LIST} holds an entry of it")
        stated[named] = reason
    return stated


def plugin_row(reader: Reader, index: dict, absent: dict, engines: dict, row_key: tuple, label: str, gate: bool,
               spans: tuple[list[tuple[str, str]], dict] | None = None) -> dict:
    """One row of panel (c): the released tool's per-feature raw bar and adjusted marker; a comparator's
    per-feature raw bar, with its adjusted marker where every entry of the row holds its adjusted counts and
    none states that the adjusted figure is not determined, or the record's reason where it lists the plugin
    as absent for the comparator. A gate row without the released tool's entry, or one that neither scores a comparator of the record per
    feature nor lists the plugin as absent for it, is refused, or noted in a draft. No per-allele figure is
    read: none is drawn.

    A row scored on real data is given `spans`: the population datasets of its assembly as (dataset, label), and
    the record's statements that a plugin has no count on a dataset. Each tool's figure on it is summed over
    the row's own entries and the entries the record holds for that tool on those datasets, and for no tool is
    a dataset read that the record does not hold for it. The released tool's row that lacks a dataset's entry
    without the record's statement is refused, or in a draft noted and drawn from what the record holds."""
    plugins = reader.chosen["plugins"]
    plugin = row_key[1].split(" (")[0]
    values = []
    sets: dict[str, int] = {}
    set_keys: dict[str, str] = {}
    for tool in concordance_tools(reader):
        engine = engines.get(tool["id"])
        feature = index.get((*row_key, engine, "feature")) if engine else None
        parts = [(row_key, feature)] if feature else []
        lacking = []
        for dataset, name in spans[0] if spans and feature else ():
            key = (POPULATION_LIST, *row_key[1:3], dataset)
            held = index.get((*key, engine, "feature"))
            if held:
                parts.append((key, held))
            elif (*row_key[1:3], dataset) not in spans[1]:
                lacking.append(name)
        summed = [entry for _, held in parts for entry in held]
        if tool["role"] == "released":
            raw = summed_counts(summed, RAW_COUNTS) if feature else None
            if raw is None:
                if gate and not reader.draft:
                    fail(f"{plugins['file']} holds no per-feature tuple of {tool['name']} on {label}; every row the record lists under rows states the released tool's figure")
                if gate:
                    reader.note(f"note: {plugins['file']} holds no per-feature tuple of {tool['name']} on {label}; the figure omits that bar")
                continue
            if lacking:
                missing = (f"{plugins['file']} holds neither a per-feature entry of {tool['name']} on {label} for {and_list(lacking)} in "
                           f"{POPULATION_LIST} nor a statement in {POPULATION_ABSENT} that the plugin has no count there")
                if not reader.draft:
                    fail(f"{missing}; a row scored on real data covers every population dataset of its assembly")
                reader.note(f"note: {missing}; the figure draws that row from the datasets the record holds")
            where = feature[0][1]
            sets[tool["name"]] = len(feature)
            set_keys[tool["name"]] = f"{entry_selector(row_key, engine, 'feature', feature)}: {SETS_KEY}"
            values.append(reader.value(tool, "bar", f1(*raw, where), plugins, entry_key(parts, engine, "feature", RAW_COUNTS), raw))
            if spans is not None:
                reader.coverage.append((tool, label, row_key[2], [key[3] for key, _ in parts[1:]]))
            if all(all(key in entry for key in ADJUSTED_COUNTS) for entry, _ in summed):
                adjusted = summed_counts(summed, ADJUSTED_COUNTS)
                if adjusted is not None:
                    values.append(reader.value(
                        tool, "marker", reader.adjusted(tool, adjusted, where), plugins,
                        entry_key(parts, engine, "feature", ADJUSTED_COUNTS), adjusted,
                    ))
                else:
                    reader.emptied_marker(tool, label, plugins, entry_key(parts, engine, "feature", ADJUSTED_COUNTS))
            else:
                reader.note(f"note: {plugins['file']} holds no adjusted counts of {tool['name']} on {label}; the figure omits that marker")
            continue
        raw = summed_counts(summed, RAW_COUNTS) if feature else None
        if raw is not None:
            where = feature[0][1]
            sets[tool["name"]] = len(feature)
            set_keys[tool["name"]] = f"{entry_selector(row_key, engine, 'feature', feature)}: {SETS_KEY}"
            values.append(reader.value(tool, "bar", f1(*raw, where), plugins, entry_key(parts, engine, "feature", RAW_COUNTS), raw))
            if spans is not None:
                reader.coverage.append((tool, label, row_key[2], [key[3] for key, _ in parts[1:]]))
            stated = [(entry[UNDETERMINED], at) for entry, at in summed if entry.get(UNDETERMINED) is not None]
            if any(count(block, UNDETERMINED_COUNT, f"{at}.{UNDETERMINED}") for block, at in stated):
                # A declared class that is not applied would set values of this tool aside: no adjusted F1 stands.
                reader.undetermined_marker(tool, label, plugins, entry_key(parts, engine, "feature", (f"{UNDETERMINED}.{UNDETERMINED_COUNT}",)))
            else:
                held = all(all(key in entry for key in ADJUSTED_COUNTS) for entry, _ in summed)
                adjusted = summed_counts(summed, ADJUSTED_COUNTS) if held else None
                key = entry_key(parts, engine, "feature", ADJUSTED_COUNTS)
                if held and adjusted is None:
                    reader.emptied_marker(tool, label, plugins, key)
                else:
                    reader.comparator_marker(values, tool, adjusted, where, plugins, key, label)
        if gate and not feature and (plugin, engine) in absent:
            values.append(reader.value(tool, "text", absent[(plugin, engine)], plugins, f"absent[plugin={plugin},engine={engine}].reason"))
        elif gate and engine and not feature:
            missing = f"{plugins['file']} neither scores {tool['name']} per feature on {label} nor lists the plugin as absent for it"
            if not reader.draft:
                fail(f"{missing}; every row the record lists under rows shows each tool the record scores")
            reader.note(f"note: {missing}; the figure draws no mark of {tool['name']} on that row")
    if len(set(sets.values())) > 1:
        # The tools' figures on this row are summed over different numbers of the sets of its own key: the texts
        # say so, and each number is kept with the entries it counts. A population dataset is no such set.
        reader.uneven.append((label, sets, set_keys))
    return {"label": label, "values": values}


def panel_plugins(reader: Reader) -> list[dict]:
    """Panel (c): the gate rows of the plugin record's `rows` in record order, grouped by the kind of data each
    is scored on, then the custom-annotation sources when the record holds them."""
    plugins = reader.chosen["plugins"]
    index = plugin_entries(plugins)
    engines = {}
    for block in block_of(plugins["data"], "engines", plugins["file"], list):
        tool = tool_of(reader.tools, block.get("id") if isinstance(block, dict) else None)
        if tool is not None:
            engines[tool["id"]] = block["id"]
    absent = {}
    for position, item in enumerate(block_of(plugins["data"], "absent", plugins["file"], list)):
        reason = item.get("reason") if isinstance(item, dict) else None
        if not isinstance(reason, str) or not reason:
            fail(f"{plugins['file']} absent[{position}] states no reason")
        named = (item.get("plugin"), item.get("engine"))
        if not all(isinstance(part, str) and part for part in named):
            fail(f"{plugins['file']} absent[{position}] names no plugin and engine")
        absent[named] = reason
    for engine in dict.fromkeys(engine for _, engine in absent):
        reasons = list(dict.fromkeys(reason for (_, named), reason in absent.items() if named == engine))
        if len(reasons) > 1:
            fail(
                f"{plugins['file']} absent states {len(reasons)} reasons for {engine} ({and_list(reasons)}); the figure draws one plain "
                "mark for a tool's absent plugins and names one reason for it"
            )
    grouped: dict[str, list[dict]] = {kind: [] for kind, *_ in PLUGIN_GROUPS}
    lists = {kind: name for kind, name, *_ in PLUGIN_GROUPS}
    listed = []
    for position, row in enumerate(block_of(plugins["data"], "rows", plugins["file"], list)):
        where = f"{plugins['file']} rows[{position}]"
        kind = row.get("scored_on") if isinstance(row, dict) else None
        if not isinstance(kind, str) or kind not in grouped:
            fail(f"{where} is scored on {kind!r}; the figure groups {', '.join(grouped)}")
        parts = tuple(row.get(key) for key in ("plugin", "assembly", "dataset_or_corpus"))
        if not all(isinstance(part, str) and part for part in parts):
            fail(f"{where} carries no plugin, assembly and dataset_or_corpus")
        listed.append((kind, parts))
    datasets = reader.population = population_datasets(reader)
    stated = population_statements(reader, index, datasets)
    for kind, parts in listed:
        label = f"{parts[0]} {parts[1]}" + (f" {SYNTHETIC_TAG}" if kind == SYNTHETIC_KIND else "")
        spans = ([(dataset, name) for dataset, (name, assembly) in datasets.items() if assembly == parts[1]], stated) if kind == REAL_KIND else None
        grouped[kind].append(plugin_row(reader, index, absent, engines, (lists[kind], *parts), label, True, spans))
    reader.row_names = list(dict.fromkeys(row["plugin"] for row in plugins["data"]["rows"]))
    labels = {kind: label for kind, _, label, _ in PLUGIN_GROUPS}
    summed = {dataset for tool, _, _, held in reader.coverage if tool["role"] == "released" for dataset in held}
    reader.real_label = labels[REAL_KIND] = and_list([labels[REAL_KIND], *dataset_families([datasets[d][0] for d in datasets if d in summed])])
    groups = [{"label": labels[kind], "rows": grouped[kind]} for kind, *_ in PLUGIN_GROUPS if grouped[kind]]
    custom = list(dict.fromkeys(key[:4] for key in index if key[0] == CUSTOM_GROUP[0]))
    if custom:
        rows = [plugin_row(reader, index, absent, engines, key, f"{key[1]} {key[2]}", False) for key in custom]
        groups.append({"label": CUSTOM_GROUP[1], "rows": rows})
    return groups


def support_blocks(reader: Reader) -> tuple[list[tuple[dict, dict, str]], list[dict]]:
    """The plugin record's `support` blocks of the released tool and the registry's comparators as (tool,
    block, where), in registry order, and the tools among them the record holds no block of. The reference has
    no block: the counts are of the plugins the released tool supports."""
    plugins = reader.chosen["plugins"]
    blocks = block_of(plugins["data"], "support", plugins["file"], list)
    found, missing = [], []
    for tool in reader.tools:
        if tool["role"] == "reference":
            continue
        held = [block for block in blocks if isinstance(block, dict) and block.get("engine") in tool["keys"]]
        found += [(tool, block, f"{plugins['file']} support[engine={block['engine']}]") for block in held]
        if not held:
            missing.append(tool)
    return found, missing


def panel_support(reader: Reader) -> list[dict]:
    """Panel (d): per tool, by its name, the plugins it reads as a bar and `<read> of <supported>` as its
    label. The panel is built when the record holds a count for the released tool and every comparator of the
    registry."""
    plugins = reader.chosen["plugins"]
    rows = []
    for tool, block, where in support_blocks(reader)[0]:
        read, of = count(block, "read", where), count(block, "of", where)
        selector = f"support[engine={block['engine']}]"
        rows.append({"label": tool["name"], "values": [
            reader.value(tool, "bar", str(read), plugins, f"{selector}.read"),
            reader.value(tool, "text", f"{read} of {of}", plugins, f"{selector}.read, of"),
        ]})
    return [{"label": "", "rows": rows}]


def genome_statements(entry: dict) -> list[str]:
    """The clauses of an `absent` entry of the genomes block, one for each of its statements: what the
    comparator's own documentation says, as a caption prints it. Empty where the entry holds no statement, or
    one that states no such clause."""
    statements = entry["stated_by"].get("statements")
    if not isinstance(statements, list) or not statements:
        return []
    says = [item.get("says") if isinstance(item, dict) else None for item in statements]
    return says if all(isinstance(clause, str) and clause for clause in says) else []


def panel_genomes(reader: Reader) -> list[dict]:
    """Panel (e): the released tool's row, one bar for the block's total of caches that read exact with
    `<exact> of <checked>` as its label, and one row for each comparator of the registry that the block's
    `absent` lists, which holds the entry's reason and no bar: the drawing sets a cross there. A comparator
    with a mark in an earlier panel that the block's `absent` does not list is refused by name, and so is an
    entry without a statement its sentence is built from; a draft notes either and leaves the row out, and
    notes an entry stated at another version than the other records state. A block without `absent` gives the
    released tool's row alone."""
    release = reader.chosen["release"]
    block, entries = reader.genomes
    where = f"{release['file']} {GENOMES_BLOCK}"
    totals = block_of(block, "totals", where)
    exact, checked = count(totals, "exact", f"{where}.totals"), count(totals, "checked", f"{where}.totals")
    if not 0 < checked >= exact:
        fail(f"{where}.totals counts {exact} exact of {checked} checked; the bar is the exact caches among one or more checked")
    rows = [{"label": reader.released["name"], "values": [
        reader.value(reader.released, "bar", str(exact), release, f"{GENOMES_BLOCK}.totals.exact"),
        reader.value(reader.released, "text", f"{exact:,} of {checked:,}", release, f"{GENOMES_BLOCK}.totals.exact, checked"),
    ]}]
    for tool in reader.tools if entries is not None else []:
        if tool["role"] != "comparator":
            continue
        entry = entries.get(tool["id"])
        if entry is None:
            if tool["id"] in reader.drawn:
                lacking = f"{where}.absent holds no entry for {tool['name']}, which has a mark in another panel"
                if not reader.draft:
                    fail(f"{lacking}; each comparator the figure draws has a row in the panel of the caches, and --draft draws a marked draft without it")
                reader.note(f"note: {lacking}; the draft leaves its row out")
            continue
        if not genome_statements(entry):
            unstated = f"{entry['where']}.stated_by holds no statement that gives what {tool['name']}'s documentation says"
            if not reader.draft:
                fail(f"{unstated}; the caption's sentence on the tool is built from its statements, and --draft draws a marked draft without its row")
            reader.note(f"note: {unstated}; the draft leaves its row out")
            continue
        version = version_text(tool, entry["stated_by"]["version"])
        others = [named for named in reader.versions[tool["id"]] if named != version]
        if others:
            # Reached in a draft alone: outside one, a tool at two versions is refused when the versions are read.
            reader.note(
                f"note: {entry['where']} is stated at {version} and the other records name {tool['name']} at {and_list(others)}; the draft "
                "draws its row and names each version"
            )
        rows.append({"label": tool["name"], "values": [
            reader.value(tool, "text", entry["reason"], release, f"{GENOMES_BLOCK}.absent[engine={entry['engine']}].reason"),
        ]})
    return [{"label": "", "rows": rows}]


def crossed_rows(panel: dict) -> list[tuple[dict, dict]]:
    """The rows of a panel of counts that hold a text and no bar, each with that text's value: the rows the
    drawing sets a cross on."""
    rows = [row for group in panel["groups"] for row in group["rows"]]
    return [
        (row, value) for row in rows if not any(held["kind"] == "bar" for held in row["values"])
        for value in row["values"] if value["kind"] == "text"
    ]


def genomes_sentences(reader: Reader, panel: dict) -> dict:
    """The sentences the caption and the alternative text take for panel (e), built once, every number, label
    and name from the block and each kept among the values file's sources: `counted`, what the bar counts,
    against the reference the block names; `sets`, the counts of the block's sets under the block's labels;
    `not_exact`, that the caches that are not exact are named on the supported-genomes page, or None where the
    block counts none; and `absent`, one sentence for each comparator with a row, from its name, the version
    its entry states and the entry's clauses."""
    if reader.genome_sentences is not None:
        return reader.genome_sentences
    release = reader.chosen["release"]
    block, entries = reader.genomes
    where = f"{release['file']} {GENOMES_BLOCK}"
    reference = reader.source(
        f"{reader.reference['name']} {image_release(block['reference'], f'{where}.reference')}", release, f"{GENOMES_BLOCK}.reference.image",
    )
    sets = []
    for position, item in enumerate(block_of(block, "sets", where, list)):
        at = f"{where}.sets[{position}]"
        label = item.get("label") if isinstance(item, dict) else None
        if not isinstance(label, str) or not label:
            fail(f"{at} carries no label")
        counts = f"{label} {count(item, 'exact', at):,} of {count(item, 'checked', at):,}"
        sets.append(reader.source(counts, release, GENOMES_SET_KEY.format(label=label)))
    if not sets:
        fail(f"{where}.sets holds no set; the texts give the counts of each")
    absent = []
    for _, value in crossed_rows(panel):
        tool = next(tool for tool in reader.tools if tool["id"] == value["tool"])
        entry = entries[tool["id"]]
        named = f"{tool['name']} {version_text(tool, entry['stated_by']['version'])}"
        says = and_list([f"that {clause}" for clause in genome_statements(entry)])
        absent.append(reader.source(GENOMES_ABSENT.format(tool=named, says=says), release, GENOMES_ABSENT_KEY.format(engine=entry["engine"])))
    reader.genome_sentences = {
        "counted": GENOMES_COUNTED.format(released=reader.released["name"], reference=reference),
        "sets": GENOMES_SETS.format(sets=and_list(sets)),
        "not_exact": GENOMES_NOT_EXACT if block["totals"]["exact"] < block["totals"]["checked"] else None,
        "absent": absent,
    }
    return reader.genome_sentences


def median_lines(panel: dict) -> tuple[list[str], bool]:
    """Of a wall-time panel, from the keys of its values: the datasets whose lines are medians, each named
    without its assembly, and whether the panel also holds lines that are sums of medians."""
    medians, sums = [], False
    for group in panel["groups"]:
        keys = [value["key"] for row in group["rows"] for value in row["values"]]
        if any(key.endswith(".median_sec") for key in keys) and group["label"].rsplit(" ", 1)[0] not in medians:
            medians.append(group["label"].rsplit(" ", 1)[0])
        sums = sums or any(key.endswith(".sum_of_cell_medians_sec") for key in keys)
    return medians, sums


def run_groups(reader: Reader) -> list[tuple[str, list[dict]]]:
    """The tools of the wall-time panel by the runs a cell of theirs holds, fewest first: each group as the
    number its tools' records state, a span where a tool's cells differ, and those tools. Empty where the panel
    holds no value, and where a record states no number for a cell of one."""
    spans: dict[tuple[int, int], list[dict]] = {}
    for tool in reader.tools:
        stated = [runs for runs, _, _ in reader.runs.get(tool["id"], [])]
        if None in stated:
            return []
        if stated:
            spans.setdefault((min(stated), max(stated)), []).append(tool)
    return [(str(low) if low == high else f"{low} to {high}", spans[low, high]) for low, high in sorted(spans)]


def runs_clause(reader: Reader) -> str:
    """The clause of the caption and the alternative text on the runs a wall-time cell holds, which names a cell
    as what it is, a dataset or a separately timed part of one: one number where every tool's cells hold as
    many, otherwise the fewest with its tools, against each other number with its tools. Empty without the
    groups of `run_groups`."""
    groups = run_groups(reader)
    if not groups:
        return ""
    (fewest, tools), others = groups[0], groups[1:]
    timed = f"{RUNS_SUBJECT} is timed in {fewest} run{'' if fewest == '1' else 's'}"
    if not others:
        return timed
    against = and_list([f"{runs} for {and_list([tool['name'] for tool in named])}" for runs, named in others])
    return f"{timed} for {and_list([tool['name'] for tool in tools])}, against {against}"


def methodology_clause(reader: Reader) -> str | None:
    """What a dataset's sum of cell medians is of, in the population record's own words: the first clause of
    its `methodology_note` that speaks of a sum of medians, kept with its key. Without such a clause the
    caption carries none, and a note is printed."""
    population = reader.chosen["population"]
    note = population["data"].get("methodology_note")
    clauses = re.split(r"[.;]\s+", note.strip().rstrip(".")) if isinstance(note, str) else []
    found = next((clause.strip() for clause in clauses if "sum of" in clause and "median" in clause), None)
    if found is None or not found.isprintable() or '"' in found:
        reader.note(f"note: {population['file']} carries no methodology note that says what a dataset's sum is of; the caption omits it")
        return None
    return reader.source(found, population, "methodology_note")


def adjusted_note(reader: Reader) -> str | None:
    """What a comparator's adjusted figure applies, in the population record's own words: the first sentence of
    its `comparator_concordance.note` that speaks of the adjusted figure, kept with its key. Without such a
    sentence the caption carries none: no sentence on it stands that no record holds."""
    population = reader.chosen["population"]
    block = population["data"].get("comparator_concordance")
    note = block.get("note") if isinstance(block, dict) else None
    sentences = re.split(r"(?<=\.)\s+", note.strip()) if isinstance(note, str) else []
    found = next((sentence for sentence in sentences if "adjusted" in sentence), None)
    if found is None or not found.isprintable() or '"' in found:
        return None
    return reader.source(found, population, "comparator_concordance.note")


def plugin_clauses(reader: Reader, panel: dict) -> list[str]:
    """The caption's clause for each group panel (c) holds, in panel order."""
    by_label = {label: (kind, clause) for kind, _, label, clause in PLUGIN_GROUPS}
    by_label[reader.real_label] = next((kind, clause) for kind, _, _, clause in PLUGIN_GROUPS if kind == REAL_KIND)
    by_label[CUSTOM_GROUP[1]] = (CUSTOM_GROUP[0], CUSTOM_GROUP[2])
    clauses = []
    for group in panel["groups"]:
        _, clause = by_label[group["label"]]
        names = list(dict.fromkeys(row["label"].rsplit(" ", 1)[0].split(" (")[0] for row in group["rows"]))
        clauses.append(clause.format(plugins=and_list(names), their="its" if len(names) == 1 else "their"))
    return clauses


def coverage_sentence(reader: Reader) -> str | None:
    """The one sentence on the datasets each tool's figure covers on the rows scored on real data, from the
    entries each row sums, the released tool first: the population datasets a tool's rows sum beside ClinVar,
    assembly by assembly and each as the population record labels it, then every row of the tool that covers
    fewer, with what it covers; a tool none of whose rows sums a population dataset covers ClinVar alone. The
    sentence is kept with its records and key; None where no row scored on real data carries a bar."""
    clauses = []
    for tool in sorted(concordance_tools(reader), key=lambda tool: tool["role"] != "released"):
        rows = [(label, assembly, held) for named, label, assembly, held in reader.coverage if named["id"] == tool["id"]]
        if not any(held for _, _, held in rows):
            clauses += [COVERAGE_ALONE.format(tool=tool["name"])] if rows else []
            continue
        found: dict[str, set[str]] = {}
        for _, assembly, held in rows:
            found.setdefault(assembly, set()).update(held)
        full = {assembly: [dataset for dataset in reader.population if dataset in held] for assembly, held in found.items() if held}

        def named(held: list[str]) -> str:
            return and_list([reader.population[dataset][0] for dataset in held])

        builds = "; ".join(f"{assembly}: {named(held)}" for assembly, held in full.items())
        fewer: dict[str, list[str]] = {}
        for label, assembly, held in rows:
            if held != full.get(assembly):
                fewer.setdefault(f"ClinVar and {named(held)}" if held else "ClinVar alone", []).append(label)
        listed = and_list([f"{and_list(labels)} ({covers})" for covers, labels in fewer.items()])
        clauses.append(COVERAGE_SUMMED.format(tool=tool["name"], builds=builds, fewer=COVERAGE_FEWER.format(rows=listed) if fewer else ""))
    if not clauses:
        return None
    whose = {"file": f"{reader.chosen['plugins']['file']} and {reader.chosen['population']['file']}"}
    return reader.source(COVERAGE_SENTENCE.format(clauses="; ".join(clauses)), whose, COVERAGE_KEY)


def stand_in_sentences(reader: Reader) -> str:
    """The caption's sentence for each reference of the plugin record that ran a plugin's module with stand-in
    files, with the package of a derived image, each kept with the keys it is drawn from; a stand-in for a
    plugin without a clause is refused by name."""
    plugins = reader.chosen["plugins"]
    sentences = []
    for index, reference in enumerate(block_of(plugins["data"], "perl_reference", plugins["file"], list)):
        stand_ins = reference.get("stand_ins") if isinstance(reference, dict) else None
        if stand_ins is None:
            continue
        plugin = stand_ins.get("plugin") if isinstance(stand_ins, dict) else None
        if plugin not in STAND_IN_CLAUSES:
            fail(f"{plugins['file']} perl_reference carries stand-ins for {plugin!r}; the caption has a clause for {', '.join(STAND_IN_CLAUSES)}")
        sentence, key = STAND_IN_CLAUSES[plugin], f"perl_reference[{index}].stand_ins.plugin"
        derived = reference.get("derived_image")
        if derived is not None:
            package = derived.get("package") if isinstance(derived, dict) else None
            if not isinstance(package, str) or not package or not package.isprintable():
                fail(f"{plugins['file']} perl_reference's derived image for {plugin} names no package: {package!r}")
            sentence += f", on the {reader.reference['name']} image with {package.split('=')[0]} added"
            key += f", perl_reference[{index}].derived_image.package"
        if f" {sentence}." not in sentences:
            sentences.append(f" {reader.source(f'{sentence}.', plugins, key)}")
    return "".join(sentences)


def circle_order(reader: Reader, panel: str, rows: list[dict]) -> str:
    """The alternative text's clause on the order of the tools' circles on a row of the wall-time panel, whose
    id is `panel`, from the left, where the shortest time stands: the panel tells its tools apart by colour
    alone, which a text read without the image cannot. One order is named for every row, and several orders
    each with the number of rows that have it; rows of one circle each have no order to name. Equal times keep
    the tools' own order. Each number of rows the clause prints is kept among the values file's sources, with
    what it counts of the panel's rows and the records of the panel's values."""
    names = {tool["id"]: tool["name"] for tool in reader.tools}
    orders: dict[tuple[str, ...], int] = {}
    for row in rows:
        timed = sorted((value for value in row["values"] if value["kind"] == "marker"), key=lambda value: Decimal(value["value"]))
        if timed:
            order = tuple(value["tool"] for value in timed)
            orders[order] = orders.get(order, 0) + 1
    if all(len(order) < 2 for order in orders):
        return ""
    named = {order: and_list([names[tool] for tool in order]) for order in orders}
    if len(orders) == 1:
        return f"; from left to right, the circles of every dataset and architecture are those of {named[next(iter(orders))]}"
    whose = {"file": " and ".join(dict.fromkeys(file for row in rows for value in row["values"] for file in value["record"].split(" and ")))}
    total = sum(orders.values())
    listed = []
    for order, number in orders.items():
        reader.source(str(number), whose, PANEL_ORDER_KEY.format(panel=panel, tools=", ".join(f"tool={tool}" for tool in order)))
        if not listed:
            # The number of rows that hold a circle is kept once, where the clause first prints it.
            reader.source(str(total), whose, PANEL_TIMED_KEY.format(panel=panel))
        listed.append(f"{named[order]} on {number} of the {total}")
    return f"; from left to right, the circles of a dataset and architecture are those of {and_list(listed)}"


def exact_sentence(exact: int, markers: int, released: str) -> str:
    """The sentence on where an F1 of exactly 1 is drawn, with how many of the released tool's adjusted markers
    stand there: its noun agrees with the number of markers and its verb with the number that stand there."""
    return EXACT_SENTENCE.format(exact=exact, one="s" if exact == 1 else "", markers=markers, plural="" if markers == 1 else "s", released=released)


def alt_sentence(reader: Reader, panel: dict) -> str:
    """The alternative text's sentence for a panel, which says what the panel holds as it is drawn: the tools
    with a mark in it and the mark each has. The wall-time panel's also names the order of the tools' circles
    on a row, from the left, and the runs a cell holds. The plugin panel's also says what its groups of rows are
    scored on, in the caption's clauses, and, for each comparator, on how many
    rows it has a figure, on how many of those its marker where that is fewer, and on how many a cross; and, as
    the panel's result, on how many rows the released tool's adjusted F1 is exactly 1. Every count is of the
    panel's own values, and each number of rows the sentence prints is kept among the values file's sources."""
    rows = [row for group in panel["groups"] for row in group["rows"]]

    def holding(kind: str, tool: dict) -> list[dict]:
        return [row for row in rows if any(value["tool"] == tool["id"] and value["kind"] == kind for value in row["values"])]

    def holders(kind: str) -> list[str]:
        return [tool["name"] for tool in reader.tools if holding(kind, tool)]

    axis = PANEL_AXES[panel["id"]]
    if axis == "seconds":
        arches = list(dict.fromkeys(row["label"] for row in rows))
        circles = f", filled for {arches[0]} and open for {arches[1]}" if len(arches) == 2 else ""
        runs = runs_clause(reader)
        return (f"Wall time per dataset and architecture of {and_list(holders('marker'))}, a circle each{circles}{circle_order(reader, panel['id'], rows)}"
                f"{'; ' if runs else ''}{runs}.")
    if panel["id"] == "e":
        label = next(value["value"] for row in rows for value in row["values"] if value["tool"] == reader.released["id"] and value["kind"] == "text")
        drawn = [GENOMES_DRAWN.format(released=reader.released["name"], label=label)]
        drawn += [GENOMES_CROSSED.format(tool=row["label"], reason=value["value"]) for row, value in crossed_rows(panel)]
        sentences = genomes_sentences(reader, panel)
        return " ".join([f"{GENOMES_CLAUSE[0].upper()}{GENOMES_CLAUSE[1:]}: {'; '.join(drawn)}.", sentences["sets"], *sentences["absent"]])
    if axis == "count":
        reads = [value["value"] for row in rows for value in row["values"] if value["kind"] == "text"]
        return f"The plugins read by {and_list(holders('bar'))}, a bar each: {and_list(reads)}."
    barred, marked = holders("bar"), holders("marker")
    groups = f" ({'; '.join(plugin_clauses(reader, panel))})" if panel["id"] == "c" else ""
    sentence = f"Raw F1 per {F1_ROWS[panel['id']]}{groups} as a bar for {'each of ' if len(barred) > 1 else ''}{and_list(barred)}"
    if marked:
        whose = ("each tool's" if len(barred) > 1 else "its") if marked == barred else f"{and_list(marked)}'s"
        sentence += f", with {whose} adjusted F1 as a marker joined to the end of the bar by a line where the two differ"
    if panel["id"] != "c":
        return f"{sentence}."
    def counted(kind: str, tool: dict) -> int:
        """The number of the panel's rows that hold a value of `kind` of `tool`, kept among the values file's
        sources as the sentence prints it, with what it counts."""
        held = len(holding(kind, tool))
        reader.source(str(held), reader.chosen["plugins"], PANEL_ROWS_KEY.format(panel=panel["id"], tool=tool["id"], kind=kind))
        return held

    total: list[str] = []

    def of_rows() -> str:
        """The words that close a count of the panel's rows, with the number of its rows, which is kept among the
        values file's sources once, where the sentence first prints it."""
        if not total:
            total.append(reader.source(str(len(rows)), reader.chosen["plugins"], PANEL_TOTAL_KEY.format(panel=panel["id"])))
        return f"of the {total[0]} row{'' if len(rows) == 1 else 's'}"

    clauses = []
    for tool in reader.tools:
        scored, crossed = holding("bar", tool), holding("text", tool)
        if tool["role"] != "comparator" or not (scored or crossed):
            continue
        clause = f"{tool['name']} has a figure on {counted('bar', tool)} {of_rows()}"
        if len(holding("marker", tool)) != len(scored):
            # The clause takes its closing comma only before the clause on crosses.
            clause += f", its marker on {counted('marker', tool)} of them" + ("," if crossed else "")
        if crossed:
            reasons = dict.fromkeys(value["value"] for row in crossed for value in row["values"] if value["tool"] == tool["id"] and value["kind"] == "text")
            clause += f" and a cross ({and_list(list(reasons))}) on {counted('text', tool)}"
        clauses.append(clause)
    ours = [value for row in holding("marker", reader.released) for value in row["values"]
            if value["tool"] == reader.released["id"] and value["kind"] == "marker"]
    if ours:
        exact = sum(not discordance(value) for value in ours)
        reader.source(str(exact), reader.chosen["plugins"], PANEL_EXACT_KEY.format(panel=panel["id"], tool=reader.released["id"]))
        clauses.append(f"{reader.released['name']}'s adjusted F1 is exactly 1 on {exact} {of_rows()}")
    return f"{sentence}: {'; '.join(clauses)}." if clauses else f"{sentence}."


def texts(reader: Reader, panels: list[dict]) -> tuple[str, str, list[str]]:
    """The caption, the alternative text and the note lines drawn inside the image, every number and version
    from the records; the panel of counts is given its note, which is drawn under it. The caption states the concordance of the released tool only over records whose every
    adjusted F1 of that tool is exactly 1; over any other it says that the render is a draft of a state that is
    not shipped."""
    by_id = {panel["id"]: panel for panel in panels}
    released = reader.released["name"]
    comparators = [reader.label(tool) for tool in reader.tools if tool["role"] == "comparator" and tool["id"] in reader.drawn]
    lead = f"{reader.label(reader.released)} against {reader.label(reader.reference)}"
    if comparators:
        lead += f", with {and_list(comparators)} as comparator{'' if len(comparators) == 1 else 's'}"
    caption = [f"**{lead}.**"]

    def marked(panel_id: str) -> list[str]:
        """The comparators with an adjusted marker in a panel."""
        return [tool["name"] for tool in reader.tools if tool["role"] == "comparator" and any(
            v["tool"] == tool["id"] and v["kind"] == "marker" for g in by_id[panel_id]["groups"] for r in g["rows"] for v in r["values"]
        )]

    adjusted = list(dict.fromkeys(name for panel_id in ("a", "c") if panel_id in by_id for name in marked(panel_id)))
    marker = "a marker is its tool's adjusted F1" if marked("a") else f"the marker is {released}'s adjusted F1"
    pooled = by_id["a"]["groups"][2]["rows"][0]["label"]
    caption.append(
        f"(a) Raw F1 per dataset on a log(1 - F1) axis; {marker}, {JOIN_CLAUSE}. {pooled} is the pool "
        f"of the {pooled.split(' ', 1)[1]} whole-genome datasets, the ClinVar and the population datasets."
    )
    medians, sums = median_lines(by_id["b"])
    statistics = []
    if medians and sums:
        # The record's key names a sum of cell medians: a cell is a separately timed part of the dataset, and
        # which parts they are, the record's methodology note says.
        clause = methodology_clause(reader)
        gloss = f' (the record\'s methodology note: "{clause}")' if clause else ""
        statistics.append(f"a {and_list(medians)} line is the median of its runs, any other the sum of the medians of {SUMMED_PARTS}{gloss}")
    runs = runs_clause(reader)
    if runs:
        statistics.append(runs)
        # Each number the clause prints, in the clause's order, with the keys it is read under, record by record.
        for number, named in run_groups(reader):
            for tool in named:
                keys: dict[str, list[str]] = {}
                for _, file, key in reader.runs[tool["id"]]:
                    keys.setdefault(file, []).append(key)
                for file, read in keys.items():
                    reader.source(number, {"file": file}, "; ".join(read))
    caption.append(f"(b) Wall time per dataset and architecture{': ' if statistics else ''}{'; '.join(statistics)}.")
    notes = []
    covered = None
    if "c" in by_id:
        covered = coverage_sentence(reader)
        whose = "its tool's" if marked("c") else f"{released}'s"
        reading = ["The bar is each tool's per-feature raw F1", f"the marker {whose} adjusted F1", f"{JOIN_CLAUSE}, as in (a)"]
        crossed = next((v for g in by_id["c"]["groups"] for r in g["rows"] for v in r["values"] if v["kind"] == "text"), None)
        absent = ""
        if crossed is not None:
            # Every cross of a tool stands for one reason, which the record states per plugin: one source names where.
            engine = crossed["key"].rsplit("engine=", 1)[1].split("]", 1)[0]
            reason = reader.source(crossed["value"], reader.chosen["plugins"], ABSENT_KEY.format(engine=engine))
            absent = f"; a cross marks a plugin for which the comparator has no counterpart (\"{reason}\")"
        caption.append(
            f"(c) Raw F1 per plugin: {'; '.join(plugin_clauses(reader, by_id['c']))}.{f' {covered}' if covered else ''} "
            f"{', '.join(reading)}{absent}.{stand_in_sentences(reader)}"
        )
    if "d" in by_id:
        plugins = reader.chosen["plugins"]
        _, block, where = next(found for found in support_blocks(reader)[0] if found[0]["role"] == "released")
        supported = count(block, "of", where)
        listed = count(block_of(plugins["data"], "not_ported", plugins["file"]), "of", f"{plugins['file']} not_ported")
        reader.source(str(listed), plugins, "not_ported.of")
        reads = and_list([
            f"{row['label']} {value['value']}" for group in by_id["d"]["groups"] for row in group["rows"] for value in row["values"] if value["kind"] == "bar"
        ])
        caption.append(
            f"(d) The plugins each tool reads, of the {supported} that {released} supports: {reads}; plugins {released} does not port are "
            "not counted and are listed by reason in the README."
        )
        by_id["d"]["note"] = COUNTS_NOTE.format(
            supported=supported, plural="" if supported == 1 else "s", listed=listed, released=released, reference=reader.reference["name"],
        )
    elif reader.chosen["plugins"] is not None:
        caption.append(SUPPORT_LEFT_OUT.format(tools=and_list([tool["name"] for tool in support_blocks(reader)[1]])))
    if "e" in by_id:
        label = next(
            value["value"] for group in by_id["e"]["groups"] for row in group["rows"] for value in row["values"]
            if value["tool"] == reader.released["id"] and value["kind"] == "text"
        )
        sentences = genomes_sentences(reader, by_id["e"])
        following = [sentences["counted"], sentences["sets"], *([sentences["not_exact"]] if sentences["not_exact"] else []), *sentences["absent"]]
        caption.append(f"(e) {GENOMES_CLAUSE[0].upper()}{GENOMES_CLAUSE[1:]}: {label}. {' '.join(following)}")
    statements = []
    if "c" in by_id:
        plugins = reader.chosen["plugins"]
        held: dict[str, list[str]] = {}
        for name in reader.row_names:
            held.setdefault(name.split(" (")[0], []).append(name)
        doubled = [f"{plugin} is scored under {and_list(names)}" for plugin, names in held.items() if len(names) > 1]
        if doubled:
            rows = sum(len(group["rows"]) for group in by_id["c"]["groups"] if group["label"] != CUSTOM_GROUP[1])
            sentence = NAMES_SENTENCE.format(
                rows=rows, names=len(reader.row_names), plugins=len(held), plural="" if len(held) == 1 else "s", doubled=and_list(doubled),
            )
            statements.append(reader.source(sentence, plugins, "rows: plugin"))
        if reader.uneven:
            listed = [f"{label} ({', '.join(f'{name} {number}' for name, number in sets.items())})" for label, sets, _ in reader.uneven]
            sentence = UNEVEN_SENTENCE.format(count=len(listed), plural="" if len(listed) == 1 else "s", rows=and_list(listed))
            statements.append(reader.source(sentence, plugins, UNEVEN_KEY))
            # Each number of sets the sentence prints, in the sentence's order, with the entries it counts.
            for _, sets, keys in reader.uneven:
                for name, number in sets.items():
                    reader.source(str(number), plugins, keys[name])
    caption += statements
    caption.append(PLACE_SENTENCE)
    caption.append(AGREE_SENTENCE)
    # A draft is not drawn from one release's final records: what a record's note says of a comparator's adjusted
    # figure may not hold for the figures a draft draws, so a draft quotes none.
    stated = adjusted_note(reader) if adjusted and not reader.draft else None
    if stated:
        caption.append(ADJUSTED_NOTE.format(sentence=stated))
    if reader.unadjusted:
        caption.append(ADJUSTED_LEFT_OUT.format(tools=and_list(reader.unadjusted)))
    for name, rows in reader.undetermined.items():
        line = ADJUSTED_UNDETERMINED.format(tool=name, rows=and_list([label for label, _ in rows]))
        caption.append(reader.source(line, reader.chosen["plugins"], "; ".join(key for _, key in rows)))
    for name, rows in reader.emptied.items():
        line = ADJUSTED_EMPTY.format(tool=name, rows=and_list([label for label, _ in rows]))
        caption.append(reader.source(line, reader.chosen["plugins"], "; ".join(key for _, key in rows)))
    ours = [
        value for panel in panels if PANEL_AXES[panel["id"]] == "f1" for group in panel["groups"] for row in group["rows"]
        for value in row["values"] if value["kind"] == "marker" and value["tool"] == reader.released["id"]
    ]
    shipped = all(matched == reference == engine for matched, reference, engine in (value["counts"] for value in ours))
    if ours:
        exact = sum(not discordance(value) for value in ours)
        sentence = exact_sentence(exact, len(ours), released)
        whose = {"file": " and ".join(dict.fromkeys(value["record"] for value in ours))}
        statements.append(reader.source(sentence, whose, EXACT_KEY.format(released=released)))
        caption.append(sentence)
    caption.append(CONCORDANCE_SENTENCE if shipped else UNSHIPPED_SENTENCE.format(released=released))
    caption.append(RECORDS_SENTENCE)

    against = f"{and_list([reader.label(reader.released), *comparators])} against {reader.label(reader.reference)}"
    sentences = " ".join(f"({panel['id']}) {alt_sentence(reader, panel)}" for panel in panels)
    first = plugin_split(by_id)
    if "c" in by_id and 0 < first < len(by_id["c"]["groups"]):
        # Read without the image, the alternative text alone says that the plugin panel runs over both columns.
        clause = CONTINUES_SENTENCE.format(group=by_id["c"]["groups"][first]["label"])
        statements.append(reader.source(clause, reader.chosen["plugins"], CONTINUES_KEY))
    alt = f"{COUNT_WORDS[len(panels)]} comparing {against}. {' '.join([sentences, AGREE_SENTENCE, *([covered] if covered else []), *statements])}"
    caption_text = " ".join(caption)
    if reader.draft:
        caption_text, alt, notes = f"{DRAFT_LINE} {caption_text}", f"{DRAFT_LINE} {alt}", [DRAFT_LINE, *notes]
    return caption_text, alt, notes


def build_values(tools: list[dict], chosen: dict, draft: bool) -> tuple[dict, list[str]]:
    """The values document and the notes to print: every plotted number with its record and key, each tool's
    name and version, the caption, the alternative text, the note lines and the record keys of what the texts
    state beyond the plotted numbers."""
    reader = Reader(tools, chosen, draft)
    read_versions(reader)
    check_one_configuration(reader)
    panel_records = {"a": chosen["release"], "b": chosen["population"], "c": chosen["plugins"], "d": chosen["plugins"], "e": chosen["release"]}
    builders = {"a": panel_datasets, "b": panel_wall_time, "c": panel_plugins, "d": panel_support, "e": panel_genomes}
    panels = []
    for panel_id, builder in builders.items():
        record = panel_records[panel_id]
        # Panel (d) compares counts: it is built when every tool it would compare has one. Panel (e) is built
        # when the release record holds the block of the genomes checked.
        if record is None or (panel_id == "d" and support_blocks(reader)[1]) or (panel_id == "e" and reader.genomes is None):
            continue
        groups = builder(reader)
        # The panel of counts names the tools it draws, each a row of it.
        named = and_list([row["label"] for group in groups for row in group["rows"]]) if panel_id == "d" else ""
        panels.append({"id": panel_id, "title": PANEL_TITLES[panel_id].format(tools=named), "groups": groups})
    for tool in tools:
        if tool["id"] in reader.drawn and tool["id"] not in reader.versions:
            fail(f"no record states a version of {tool['name']}; every mark names its tool's version")
    caption, alt, notes = texts(reader, panels)
    document = {
        "image": IMAGE,
        "canvas": {"width_px": CANVAS[0], "height_px": CANVAS[1]},
        "alt_text": alt,
        "caption": caption,
        "notes": notes,
        "sources": reader.sources,
        "tools": [
            {"id": tool["id"], "name": tool["name"], "version": " and ".join(reader.versions[tool["id"]]),
             "colour": tool["fill"]}
            for tool in tools if tool["id"] in reader.drawn
        ],
        "records": [{"file": record["file"], "sha256": record["sha256"]} for record in chosen["read"]],
        "panels": panels,
    }
    return document, reader.printed


def values_bytes(document: dict) -> bytes:
    return (json.dumps(document, ensure_ascii=False, indent=2) + "\n").encode("utf-8")


def tool_label(tool: dict) -> str:
    return f"{tool['name']} {tool['version']}"


def lane_tools(document: dict, panel: dict) -> list[str]:
    """The ids of the tools with a mark in `panel`, in the order of the document's tools: a row's lanes from
    the top."""
    present = {value["tool"] for group in panel["groups"] for row in group["rows"] for value in row["values"]}
    return [tool["id"] for tool in document["tools"] if tool["id"] in present]


def axis_lines(document: dict, tools: list[dict], panel: dict) -> list[str]:
    """The lines of a panel's value-axis label. An F1 axis names the reference tool, whose version the legend
    states, and no kind of F1: its bars are raw figures and its markers adjusted ones, as the legend says; a
    wall-time axis says which lines are medians and which are sums of medians. An F1 axis of a document whose
    tools hold no reference is refused by name: the legend states a version only of a tool with a mark."""
    axis = PANEL_AXES[panel["id"]]
    if axis == "f1":
        listed = [row for row in tools if row["role"] == "reference"]
        reference = next((tool for tool in document["tools"] if any(row["id"] == tool["id"] for row in listed)), None)
        if reference is None:
            fail(
                f"the values file's tools hold no {and_list([row['name'] for row in listed])}, the registry's reference, which an F1 "
                "axis names and whose version the legend states: a tool stands among them with a mark of its own, and the reference's "
                "marks are its wall times, read from the population record's wall_time.<key>.per_dataset and a wall-time record's "
                f"per_cell_aggregates[engine=<key>] under its record keys ({', '.join(key for row in listed for key in row['keys'])}) "
                "where the registry admits it to the speed panel"
            )
        return [f"F1 against {reference['name']}", F1_AXIS_SCALE]
    if axis == "seconds":
        medians, sums = median_lines(panel)
        return [SECONDS_AXIS, SECONDS_STATISTICS.format(medians=and_list(medians))] if medians and sums else [SECONDS_AXIS]
    return [GENOMES_AXIS if panel["id"] == "e" else COUNT_AXIS]


def discordance(value: dict) -> Fraction:
    """An F1 value's 1 - F1 from its integer counts, exact where the printed figure is rounded: zero only
    when the three counts are equal."""
    matched, reference, engine = value["counts"]
    return 1 - Fraction(2 * matched, reference + engine)


def marker_shape(value: dict) -> str:
    """How a marker of an F1 row is drawn: `exact`, a triangle at the floor, for an F1 of exactly 1, at the
    end of a full bar too; `adjusted`, a diamond, for one below it."""
    return "adjusted" if discordance(value) else "exact"


def legend_entries(document: dict) -> list[dict]:
    """The legend, one entry per tool and per class of mark the panels draw: {mark, tool, label}; `tool` is
    None for a mark that is no one tool's. A tool's entry is a fill when a panel holds a bar of it, a dot when
    it has circles alone and a fill with a dot beside it when it has both, bars in an F1 or count panel and
    circles in the wall-time panel; it names the tool with its version and nothing else: its colour says which
    bar and which circle are its own."""
    panels = {panel["id"]: panel for panel in document["panels"]}
    entries = []
    for tool in document["tools"]:
        held = {(PANEL_AXES[p["id"]], v["kind"]) for p in document["panels"] for g in p["groups"] for r in g["rows"] for v in r["values"] if v["tool"] == tool["id"]}
        barred, circled = any(kind == "bar" for _, kind in held), ("seconds", "marker") in held
        entries.append({"mark": "fill_dot" if barred and circled else "fill" if barred else "dot", "tool": tool["id"], "label": tool_label(tool)})
    f1_values = [v for p in document["panels"] if PANEL_AXES[p["id"]] == "f1" for g in p["groups"] for r in g["rows"] for v in r["values"]]
    if any(v["kind"] == "bar" for v in f1_values):
        # What a bar of an F1 panel states, for every tool: the axis names no kind of F1.
        entries.append({"mark": "raw", "tool": None, "label": RAW_LEGEND})
    markers = [v for v in f1_values if v["kind"] == "marker"]
    shapes = {marker_shape(v) for v in markers}
    # The markers of one tool are named as that tool's and wear its colour; those of several are named for
    # what they state, in no tool's colour.
    marked = [tool for tool in document["tools"] if any(v["tool"] == tool["id"] for v in markers)]
    whose, lead = (marked[0]["id"], f"{marked[0]['name']} adjusted F1") if len(marked) == 1 else (None, "Adjusted F1")
    entries += [
        {"mark": mark, "tool": whose, "label": f"{lead}, {words}"}
        for mark, words in (("exact", "exactly 1"), ("adjusted", "below 1")) if mark in shapes
    ]
    for tool in document["tools"]:
        for reason in dict.fromkeys(v["value"] for v in f1_values if v["kind"] == "text" and v["tool"] == tool["id"]):
            entries.append({"mark": "absent", "tool": tool["id"], "label": f"{tool['name']}: {reason}"})
    if "b" in panels:
        labels = list(dict.fromkeys(row["label"] for group in panels["b"]["groups"] for row in group["rows"]))
        entries += [{"mark": mark, "tool": None, "label": label} for mark, label in zip(("filled", "open"), labels)]
    return entries


def legend_rows(entries: list[dict]) -> list[list[dict]]:
    """The legend's entries by kind, a row for each kind present: the tools, then the marks of an F1, what a
    bar states first and the two adjusted entries side by side, then the architectures. The dashed line from
    a bar's end to its marker has no entry: the caption says what it is."""
    rows = [[entry for entry in entries if entry["mark"] in kinds] for kinds in LEGEND_KINDS]
    return [row for row in rows if row]


def drawn_texts(document: dict, tools: list[dict]) -> list[str]:
    """Every string the drawing layer places, apart from the tick labels of the value axes."""
    out = [entry["label"] for entry in legend_entries(document)] + list(document["notes"])
    for panel in document["panels"]:
        out += [f"({panel['id']})", panel["title"], *axis_lines(document, tools, panel)] + ([panel["note"]] if panel.get("note") else [])
        for group in panel["groups"]:
            out += [group["label"]] if group["label"] else []
            for row in group["rows"]:
                out.append(row["label"])
                # A reason of an F1 panel is drawn as a plain mark, which the legend names; a count's text as itself.
                if PANEL_AXES[panel["id"]] != "f1":
                    out += [value["value"] for value in row["values"] if value["kind"] == "text"]
    return out


def png_chunks(data: bytes):
    """A PNG's chunks as (type, body); a file that is not a PNG is refused."""
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        fail("the image is not a PNG")
    position = 8
    while position + 8 <= len(data):
        (length,), kind = struct.unpack(">I", data[position:position + 4]), data[position + 4:position + 8]
        yield kind, data[position + 8:position + 8 + length]
        position += 12 + length


def png_size(data: bytes) -> tuple[int, int]:
    """A PNG's width and height in pixels, from its header chunk."""
    kind, body = next(png_chunks(data))
    if kind != b"IHDR":
        fail("the PNG does not open with its header chunk")
    return struct.unpack(">II", body[:8])


def png_text(data: bytes) -> dict[str, str]:
    """A PNG's text chunks as {keyword: text}."""
    out = {}
    for kind, body in png_chunks(data):
        if kind == b"tEXt":
            keyword, _, value = body.partition(b"\x00")
            out[keyword.decode("latin-1")] = value.decode("latin-1")
        elif kind == b"zTXt":
            keyword, _, rest = body.partition(b"\x00")
            out[keyword.decode("latin-1")] = zlib.decompress(rest[1:]).decode("latin-1")
        elif kind == b"iTXt":
            keyword, _, rest = body.partition(b"\x00")
            compressed, rest = rest[0], rest[2:]
            _, _, rest = rest.partition(b"\x00")
            _, _, value = rest.partition(b"\x00")
            out[keyword.decode("latin-1")] = (zlib.decompress(value) if compressed else value).decode("utf-8")
    return out


DPI = 150
INK, INK_SECONDARY = "#0b0b0b", "#52514e"
GROUND, MIN_CONTRAST = "#ffffff", 4.5
# Every text is set in one of these inks on the ground, which fills the canvas and every plot. The contrast of
# text is measured on the colours as set, a text's and its ground's, never on pixels softened at a glyph's edge.
TEXT_INKS = (INK, INK_SECONDARY)
GRID = "#D1D5DB"
# F1 is plotted as 1 - F1 on a logarithmic axis, reversed, so further right is higher F1: on a linear axis
# every bar near 1 is full length and identical. log(1 - F1) needs a clip at one end only, the floor, where
# an F1 of exactly 1 is drawn. The floor stands no nearer than this, and at the first decade under a value
# that lies past it.
DISCORD_FLOOR = 1e-7
# The worse end of the scale: 1 - F1 = 1 is F1 = 0, where the tick that reads 0 stands.
CONCORDANCE_ANCHOR = 1.0
# The axis runs past the floor to this share of it, so a mark drawn at the floor is whole and not centred on
# the plot's edge.
FLOOR_MARGIN = 0.45
MARGIN = 16
GUTTER = 28
TITLE_HEIGHT = 26
AXIS_HEIGHT = 66
PANEL_GAP = 14
LABEL_GAP = 6
PLOT_PAD = 8
MIN_PLOT = 120
# The pixels a row takes at most and at least. Under the least pitch a dashed line would lie over the tip of the
# outline of a marker on the row's other lane: rows are never set tighter, and a figure with more rows than
# fit at the least pitch takes a taller canvas.
MAX_PITCH, MIN_PITCH = 24.0, 18.0
COUNT_PITCH = 26.0
MIN_POINTS = 7.0
# A wall-time circle's diameter and edge in points: about ten pixels across with its edge.
CIRCLE, CIRCLE_EDGE = 3.8, 0.8
# The pixels kept clear between two wall-time circles.
CIRCLE_GAP = 1.0
# The adjusted-F1 markers as (shape, size in points): a triangle for an F1 of exactly 1, a diamond below it.
MARKERS = {"exact": (">", 4.6), "adjusted": ("D", 4.2)}
# The orders the marks of an F1 row are drawn in, a marker over every bar and cross, and a marker's outline in
# points, which is drawn in the ground's colour.
BAR_ORDER, ABSENT_ORDER, MARKER_ORDER = 2, 3, 4
MARKER_EDGE = 0.5
# A marker's reach to either side of its place in units of its size: the triangle's from its centre to its tip,
# the diamond's half diagonal.
MARKER_SPANS = {">": 0.5, "D": math.sqrt(0.5)}
# The part of its outline's width that the outline covers of a marker along its lane, on the side of F1 = 0 and
# on the side of F1 = 1: half the width on a flat side, and as far as the mitre reaches inward at a corner. The
# marker's body is what its outline leaves of it.
MARKER_INSETS = {">": (0.5, math.sqrt(1.25)), "D": (math.sqrt(0.5), math.sqrt(0.5))}
# The pixels kept clear between the left edge of an F1 plot and a marker that stands at F1 = 0: the room they
# leave holds the cross of an absent plugin clear of the axis line, with ground before the tick that reads 0.
ROOM_CLEAR = 4.0
# The plain mark of a plugin a comparator has no counterpart for: its shape and its size and stroke in points,
# large enough to read as a cross at the README's width, where the image is drawn at 830 pixels of 1080. It
# stands in the room left of the tick that reads 0, where no value has a place, clear of the axis line and of
# the bar on the lane above.
ABSENT_MARK, ABSENT_SIZE, ABSENT_STROKE = "x", 3.3, 1.0
# The pixel columns of ground between a cross and the gridline of the tick that reads 0.
ABSENT_CLEAR = 1
# The dashed line from a bar's end to its tool's marker: its width in points, its dashes and gaps in widths,
# and its order, over the bars and the markers' outlines. The line starts with a dash at the bar's end, so
# the shortest stretch between a bar's end and a marker carries ink; a dash is over half again as long as the
# line is wide, so it reads as a dash and not as a dot, and a dash with its gap is short, so a line of eight
# pixels still reads as dashed. The line runs through its marker's outline to the marker's body.
JOIN_WIDTH, JOIN_DASHES, JOIN_ORDER = 0.9, (1.6, 0.9), 4.5
# The least length in widths that a dash is fitted to, where a line's dashes and gaps are fitted to its length:
# half again the line's width, under which a dash reads as a dot.
JOIN_LEAST = 1.5
# The least length in pixels at which a dash shows as one at the image's size: a join too short for a first dash
# of that length before a gap and a last dash of their least lengths is drawn solid.
JOIN_SHOWN = 1.5
# The pixels the line is drawn on into its marker's body, where it is of the body's colour, and before a marker's
# outline that lies over the bar's end, where it is of the bar's colour: a marker and a line's ends are each set
# on whole pixels, up to half a pixel from their places, and the line still meets the body and the bar.
JOIN_LAP = 1.0
# How far the edge of a marker's body, on the side of its bar, stands further from the bar per pixel off the
# lane's line: not at all on the triangle's flat side, by as much on the diamond's sides. A line is set on a
# pixel row up to half a pixel from its lane's, and its ink reaches half its width to either side.
MARKER_TAPERS = {">": 0.0, "D": 1.0}
# The pixels a marker's body stands over its bar's end where the two touch and no line is drawn. A marker is set
# on whole pixels, up to half a pixel from its place: a body that reaches less far over the bar's end may leave
# a column of ground before it, which the line covers.
JOIN_TOUCH = 0.5
# The break mark on an F1 axis before the tick that reads 1, in pixels: its centre's distance before that
# tick, the width of the cut in the axis line, a stroke's half height and its lean along the axis. Then the
# widths in points of the line that cuts the axis and of a stroke, and the order the cut is drawn in, over
# the axis line.
BREAK_PLACE, BREAK_WIDTH, BREAK_REACH, BREAK_LEAN = 6.5, 4.5, 4.0, 2.5
BREAK_CUT, BREAK_STROKE, BREAK_ORDER = 1.6, 0.7, 5
# The pixels between the legend's frame and its rows, and between two rows; then the width of an entry's mark,
# of a swatch with a circle beside it and of the space between two entries, in units of the legend's type size.
LEGEND_PAD, LEGEND_ROW_GAP = 7.0, 5.0
LEGEND_HANDLE, LEGEND_PAIR, LEGEND_SPACING = 1.2, 2.6, 1.4
# Lines of room under a panel's last row, so a text on its lower lane clears the axis.
ROW_FOOT = 0.2
# The pixels a line of a note takes, and the pixels between a note and what stands over it.
NOTE_PITCH, NOTE_GAP = 17, 8
# The lines the split of the plugin panel allows the note under the panel of counts.
NOTE_LINES = 2
STYLE = {
    "font.family": "serif",
    "font.serif": ["STIXGeneral"],
    "mathtext.fontset": "stix",
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "figure.facecolor": GROUND,
    "axes.facecolor": GROUND,
    "legend.fontsize": 7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.top": False,
    "ytick.right": False,
    "xtick.minor.visible": False,
    "ytick.minor.visible": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.linewidth": 0.6,
    "axes.grid": True,
    # A solid pale grid: an alpha'd dark one prints heavier wherever transparency is flattened.
    "grid.alpha": 1.0,
    "grid.color": GRID,
    "grid.linewidth": 0.4,
    "legend.frameon": True,
}


def decade_of(distance: Fraction) -> int:
    """The decade a positive 1 - F1 lies in: the largest k with 10^k at or under it, found exactly."""
    exponent = math.floor(math.log10(distance))
    while Fraction(10) ** exponent > distance:
        exponent -= 1
    while Fraction(10) ** (exponent + 1) <= distance:
        exponent += 1
    return exponent


def discord_floor(values: list[dict]) -> float:
    """Where an F1 of exactly 1 is drawn: the first decade under the smallest 1 - F1 among `values` that is
    not zero, the decade below for one that lies exactly on a decade, and no nearer than DISCORD_FLOOR. An F1
    below 1 so stands short of the floor, though it print as 1 at nine places."""
    below = [discordance(value) for value in values if discordance(value)]
    if not below:
        return DISCORD_FLOOR
    decade = decade_of(min(below))
    return min(DISCORD_FLOOR, float(f"1e{decade - (Fraction(10) ** decade == min(below))}"))


def discord(value: dict, floor: float) -> float:
    """An F1 value's place on the axis, from its counts and never from its printed figure: 1 - F1, and the
    floor for an F1 of exactly 1."""
    distance = discordance(value)
    return float(distance) if distance else floor


def f1_ticks(floor: float) -> list[float]:
    """The ticks of an F1 axis as 1 - F1: one per decade from F1 = 0 to the decade before the floor, then
    the floor."""
    return [float(f"1e-{k}") for k in range(round(-math.log10(floor)))] + [floor]


def f1_tick_labels(positions: list[float], floor: float = DISCORD_FLOOR) -> list[str]:
    """Positions are 1 - F1; labels are the F1 they stand for, built from the exponent: 1e-k reads "0." and k
    nines, since `1 - 1e-7` formatted naively reads "1.0". The floor, where an F1 of exactly 1 is drawn, reads 1."""
    return ["0" if pos >= 1.0 else "1" if pos <= floor else "0." + "9" * round(-math.log10(pos)) for pos in positions]


def marker_reach(shape: str) -> float:
    """The pixels a marker of `shape` reaches to either side of its place, its outline included as far as a
    mitred corner carries it."""
    marker, size = MARKERS[shape]
    return (MARKER_SPANS[marker] * size + math.sqrt(0.5) * MARKER_EDGE) * DPI / 72


def f1_room() -> float:
    """The pixels an F1 plot's left edge stands left of the tick that reads 0: the widest marker's half width
    and ROOM_CLEAR more. No value has a place in this room."""
    return max(marker_reach(shape) for shape in MARKERS) + ROOM_CLEAR


def absent_reach() -> float:
    """The pixels the cross of an absent plugin reaches to either side of its centre, its stroke included."""
    return (ABSENT_SIZE / 2 + ABSENT_STROKE * math.sqrt(0.125)) * DPI / 72


def absent_columns() -> int:
    """The pixel columns the cross of an absent plugin takes to either side of its centre: the drawing library
    sets the cross's centre between two pixels and the ends of its strokes on whole pixels, and a stroke's
    corners reach past its end."""
    ends = math.floor(ABSENT_SIZE / 2 * DPI / 72 + 0.5)
    return math.ceil(ends + ABSENT_STROKE * math.sqrt(0.125) * DPI / 72)


def absent_place() -> float:
    """The pixels the centre of an absent plugin's cross stands right of an F1 plot's left edge: a whole number
    of pixels left of the tick that reads 0. A gridline and a marker are each set on whole pixels, by one
    rounding, so wherever the plot's edge falls the cross's last column lies ABSENT_CLEAR columns of ground
    before the gridline's."""
    return f1_room() - (absent_columns() + ABSENT_CLEAR)


def f1_scale(plot_width: float, floor: float) -> float:
    """The pixels a decade takes along an F1 axis of `plot_width` pixels, which runs from the room left of
    F1 = 0 to a little past `floor`."""
    return (plot_width - f1_room()) / (math.log10(CONCORDANCE_ANCHOR) - math.log10(floor * FLOOR_MARGIN))


def marker_body(shape: str, side: int) -> float:
    """The pixels from a marker's place to the edge of its body along its lane, inside its outline: on the
    side of F1 = 0 for `side` 0 and on the side of F1 = 1 for `side` 1."""
    marker, size = MARKERS[shape]
    return (MARKER_SPANS[marker] * size - MARKER_INSETS[marker][side] * MARKER_EDGE) * DPI / 72


def join_end(tip: float, place: float, shape: str, scale: float) -> float | None:
    """Where the dashed line from a bar's end at `tip` to its marker of `shape` at `place` ends, as a value of
    1 - F1 at `scale` pixels a decade: the edge of the marker's body on the side of the bar's end, inside the
    marker's outline. None where the body stands JOIN_TOUCH pixels or more over the bar's end: the two touch,
    and no line is drawn. None too where the marker stands left of its bar's end, an adjusted F1 under the
    raw one: a line there would lie over the bar, and the marker stands alone."""
    apart = scale * math.log10(tip / place)
    reach = marker_body(shape, 0)
    if apart <= reach - JOIN_TOUCH:
        return None
    return place * 10 ** (reach / scale)


def join_drawn(tip: float, place: float, shape: str, scale: float) -> tuple[float, float] | None:
    """The values of 1 - F1 between which the dashed line from a bar's end at `tip` to its marker of `shape`
    at `place` is drawn, at `scale` pixels a decade; None where the two touch and where the marker stands
    left of the bar's end. The line is drawn on JOIN_LAP
    pixels into the marker's body, short of the marker's place. It is drawn from the bar's end; where the
    marker's outline lies over the bar's end, or within JOIN_LAP pixels of it, it is drawn from JOIN_LAP
    pixels before the outline, over the bar and in the bar's colour, so that it crosses the whole outline."""
    if join_end(tip, place, shape, scale) is None:
        return None
    apart = scale * math.log10(tip / place)
    before = marker_reach(shape) + JOIN_LAP
    stop = marker_body(shape, 0) - JOIN_LAP
    return place * 10 ** (before / scale) if apart < before else tip, place * 10 ** (stop / scale)


def join_foot(place: float, shape: str, scale: float) -> tuple[float, float] | None:
    """The values of 1 - F1 between which a join is drawn on past its end, solid, into the body of its marker
    of `shape` at `place`, at `scale` pixels a decade; None for a marker whose body's edge is as near the bar
    on every pixel row. The join ends JOIN_LAP pixels inside the body on the lane's line. A body that narrows
    towards the bar begins further on along a pixel row off that line, so the join is drawn on by as much as
    the body's edge stands further on the outermost row the line inks: it ends on the body on every such
    row, inside the body and short of the marker's place. The stretch starts half of JOIN_LAP before the
    body's edge, inside the join's last dash: a dashed join's dashes are set from its start, which is set on
    a whole pixel, so its last dash may end up to half a pixel short of the join's end."""
    marker, _ = MARKERS[shape]
    further = MARKER_TAPERS[marker] * (JOIN_WIDTH * DPI / 72 / 2 + 0.5)
    if not further:
        return None
    body = marker_body(shape, 0)
    return place * 10 ** ((body + JOIN_LAP / 2) / scale), place * 10 ** ((body - JOIN_LAP - further) / scale)


def join_dashes(length: float) -> str | tuple[float, tuple[float, float]]:
    """The line style of a join that is drawn `length` pixels on from its bar's end. A line that holds two
    dashes of JOIN_LEAST and the gap between them holds as many dashes as fit, the dashes and the gaps
    between them fitted alike from JOIN_DASHES and no dash shorter than JOIN_LEAST, so that the line starts
    with a dash at the bar's end and ends with one in the marker's body. A shorter line takes a finer
    pattern of two dashes: the last dash and the gap before it keep their least lengths, so that the gap
    shows and the last dash still reaches the body where a marker and the line's start are each set on a
    whole pixel, and the first dash is shortened. Where the first dash would be under JOIN_SHOWN pixels
    the line cannot show two dashes at the image's size: it is a solid stub. The style's lengths are in
    line widths."""
    dash, gap = JOIN_DASHES
    widths = length / (JOIN_WIDTH * DPI / 72)
    gaps = math.floor((widths * dash / JOIN_LEAST - dash) / (dash + gap))
    if gaps < 1:
        least_gap = JOIN_LEAST * gap / dash
        first = widths - JOIN_LEAST - least_gap
        if first * JOIN_WIDTH * DPI / 72 < JOIN_SHOWN:
            return "solid"
        return (JOIN_LEAST - first, (JOIN_LEAST, least_gap))
    fit = widths / ((gaps + 1) * dash + gaps * gap)
    return (0, (dash * fit, gap * fit))


def lane_offsets(lanes: list[str], pitch: float) -> tuple[float, dict[str, float]]:
    """A bar's height and each lane's offset from its row's line, in lines of `pitch` pixels. The lanes of a
    row stand a bar's height and a gap apart. Where that would bring a dashed line on one lane onto the
    vertex of a marker on the next, as it would at a pitch under the least, they stand a marker's half height
    and half the line's width apart, as far as the row holds their bars. A marker is as tall as it reaches along its
    lane: the diamond is square, and the triangle's base is as long as the triangle."""
    bar = min(0.3, 0.72 / max(len(lanes), 1))
    clear = (max(MARKER_SPANS[marker] * size for marker, size in MARKERS.values()) + JOIN_WIDTH / 2) * DPI / 72
    step = max(bar + 0.14, min(clear / pitch, (1 - bar) / max(len(lanes) - 1, 1)))
    return bar, {tool: (index - (len(lanes) - 1) / 2) * step for index, tool in enumerate(lanes)}


def row_pitch(columns: list[tuple[float, int]]) -> float:
    """The pixels a row takes in every panel of rows: the widest pitch, in half pixels and no wider than
    MAX_PITCH, at which each column's lines fit the pixels the column has for them. `columns` holds, per
    column, those pixels and the number of its lines. Rows are never set tighter than MIN_PITCH: a figure
    whose rows would need a tighter pitch is refused by name, since it takes a taller canvas."""
    pitch = math.floor(min(MAX_PITCH, *(room / lines for room, lines in columns)) * 2) / 2
    if pitch < MIN_PITCH:
        fail(
            f"the figure's rows do not fit the canvas: {max(lines for _, lines in columns)} lines in a column at {pitch:g} pixels each, "
            f"under the least pitch of {MIN_PITCH:g}; a figure with more rows takes a taller canvas, never a tighter pitch"
        )
    return pitch


def f1_edge(scale: float) -> float:
    """The left edge of an F1 plot as a value of 1 - F1, at `scale` pixels a decade: the room left of F1 = 0,
    past any value an F1 has. Every bar starts here."""
    return CONCORDANCE_ANCHOR * 10 ** (f1_room() / scale)


def fitted_labels(places: list[float], widths: list[float], gap: float) -> list[bool]:
    """Which tick labels of an F1 axis are set, given their centres and widths in pixels, the first being the
    one that reads 0 and the last the one that reads 1. The label that reads 1 is set, then the one that reads
    0, so the axis says where an F1 of 0 stands; then, taken from the end that reads 1, each label between
    them. A label is set when it stands `gap` pixels clear of every label set before it: of two labels between
    the ends that would stand closer than `gap` the one nearer 1 is kept, and a label too close to the one
    that reads 0 gives way to it."""
    kept: list[int] = []
    # The two ends first, the one that reads 1 before the one that reads 0, then the rest from the end that reads 1.
    order = dict.fromkeys([len(places) - 1, 0, *reversed(range(len(places)))]) if places else {}
    for index in order:
        low, high = places[index] - widths[index] / 2, places[index] + widths[index] / 2
        if all(high + gap <= places[k] - widths[k] / 2 or low - gap >= places[k] + widths[k] / 2 for k in kept):
            kept.append(index)
    return [index in kept for index in range(len(places))]


def f1_axis_labels(sheet, plot_widths: list[float], floor: float) -> list[str]:
    """The tick labels of the F1 axes, one set for every F1 plot, given the plots' widths in pixels: a label
    per tick, blank where it would stand closer than a label's own height to one nearer 1, or to the one that
    reads 0, on the narrowest plot. The panels of one scale so label the same decades, whatever their widths,
    and every axis labels its ticks that read 0 and 1."""
    ticks = f1_ticks(floor)
    labels = f1_tick_labels(ticks, floor)
    scale = f1_scale(min(plot_widths), floor)
    shown = fitted_labels(
        [scale * -math.log10(tick) for tick in ticks], [sheet.measure(label) for label in labels], max(sheet.rise(label) for label in labels),
    )
    return [label if keep else "" for label, keep in zip(labels, shown)]


def seconds_limits(numbers: list[float], where: str = "the wall-time panel") -> tuple[float, float]:
    """The limits of a wall-time axis: from the decade at or under the smallest value, so that no value
    lies before the axis's lowest tick, to a fifth of a decade past the largest. `where` names the panel whose
    values the numbers are; a panel that holds no value is refused by that name."""
    if not numbers:
        fail(f"{where} holds no value; a wall-time axis starts at the decade under its smallest value")
    return float(f"1e{Decimal(str(min(numbers))).adjusted()}"), max(numbers) * 10.0**0.2


def time_label(value: float) -> str:
    """A wall-time tick in short form: 1, 20, 500, 1K, 20K, 1M; a time below one as a decimal."""
    for power, suffix in ((6, "M"), (3, "K")):
        if value >= 10**power:
            return f"{value / 10**power:g}{suffix}"
    return f"{value:g}"


def seconds_ticks(low: float, high: float) -> tuple[list[float], list[str]]:
    """The ticks of a wall-time axis and their labels: the decades inside its limits; where fewer than two fall
    inside, the 1, 2 and 5 of each decade; and where fewer than two of those do, the axis's start and the
    largest value at two figures. An axis of values that lie close together has ticks to read them against."""
    decades = range(math.floor(math.log10(low)), math.ceil(math.log10(high)) + 1)
    ticks = [10.0**decade for decade in decades if low <= 10.0**decade <= high]
    if len(ticks) < 2:
        ticks = [step * 10.0**decade for decade in decades for step in (1, 2, 5) if low <= step * 10.0**decade <= high]
    if len(ticks) < 2:
        ticks = sorted({low, float(f"{high / 10.0**0.2:.2g}")})
    return ticks, [time_label(tick) for tick in ticks]


def set_time_axis(axes, numbers: list[float], where: str = "the wall-time panel") -> list[float]:
    """Give `axes` the logarithmic wall-time axis of `numbers`, the values of the panel `where` names: its
    limits, its ticks with their labels, and no minor ticks, which a logarithmic scale draws unasked. Returns
    the ticks."""
    low, high = seconds_limits(numbers, where)
    ticks, labels = seconds_ticks(low, high)
    axes.set_xscale("log")
    axes.set_xlim(low, high)
    axes.set_xticks(ticks)
    axes.set_xticklabels(labels)
    axes.minorticks_off()
    return ticks


def circle_heights(circles: list[tuple[float, float, float, str]], clear: float, reach: float) -> list[float]:
    """The height each wall-time circle is set at, in pixels down the plot, given the circles in drawing order
    as (place along the axis, lane's height, line's centre, name). A circle keeps its lane when it stands
    `clear` pixels, centre to centre, from every circle set before it. One that would meet another is set at
    the nearest height at which it stands clear, the one further from its line's centre when two are as near
    and never more than `reach` from that centre; a circle with no such height is refused by name."""
    placed: list[tuple[float, float]] = []
    for x, lane, line, name in circles:
        near = [(height, math.sqrt(clear**2 - (place - x) ** 2)) for place, height in placed if abs(place - x) < clear]
        candidates = [lane] + [height + side * rise for height, rise in near for side in (-1, 1)]
        free = [
            y for y in candidates
            if abs(y - line) <= reach and all(abs(y - height) >= rise - 1e-9 for height, rise in near)
        ]
        if not free:
            fail(f"the drawn figure has no room on its line to set the circle of {name} clear of the circles it meets")
        placed.append((x, min(free, key=lambda y: (round(abs(y - lane), 6), -abs(y - line)))))
    return [height for _, height in placed]


def split_tag(label: str) -> tuple[str, str]:
    """A row label as (name, tag): the tag is the word a synthetic-fixture row's label ends in."""
    if label.endswith(f" {SYNTHETIC_TAG}"):
        return label[: -len(SYNTHETIC_TAG) - 1], SYNTHETIC_TAG
    return label, ""


def lines_of(groups: list[dict]) -> int:
    """The lines a panel's groups take: one per row and one per group label."""
    return sum(len(group["rows"]) + bool(group["label"]) for group in groups)


def split_columns(groups: list[dict], above: int, below: int) -> int:
    """How many leading groups of the plugin panel stand in the first column: the split that keeps the taller
    column shortest, with `above` lines over both parts and `below` lines under the second."""
    best, best_height = len(groups), None
    for first in range(1, len(groups) + 1):
        height = max(above + lines_of(groups[:first]), above + lines_of(groups[first:]) + below)
        if best_height is None or height < best_height:
            best, best_height = first, height
    return best


def count_panel_height(panel: dict) -> float:
    """The pixels a panel of counts is allowed in its column: its title, its rows, its axis and the foot of a
    row."""
    return TITLE_HEIGHT + AXIS_HEIGHT + ROW_FOOT * MAX_PITCH - 17 + sum(len(group["rows"]) for group in panel["groups"]) * COUNT_PITCH


def plugin_split(panels: dict[str, dict]) -> int:
    """How many leading groups of the plugin panel stand in the first column, under the dataset panel; the
    others continue the panel in the second column, over the panels of counts: (d) with its note and, under
    it, (e). A function of the panels by id alone, so the texts can say where the panel continues: the note
    under panel (d) is allowed NOTE_LINES lines. 0 without a plugin panel."""
    groups = panels["c"]["groups"] if "c" in panels else []
    if not groups:
        return 0
    under = 0.0
    if "d" in panels:
        under = count_panel_height(panels["d"])
        under += NOTE_LINES * NOTE_PITCH + NOTE_GAP if panels["d"].get("note") else 0
    if "e" in panels:
        under += (PANEL_GAP if under else 0) + count_panel_height(panels["e"])
    return split_columns(groups, lines_of(panels["a"]["groups"]), round(under / MIN_PITCH))


def column_loads(panels: dict[str, dict], first: int, note_lines: int) -> list[tuple[float, int]]:
    """Per column, the pixels its titles, its axes and its panels of counts take and the lines of its rows,
    with the plugin panel's `first` leading groups in the first column and `note_lines` lines of note under
    panel (d): the first column holds (a) and the plugin panel's first part, the second (b), the plugin
    panel's other part, (d) with its note and (e)."""
    groups = panels["c"]["groups"] if "c" in panels else []
    dataset_lines = lines_of(panels["a"]["groups"])
    fixed = TITLE_HEIGHT + AXIS_HEIGHT + ROW_FOOT * MAX_PITCH
    under = 0.0
    if "d" in panels:
        under += PANEL_GAP + count_panel_height(panels["d"]) + (note_lines * NOTE_PITCH + NOTE_GAP if note_lines else 0)
    if "e" in panels:
        under += PANEL_GAP + count_panel_height(panels["e"])
    return [
        (fixed + (PANEL_GAP + fixed if part else 0) + (under if index == 1 else 0), dataset_lines + lines_of(part))
        for index, part in enumerate((groups[:first], groups[first:]))
    ]


class Sheet:
    """The canvas in pixels from its top-left corner."""

    def __init__(self, figure):
        self.figure = figure
        self.renderer = figure.canvas.get_renderer()
        self.width, self.height = CANVAS

    def axes(self, left: float, top: float, width: float, height: float):
        return self.figure.add_axes([left / self.width, 1 - (top + height) / self.height, width / self.width, height / self.height])

    def text(self, x: float, y: float, content: str, **style):
        style.setdefault("color", INK)
        style.setdefault("fontsize", MIN_POINTS)
        style.setdefault("va", "center")
        return self.figure.text(x / self.width, 1 - y / self.height, content, **style)

    def measure(self, content: str, **style) -> float:
        probe = self.text(0, 0, content, **style)
        width = probe.get_window_extent(self.renderer).width
        probe.remove()
        return width

    def rise(self, content: str, **style) -> float:
        """The height of a text's glyphs in pixels; a text's extent is a whole line high whatever it holds."""
        probe = self.text(0, 0, content, **style)
        _, height, _ = self.renderer.get_text_width_height_descent(content, probe.get_fontproperties(), ismath=False)
        probe.remove()
        return height

    def wrap(self, content: str, limit: float, **style) -> list[str]:
        """Greedy word wrap by measured width."""
        lines: list[str] = []
        for word in content.split():
            if lines and self.measure(f"{lines[-1]} {word}", **style) <= limit:
                lines[-1] += f" {word}"
            else:
                lines.append(word)
        return lines


def part_heading(panel: dict, continued: bool) -> tuple[str | None, str]:
    """The letter and the title over a part of a panel: the panel's own over its first part, or its only one;
    over a part that continues it, no letter and the title as the first part's continued."""
    return (None, f"{panel['title']}{CONTINUED}") if continued else (panel["id"], panel["title"])


def draw_title(sheet, left: float, top: float, panel_id: str | None, title: str) -> None:
    """A panel's letter and title over its plot; a part that continues a panel has no letter, and its title
    starts where the letter would."""
    if panel_id is None:
        sheet.text(left, top + TITLE_HEIGHT / 2 - 2, title, fontsize=8, ha="left")
        return
    letter = f"({panel_id})"
    sheet.text(left, top + TITLE_HEIGHT / 2 - 2, letter, fontsize=9, fontweight="bold", ha="left")
    sheet.text(left + sheet.measure(letter, fontsize=9, fontweight="bold") + 6, top + TITLE_HEIGHT / 2 - 2, title, fontsize=8, ha="left")


def draw_axis_label(sheet, centre: float, top: float, lines: list[str]) -> None:
    """The value-axis label under a plot, centred on the plot and kept on the canvas."""
    for index, line in enumerate(lines):
        half = sheet.measure(line, fontsize=8) / 2
        x = min(max(centre, MARGIN / 2 + half), sheet.width - MARGIN / 2 - half)
        sheet.text(x, top + 36 + index * 17, line, fontsize=8, ha="center")


def label_column(sheet, groups: list[dict]) -> float:
    """The width a panel's labels take left of its plot: the widest row label or group label of `groups`, with
    its gap to the plot. A group label stands in this column too, so it ends before the plot begins."""
    widths = [sheet.measure(row["label"]) for group in groups for row in group["rows"]]
    widths += [sheet.measure(group["label"], fontstyle="italic") for group in groups if group["label"]]
    return (max(widths) if widths else 0) + LABEL_GAP + 4


def column_widths(sheet, content: float, parts: tuple[list[dict], list[dict]]) -> tuple[float, float]:
    """The widths of the two columns across `content` pixels. With the plugin panel in both, the columns are
    sized so that its two plots have one width, and so one scale; otherwise the columns are equal."""
    room = content - GUTTER
    if not (parts[0] and parts[1]):
        return room / 2, room / 2
    labels = [label_column(sheet, part) for part in parts]
    plot = (room - sum(labels) - 2 * PLOT_PAD) / 2
    return labels[0] + plot + PLOT_PAD, labels[1] + plot + PLOT_PAD


def plot_box(left: float, width: float, label_width: float, title: str) -> tuple[float, float]:
    """A plot's left edge and width inside its column; a plot narrower than MIN_PLOT is refused by name, one
    whose labels take the whole column as left no room."""
    plot_width = width - label_width - PLOT_PAD
    if plot_width < MIN_PLOT:
        room = f"{plot_width:.0f} pixels" if plot_width >= 0.5 else "no room"
        fail(f"the drawn figure leaves the plot of {title!r} {room} beside its labels, under the minimum of {MIN_PLOT}")
    return left + label_width, plot_width


def draw_rows(sheet: Sheet, axes, left: float, plot_left: float, plot_top: float, pitch: float, groups: list[dict]) -> dict:
    """The group labels and row labels of a panel whose rows run down the page; returns each row's line index
    by (group index, row index)."""
    lines, line = {}, 0
    for g, group in enumerate(groups):
        if group["label"]:
            y = plot_top + (line + 0.5) * pitch
            sheet.text(left, y, group["label"], color=INK_SECONDARY, fontstyle="italic", ha="left")
            line += 1
        for r, row in enumerate(group["rows"]):
            y = plot_top + (line + 0.5) * pitch
            name, tag = split_tag(row["label"])
            x = plot_left - LABEL_GAP
            if tag:
                sheet.text(x, y, tag, color=INK_SECONDARY, ha="right")
                x -= sheet.measure(f" {tag}")
            sheet.text(x, y, name, ha="right")
            lines[(g, r)] = line
            line += 1
    axes.set_ylim(line - 0.5 + ROW_FOOT, -0.5)
    axes.set_yticks([])
    axes.grid(axis="y", visible=False)
    axes.tick_params(axis="x", length=2.5, width=0.5, pad=2)
    return lines


def draw_f1_marks(axes, rows: list[tuple[float, dict]], offsets: dict[str, float], fills: dict[str, str], bar: float,
                  floor: float, scale: float) -> None:
    """The marks of an F1 panel's rows, each row given with its line: a bar per tool from the plot's left
    edge and the adjusted F1 as a marker, each on its tool's lane; a marker of exactly 1 stands at the floor,
    at the end of a full bar too. A marker is drawn whole and over every bar, with a
    thin outline in the ground's colour: no gap is cut into it where it stands over a bar's tip. A thin dashed
    line in the tool's colour joins the end of its bar to its marker wherever the marker stands right of the
    bar's end and the two do not touch, and a marker left of its bar's end has none, by one
    rule for every tool: a marker touches its bar where its body, inside its outline, stands half a pixel or
    more over the bar's end, and there no line is drawn; the line is drawn at any stretch between the two, a solid stub where the
    stretch is short, and it starts with a dash at the bar's end and ends with one at the marker's body,
    through the marker's outline, so that no ground stands before the body; into a diamond, whose body
    narrows towards the bar, it is drawn on solid so that it ends on the body on every pixel row it inks. A bar ends at its counts' place, however
    near 0 that is, and is seen left of its own marker there: the edge it starts at stands a marker's half
    width and more left of F1 = 0. `scale` is the pixels a decade takes along the axis."""
    edge = f1_edge(scale)
    for line, row in rows:
        tips = {value["tool"]: discord(value, floor) for value in row["values"] if value["kind"] == "bar"}
        for value in row["values"]:
            if value["kind"] == "text":
                continue
            x, y = discord(value, floor), line + offsets[value["tool"]]
            if value["kind"] == "bar":
                # A negative width with `left` at the edge keeps both endpoints explicit, which a
                # logarithmic scale requires.
                axes.barh(y, tips[value["tool"]] - edge, left=edge, height=bar, color=fills[value["tool"]], linewidth=0,
                          zorder=BAR_ORDER)
            else:
                shape = marker_shape(value)
                tip = tips.get(value["tool"])
                drawn = join_drawn(tip, x, shape, scale) if tip is not None else None
                if drawn is not None:
                    axes.plot(list(drawn), [y, y], color=fills[value["tool"]], linewidth=JOIN_WIDTH,
                              linestyle=join_dashes(scale * abs(math.log10(tip / drawn[1]))), dash_capstyle="butt", solid_capstyle="butt",
                              zorder=JOIN_ORDER)
                    foot = join_foot(x, shape, scale)
                    if foot is not None:
                        axes.plot(list(foot), [y, y], color=fills[value["tool"]], linewidth=JOIN_WIDTH, linestyle="solid", solid_capstyle="butt",
                                  zorder=JOIN_ORDER)
                axes.plot([x], [y], marker=MARKERS[shape][0], markersize=MARKERS[shape][1], linestyle="none",
                          markerfacecolor=fills[value["tool"]], markeredgecolor=GROUND, markeredgewidth=MARKER_EDGE, zorder=MARKER_ORDER,
                          clip_on=False)


def draw_absent(axes, across, rows: list[tuple[float, dict]], offsets: dict[str, float], fills: dict[str, str], width: float) -> None:
    """The plain mark of each row a record lists as absent for a tool: a cross in the tool's colour on the
    tool's lane, in the room left of the tick that reads 0 with ABSENT_CLEAR pixel columns of ground before that
    tick's gridline, so the lane is seen to be the tool's and the cross stands where no figure can. `across`
    places a point by its share of the plot's width, which is `width` pixels, and its line."""
    for line, row in rows:
        for value in row["values"]:
            if value["kind"] == "text":
                axes.plot([absent_place() / width], [line + offsets[value["tool"]]], transform=across, marker=ABSENT_MARK, markersize=ABSENT_SIZE,
                          markeredgewidth=ABSENT_STROKE, linestyle="none", color=fills[value["tool"]], zorder=ABSENT_ORDER, clip_on=False)


def draw_axis_break(axes, along, floor: float, scale: float, height: float) -> None:
    """The break mark on an F1 axis before the tick that reads 1, which has no place on an axis of 1 - F1 in
    decades: the axis line is cut and a leaning stroke marks each end of the cut. `along` places a point by
    its value along the axis and its share of the plot's height, `scale` is the pixels a decade takes and
    `height` the plot's height in pixels."""
    def before(pixels: float) -> float:
        return floor * 10 ** (pixels / scale)

    ends = (BREAK_PLACE + BREAK_WIDTH / 2, BREAK_PLACE - BREAK_WIDTH / 2)
    axes.plot([before(end) for end in ends], [0, 0], transform=along, color=GROUND, linewidth=BREAK_CUT, solid_capstyle="butt",
              zorder=BREAK_ORDER, clip_on=False)
    for end in ends:
        axes.plot([before(end + BREAK_LEAN), before(end - BREAK_LEAN)], [-BREAK_REACH / height, BREAK_REACH / height], transform=along,
                  color=INK, linewidth=BREAK_STROKE, solid_capstyle="butt", zorder=BREAK_ORDER + 1, clip_on=False)


def draw_f1(sheet: Sheet, document: dict, tools: list[dict], panel: dict, groups: list[dict], heading: tuple[str | None, str],
            left: float, top: float, width: float, pitch: float, floor: float, tick_labels: list[str]) -> float:
    """A panel of raw F1 bars per row, or one part of it under `heading`, its letter and title: one lane per
    tool, the adjusted F1 as a marker on its tool's lane and a plugin a tool has no counterpart for as a
    plain cross. An F1 of exactly
    1 is drawn at `floor`, whose tick stands behind a break mark; `tick_labels` are the labels of the axis's
    ticks, the same in every F1 panel. Returns the height taken."""
    from matplotlib.transforms import blended_transform_factory

    fills = {tool["id"]: tool["colour"] for tool in document["tools"]}
    lanes = [t for t in lane_tools(document, panel) if any(r["id"] == t and r["role"] != "reference" for r in tools)]
    letter, title = heading
    draw_title(sheet, left, top, letter, title)
    (plot_left, plot_width), plot_top = plot_box(left, width, label_column(sheet, groups), title), top + TITLE_HEIGHT
    plot_height = (lines_of(groups) + ROW_FOOT) * pitch
    axes = sheet.axes(plot_left, plot_top, plot_width, plot_height)
    lines = draw_rows(sheet, axes, left, plot_left, plot_top, pitch, groups)
    axes.set_xscale("log")
    scale = f1_scale(plot_width, floor)
    # Descending limits reverse the axis, which starts the room left of F1 = 0.
    axes.set_xlim(f1_edge(scale), floor * FLOOR_MARGIN)
    bar, offsets = lane_offsets(lanes, pitch)
    rows = [(lines[(g, r)], row) for g, group in enumerate(groups) for r, row in enumerate(group["rows"])]
    draw_absent(axes, blended_transform_factory(axes.transAxes, axes.transData), rows, offsets, fills, plot_width)
    draw_f1_marks(axes, rows, offsets, fills, bar, floor, scale)
    # A tick per decade and one at the floor.
    axes.set_xticks(f1_ticks(floor))
    axes.set_xticklabels(tick_labels)
    draw_axis_break(axes, blended_transform_factory(axes.transData, axes.transAxes), floor, scale, plot_height)
    draw_axis_label(sheet, plot_left + plot_width / 2, plot_top + plot_height, axis_lines(document, tools, panel))
    return TITLE_HEIGHT + plot_height + AXIS_HEIGHT


def draw_seconds(sheet: Sheet, document: dict, tools: list[dict], panel: dict, left: float, top: float, width: float,
                 plot_height: float) -> float:
    """A panel of wall times: one line per dataset, its rows as lanes (the first filled, the second open), one
    marker per tool in the tool's colour, which the legend names, and circles that would meet set apart along
    the line's height. No tool is named in the panel. A panel that holds no group, or no value, is refused by
    name. Returns the height taken."""
    fills = {tool["id"]: tool["colour"] for tool in document["tools"]}
    names = {tool["id"]: tool["name"] for tool in document["tools"]}
    groups = panel["groups"]
    where = f"the values file's panels[id={panel['id']}]"
    if not groups:
        fail(f"{where} holds no group; the wall-time panel draws one line per dataset")
    draw_title(sheet, left, top, panel["id"], panel["title"])
    label_width = max(sheet.measure(group["label"]) for group in groups) + LABEL_GAP + 4
    (plot_left, plot_width), plot_top = plot_box(left, width, label_width, panel["title"]), top + TITLE_HEIGHT
    pitch = plot_height / len(groups)
    axes = sheet.axes(plot_left, plot_top, plot_width, plot_height)
    axes.set_ylim(len(groups) - 0.5, -0.5)
    axes.set_yticks([])
    axes.grid(axis="y", visible=False)
    axes.tick_params(axis="x", length=2.5, width=0.5, pad=2)
    set_time_axis(axes, [float(v["value"]) for group in groups for row in group["rows"] for v in row["values"]], where)
    for g, group in enumerate(groups):
        sheet.text(plot_left - LABEL_GAP, plot_top + (g + 0.5) * pitch, group["label"], ha="right")
    marks = [(g, r, value) for g, group in enumerate(groups) for r, row in enumerate(group["rows"]) for value in row["values"]]
    # Two tools whose times lie close on one lane would stand one over the other: such circles are set apart
    # along the line's height, each still at its time.
    heights = circle_heights([
        (axes.transData.transform((float(value["value"]), 0))[0], (g + 0.5 + (r - (len(groups[g]["rows"]) - 1) / 2) * 0.36) * pitch,
         (g + 0.5) * pitch, f"{names[value['tool']]} on {groups[g]['label']} ({groups[g]['rows'][r]['label']})")
        for g, r, value in marks
    ], (CIRCLE + CIRCLE_EDGE) * DPI / 72 + CIRCLE_GAP, pitch / 2)
    for (_, r, value), height in zip(marks, heights):
        colour, y = fills[value["tool"]], height / pitch - 0.5
        # A ring in the ground's colour under the circles keeps each off the gridlines.
        axes.plot([float(value["value"])], [y], marker="o", markersize=CIRCLE + CIRCLE_EDGE + 1.4, linestyle="none",
                  markerfacecolor=GROUND, markeredgewidth=0, zorder=2.5, clip_on=False)
        axes.plot([float(value["value"])], [y], marker="o", markersize=CIRCLE, linestyle="none", markeredgewidth=CIRCLE_EDGE,
                  markerfacecolor=colour if r == 0 else GROUND, markeredgecolor=colour, zorder=3, clip_on=False)
    draw_axis_label(sheet, plot_left + plot_width / 2, plot_top + plot_height, axis_lines(document, tools, panel))
    return TITLE_HEIGHT + plot_height + AXIS_HEIGHT


def count_ticks(limit: int) -> list[int]:
    """The ticks of a count axis that runs from 0 to `limit`: the multiples, from 0, of the least step of 1, 2
    or 5 times a power of ten that gives seven ticks or fewer."""
    steps = sorted(factor * 10**power for power in range(len(str(limit)) + 1) for factor in (1, 2, 5))
    step = next(step for step in steps if limit // step < 7)
    return list(range(0, limit + 1, step))


def draw_count_cross(axes, across, line: float, row: dict, fills: dict[str, str], width: float) -> None:
    """The cross of a row of a panel of counts that holds a text and no bar, and that text beside it. The
    cross is the plugin panel's, as many pixels right of the plot's left edge, where a bar would start: it
    stands clear of the axis line, in its tool's colour, and the text starts a cross's reach and LABEL_GAP
    pixels right of its centre. `across` places a point by its share of the plot's width, which is `width`
    pixels, and its line."""
    draw_absent(axes, across, [(line, row)], {value["tool"]: 0.0 for value in row["values"]}, fills, width)
    for value in row["values"]:
        if value["kind"] == "text":
            axes.annotate(value["value"], xy=(absent_place() / width, line), xycoords=across, xytext=((absent_reach() + LABEL_GAP) * 72 / DPI, 0),
                          textcoords="offset points", fontsize=MIN_POINTS, color=INK, ha="left", va="center", annotation_clip=False)


def draw_count(sheet: Sheet, document: dict, tools: list[dict], panel: dict, left: float, top: float, width: float) -> float:
    """A panel of counts: one solid bar per row, labelled at its tip, and a cross with its text on a row
    that holds a text and no bar. A count of a thousand or more is written with thousands separators, on the
    axis as in a label. Returns the height taken."""
    fills = {tool["id"]: tool["colour"] for tool in document["tools"]}
    rows = [row for group in panel["groups"] for row in group["rows"]]
    crossed = {id(row): value["value"] for row, value in crossed_rows(panel)}
    draw_title(sheet, left, top, panel["id"], panel["title"])
    label_width = max(sheet.measure(row["label"]) for row in rows) + LABEL_GAP + 4
    tip_labels = [value["value"] for row in rows if id(row) not in crossed for value in row["values"] if value["kind"] == "text"]
    tip_width = max((sheet.measure(label) for label in tip_labels), default=0) + 10
    (plot_left, plot_width), plot_top = plot_box(left, width - tip_width, label_width, panel["title"]), top + TITLE_HEIGHT
    plot_height = len(rows) * COUNT_PITCH
    axes = sheet.axes(plot_left, plot_top, plot_width, plot_height)
    axes.set_ylim(len(rows) - 0.5, -0.5)
    axes.set_yticks([])
    axes.grid(axis="y", visible=False)
    axes.tick_params(axis="x", length=2.5, width=0.5, pad=2)
    counts = [int(value["value"]) for row in rows for value in row["values"] if value["kind"] == "bar"]
    totals = [int(found.group(1).replace(",", "")) for label in tip_labels if (found := re.fullmatch(r"[\d,]+ of ([\d,]+)", label))]
    axes.set_xlim(0, max(counts + totals))
    if max(counts + totals) >= 1000:
        ticks = count_ticks(max(counts + totals))
        axes.set_xticks(ticks, [f"{tick:,}" for tick in ticks])
    for text in crossed.values():
        if sheet.measure(text) > plot_width + tip_width - (absent_place() + absent_reach() + LABEL_GAP):
            fail(f"the drawn figure leaves {text!r} no room beside its cross in the plot of {panel['title']!r}")
    for index, row in enumerate(rows):
        sheet.text(plot_left - LABEL_GAP, plot_top + (index + 0.5) * COUNT_PITCH, row["label"], ha="right")
        if id(row) in crossed:
            draw_count_cross(axes, axes.get_yaxis_transform(), index, row, fills, plot_width)
            continue
        for value in row["values"]:
            if value["kind"] == "bar":
                axes.barh(index, int(value["value"]), height=0.56, color=fills[value["tool"]], linewidth=0, zorder=2)
                tip = int(value["value"])
        for value in row["values"]:
            if value["kind"] == "text":
                axes.annotate(value["value"], xy=(tip, index), xytext=(4, 0), textcoords="offset points", fontsize=MIN_POINTS,
                              color=INK, ha="left", va="center", annotation_clip=False)
    draw_axis_label(sheet, plot_left + plot_width / 2, plot_top + plot_height, axis_lines(document, tools, panel))
    return TITLE_HEIGHT + plot_height + AXIS_HEIGHT - 17


def legend_handle(entry: dict, colour: str):
    """The mark a legend entry shows: a swatch for a tool with bars and for what a bar states, a circle for a
    tool with circles alone or for an architecture, a swatch with a circle in the same colour beside it for a
    tool with both, and the marker or the cross of a class of mark."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    circle = Line2D([], [], marker="o", markersize=CIRCLE, linestyle="none", markeredgewidth=CIRCLE_EDGE, markeredgecolor=colour,
                    markerfacecolor=GROUND if entry["mark"] == "open" else colour, label=entry["label"])
    if entry["mark"] == "fill_dot":
        return (Patch(facecolor=colour, edgecolor="none"), circle)
    if entry["mark"] in ("fill", "raw"):
        return Patch(facecolor=colour, edgecolor="none", label=entry["label"])
    if entry["mark"] in MARKERS:
        marker, size = MARKERS[entry["mark"]]
        return Line2D([], [], marker=marker, markersize=size, linestyle="none", markerfacecolor=colour, markeredgecolor=GROUND,
                      markeredgewidth=MARKER_EDGE, label=entry["label"])
    if entry["mark"] == "absent":
        return Line2D([], [], marker=ABSENT_MARK, markersize=ABSENT_SIZE, markeredgewidth=ABSENT_STROKE, linestyle="none", color=colour,
                      label=entry["label"])
    return circle


def legend_runs(row: list[dict]) -> list[list[dict]]:
    """A legend row's entries as runs of neighbours whose marks take one width: the entries with a swatch and a
    circle side by side, which take two marks' width, apart from those with one mark. Each run is set with its
    own mark width, so the gap from an entry's mark to the entry before it is the same along the row."""
    runs: list[list[dict]] = []
    for entry in row:
        if runs and (runs[-1][-1]["mark"] == "fill_dot") == (entry["mark"] == "fill_dot"):
            runs[-1].append(entry)
        else:
            runs.append([entry])
    return runs


def draw_legend(sheet: Sheet, document: dict, tools: list[dict], top: float) -> float:
    """The legend across the top of the canvas, in one frame: a row per kind of entry, each centred, so no
    entry of one kind stands in a column under an entry of another. A row is set as its runs side by side, a
    column's spacing apart. Returns the height taken."""
    from matplotlib.legend_handler import HandlerTuple
    from matplotlib.patches import Rectangle

    fills = {tool["id"]: tool["colour"] for tool in document["tools"]}
    below, widest = top + LEGEND_PAD, 0.0
    spacing = LEGEND_SPACING * MIN_POINTS * DPI / 72
    for row in legend_rows(legend_entries(document)):
        legends = [
            sheet.figure.legend(
                handles=[legend_handle(entry, fills.get(entry["tool"], "0.35")) for entry in run], labels=[entry["label"] for entry in run],
                handler_map={tuple: HandlerTuple(ndivide=None, pad=0)}, loc="upper left", bbox_to_anchor=(0, 1 - below / sheet.height),
                ncol=len(run), fontsize=MIN_POINTS, labelcolor=INK, frameon=False, borderpad=0, borderaxespad=0,
                handlelength=LEGEND_PAIR if run[0]["mark"] == "fill_dot" else LEGEND_HANDLE, handletextpad=0.5, columnspacing=LEGEND_SPACING,
            ) for run in legend_runs(row)
        ]
        sheet.figure.canvas.draw()
        boxes = [legend.get_window_extent(sheet.renderer) for legend in legends]
        width = sum(box.width for box in boxes) + spacing * (len(boxes) - 1)
        left = (sheet.width - width) / 2
        for legend, box in zip(legends, boxes):
            # In the figure's own coordinates, as the frame is set: an anchor given without a transform keeps to
            # the uncropped canvas, and the run stands apart from its frame once the save centres its window.
            legend.set_bbox_to_anchor((left / sheet.width, 1 - below / sheet.height), transform=sheet.figure.transFigure)
            left += box.width + spacing
        below, widest = below + max(box.height for box in boxes) + LEGEND_ROW_GAP, max(widest, width)
    height = below - LEGEND_ROW_GAP + LEGEND_PAD - top
    sheet.figure.add_artist(Rectangle(
        (((sheet.width - widest) / 2 - LEGEND_PAD) / sheet.width, 1 - (top + height) / sheet.height),
        (widest + 2 * LEGEND_PAD) / sheet.width, height / sheet.height, transform=sheet.figure.transFigure, facecolor="none",
        edgecolor=GRID, linewidth=0.5,
    ))
    return height


def hidden_tick_labels(axes_list: list) -> set[int]:
    """The ids of the tick labels that lie outside their axis's view: an axis keeps a label for every tick its
    locator makes and draws those inside the view alone."""
    hidden = set()
    for axes in axes_list:
        for axis in (axes.xaxis, axes.yaxis):
            low, high = sorted(axis.get_view_interval())
            for tick in (*axis.get_major_ticks(), *axis.get_minor_ticks()):
                location = tick.get_loc()
                if not (low <= location <= high or math.isclose(location, low, rel_tol=1e-9) or math.isclose(location, high, rel_tol=1e-9)):
                    hidden.update((id(tick.label1), id(tick.label2)))
    return hidden


def placed_texts(texts: list, axes_list: list, renderer) -> list[tuple]:
    """The texts a reader sees, each with its extent: visible, not blank, and no label of a tick outside its
    axis's view."""
    hidden = hidden_tick_labels(axes_list)
    return [
        (artist, artist.get_window_extent(renderer))
        for artist in texts if artist.get_visible() and artist.get_text().strip() and id(artist) not in hidden
    ]


def ink_contrast(colour: str, ground: str = GROUND) -> float:
    """The contrast of a colour against its ground, each as #rrggbb: the lighter relative luminance plus 0.05
    over the darker plus 0.05, from 1 for one colour to 21 for black on white."""
    def luminance(shade: str) -> float:
        channels = [int(shade[index:index + 2], 16) / 255 for index in (1, 3, 5)]
        red, green, blue = (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in channels)
        return 0.2126 * red + 0.7152 * green + 0.0722 * blue

    darker, lighter = sorted((luminance(colour), luminance(ground)))
    return (lighter + 0.05) / (darker + 0.05)


def check_texts(texts: list[tuple]) -> None:
    """Refuse type under the minimum size, a text under the contrast of text against the ground, and two texts
    that overlap, given each text a reader sees as (content, size in points, colour as #rrggbb, extent). The
    contrast is that of the colours as set, the text's and the ground's."""
    for content, size, colour, _ in texts:
        if size < MIN_POINTS:
            fail(f"the drawn figure sets {content!r} at {size:g} points, under the minimum of {MIN_POINTS:g}")
        if ink_contrast(colour) < MIN_CONTRAST:
            fail(
                f"the drawn figure sets {content!r} at a contrast of {ink_contrast(colour):.1f} to 1 against its ground, under the "
                f"minimum of {MIN_CONTRAST:g}"
            )
    for index, (content, _, _, box) in enumerate(texts):
        for other, _, _, other_box in texts[index + 1:]:
            if min(box.x1, other_box.x1) - max(box.x0, other_box.x0) > 1 and min(box.y1, other_box.y1) - max(box.y0, other_box.y0) > 1:
                fail(f"the drawn figure sets {content!r} over {other!r}")


def check_letters(strings: list[str], letters: list[str]) -> None:
    """Refuse a drawing that sets one of the panels' `letters` more than once or not at all, given every string
    it sets: a panel carries its letter once, over its first part."""
    for letter in letters:
        if strings.count(letter) != 1:
            fail(f"the drawn figure sets the letter {letter} {strings.count(letter)} times; a panel carries its letter once")


def check_drawing(sheet: Sheet, right_spine: float, note_artists: list, expected: list[str], letters: list[str]) -> None:
    """Refuse a drawing with type under the minimum size or the contrast of text, two texts that overlap, a
    note that ends past the plots' right edge, one without one of the `expected` strings, or one that does not
    set each of the panels' `letters` once."""
    from matplotlib.colors import to_hex
    from matplotlib.text import Text

    placed = placed_texts(sheet.figure.findobj(Text), sheet.figure.axes, sheet.renderer)
    check_texts([(artist.get_text(), artist.get_fontsize(), to_hex(artist.get_color()), box) for artist, box in placed])
    for artist in note_artists:
        if artist.get_window_extent(sheet.renderer).x1 > right_spine + 0.5:
            fail(f"the drawn figure's note {artist.get_text()!r} ends past the plots' right edge")
    check_letters([artist.get_text() for artist, _ in placed], letters)
    # A row label ending in a tag is set as two texts, and a note as its wrapped lines.
    strings = {artist.get_text() for artist, _ in placed}
    wrapped = " ".join(artist.get_text() for artist in note_artists)
    for content in expected:
        name, tag = split_tag(content)
        if content not in strings and not (tag and {name, tag} <= strings) and content not in wrapped:
            fail(f"the drawn figure lacks {content!r}")


def draw(document: dict, tools: list[dict], png: Path, digest: str) -> None:
    """Draw the values document on the canvas and write the PNG with `digest` in a text chunk. The layout is
    two columns: (a) over the plugin panel's first groups, (b) over its other groups, (d) with its note
    under it and (e) under that, any other note under both columns and the legend above. The columns are sized so that the plugin panel's two plots share one scale; the
    tight crop is padded to the canvas."""
    try:
        import matplotlib
    except ImportError:
        fail("the drawing layer needs the packages of scripts/docs/requirements-figure.txt; --values-only writes the values file without them")
    matplotlib.use("Agg")
    import warnings

    from matplotlib import rc_context
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.transforms import Bbox

    width, height = CANVAS
    panels = {panel["id"]: panel for panel in document["panels"]}
    with rc_context(STYLE), warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*[Gg]lyph.*")
        figure = Figure(figsize=(width / DPI, height / DPI), dpi=DPI, facecolor=GROUND)
        FigureCanvasAgg(figure)
        sheet = Sheet(figure)
        plugin_groups = panels["c"]["groups"] if "c" in panels else []
        dataset_lines = lines_of(panels["a"]["groups"])
        panel_note = panels["d"].get("note") if "d" in panels else None

        # The plugin panel is split by the panels alone, as the alternative text says where it continues.
        first = plugin_split(panels)
        parts = (plugin_groups[:first], plugin_groups[first:])
        columns = column_widths(sheet, width - 2 * MARGIN, parts)
        panel_note_lines = sheet.wrap(panel_note, columns[1] - PLOT_PAD, color=INK_SECONDARY) if panel_note else []
        lefts = (MARGIN, MARGIN + columns[0] + GUTTER)
        right_spine = lefts[1] + columns[1] - PLOT_PAD

        top = float(MARGIN)
        notes = list(document["notes"])
        if notes and notes[0] == DRAFT_LINE:
            sheet.text(MARGIN, top + 9, notes.pop(0), fontsize=9, fontweight="bold", ha="left")
            top += 26
        top += draw_legend(sheet, document, tools, top) + 12
        note_lines = [line for note in notes for line in sheet.wrap(note, right_spine - MARGIN, color=INK_SECONDARY)]
        bottom = height - MARGIN - (len(note_lines) * NOTE_PITCH + NOTE_GAP if note_lines else 0)

        # Per column, the pixels left for its rows under the titles, the axes and the panels of counts.
        pitch = row_pitch([(bottom - top - used, lines) for used, lines in column_loads(panels, first, len(panel_note_lines))])

        floor = discord_floor([
            value for panel in document["panels"] if PANEL_AXES[panel["id"]] == "f1"
            for group in panel["groups"] for row in group["rows"] for value in row["values"] if value["kind"] != "text"
        ])
        # One set of tick labels serves every F1 plot: the set that fits the narrowest.
        f1_plots = [(panels["a"]["groups"], panels["a"]["title"], 0)] + [(part, panels["c"]["title"], index) for index, part in enumerate(parts) if part]
        tick_labels = f1_axis_labels(
            sheet, [plot_box(lefts[index], columns[index], label_column(sheet, groups), title)[1] for groups, title, index in f1_plots], floor,
        )
        tops = [top, top]
        tops[0] += draw_f1(
            sheet, document, tools, panels["a"], panels["a"]["groups"], part_heading(panels["a"], False), lefts[0], top, columns[0], pitch, floor,
            tick_labels,
        )
        tops[1] += draw_seconds(sheet, document, tools, panels["b"], lefts[1], top, columns[1], (dataset_lines + ROW_FOOT) * pitch)
        for index, part in enumerate(parts):
            if part:
                # The second part continues the first: it carries the title as continued, and no second letter.
                heading = part_heading(panels["c"], bool(index and parts[0]))
                tops[index] += PANEL_GAP
                tops[index] += draw_f1(
                    sheet, document, tools, panels["c"], part, heading, lefts[index], tops[index], columns[index], pitch, floor, tick_labels,
                )
        note_artists = []
        if "d" in panels:
            tops[1] += PANEL_GAP
            tops[1] += draw_count(sheet, document, tools, panels["d"], lefts[1], tops[1], columns[1])
            note_artists += [
                sheet.text(lefts[1], tops[1] + NOTE_GAP + index * NOTE_PITCH, line, color=INK_SECONDARY, ha="left")
                for index, line in enumerate(panel_note_lines)
            ]
            tops[1] += len(panel_note_lines) * NOTE_PITCH + NOTE_GAP if panel_note_lines else 0
        if "e" in panels:
            tops[1] += PANEL_GAP
            tops[1] += draw_count(sheet, document, tools, panels["e"], lefts[1], tops[1], columns[1])
        note_top = max(tops) + NOTE_GAP
        note_artists += [
            sheet.text(MARGIN, note_top + index * NOTE_PITCH + NOTE_GAP, line, color=INK_SECONDARY, ha="left") for index, line in enumerate(note_lines)
        ]

        figure.canvas.draw()
        check_drawing(sheet, right_spine, note_artists, drawn_texts(document, tools), [f"({panel['id']})" for panel in document["panels"]])
        tight = figure.get_tightbbox(sheet.renderer)
        if tight.width > width / DPI or tight.height > height / DPI:
            fail(f"the drawn figure is {tight.width * DPI:.0f} by {tight.height * DPI:.0f} pixels, past the canvas of {width} by {height}")
        # The tight crop is padded to the canvas: the window is the canvas's size, centred on the crop. The
        # quarter pixel keeps the pixel count from rounding down.
        window = Bbox.from_bounds(
            (tight.x0 + tight.x1) / 2 - width / DPI / 2, (tight.y0 + tight.y1) / 2 - height / DPI / 2,
            (width + 0.25) / DPI, (height + 0.25) / DPI,
        )
        scratch = png.with_name(png.name + ".part")
        figure.savefig(scratch, format="png", dpi=DPI, bbox_inches=window, facecolor=GROUND, metadata={"Software": None, STAMP_KEY: digest})
    written = scratch.read_bytes()
    if png_size(written) != CANVAS or png_text(written).get(STAMP_KEY) != digest:
        scratch.unlink()
        fail(f"the drawn PNG is {png_size(written)[0]} by {png_size(written)[1]} pixels or lacks its stamp; the canvas is {width} by {height}")
    scratch.replace(png)


def under_assets(path: Path) -> bool:
    """Whether `path` lies under a `docs/assets/` directory."""
    parents = path.resolve().parents
    return any(parent.name == "assets" and parent.parent.name == "docs" for parent in parents)


def main(argv: list[str] | None = None) -> None:
    here = Path(__file__).resolve()
    ap = argparse.ArgumentParser(description="Render the release's comparison figure and its values file from the provenance records.")
    ap.add_argument("--repo", type=Path, default=here.parents[2], help="repository root")
    ap.add_argument("--records", type=Path, nargs="+", help="record files or directories; default: the README's release under docs/concordance-provenance")
    ap.add_argument("--registry", type=Path, default=here.with_name("comparators.tsv"), help="the tools' registry")
    ap.add_argument("--values", type=Path, help=f"values file; default {VALUES} under --repo")
    ap.add_argument("--png", type=Path, help=f"image; default {IMAGE} under --repo")
    ap.add_argument("--values-only", action="store_true", help="write the values file and draw nothing")
    ap.add_argument("--draft", action="store_true", help=DRAFT_HELP)
    args = ap.parse_args(argv)
    values_path = args.values or args.repo / VALUES
    png_path = args.png or args.repo / IMAGE
    outputs = [values_path] if args.values_only else [values_path, png_path]
    if args.draft:
        for path in outputs:
            if under_assets(path):
                fail(f"--draft writes no file under docs/assets/: {path}")

    tools = read_registry(args.registry)
    chosen = select(discover(args.repo, args.records), tools, args.draft)
    document, printed = build_values(tools, chosen, args.draft)
    body = values_bytes(document)
    digest = hashlib.sha256(body).hexdigest()
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    if not args.values_only:
        draw(document, tools, png_path, digest)
    values_path.write_bytes(body)
    for note in printed:
        print(note)
    print(f"wrote {values_path} ({len(body)} bytes, sha256 {digest})")
    if not args.values_only:
        print(f"wrote {png_path} ({png_path.stat().st_size} bytes, {CANVAS[0]} by {CANVAS[1]} pixels)")


if __name__ == "__main__":
    main()
