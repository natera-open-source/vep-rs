#!/usr/bin/env python3
r"""Plugin tuple counts of scored engines against one Ensembl VEP default-format output.

A plugin tuple is one populated plugin field of one output row: (Location, Allele, Feature,
field name, value). A field is populated when its value is neither empty nor ``-``, so a row
where neither side writes the field contributes nothing.

Sources. ``--plugins A,B`` names the plugins the scored run loaded, by their names in the
registry of ``compare_plugin_outputs.py`` or by the names the registry lists as their
aliases. A field belongs to the plugin the registry's ``plugin_of_field`` returns among
those plugins, never to a plugin the run did not load, so a field name two plugins write is
counted once, for the run's own. ``--run-command`` hands the registry the words of the
scored run's command line, and the registry refuses a run whose flags make the name of a
loaded plugin's column a key of Ensembl VEP's own as well, where both engines would score
that key as the plugin's. ``--custom-source NAME=PREFIX`` groups the fields of a ``--custom``
source: those starting with the prefix that no plugin of the run claims. A run names its
plugins, its custom sources or both. Every populated Extra key that no source claims and
that is not a column of Ensembl VEP's own (``CORE_COLUMNS`` of ``compare_vep_formats.py``)
is listed per input under ``unclaimed_keys`` with the number of times the input writes it.
So is a field two plugins of the run both write, where no custom source claims it.

Usage:
    python3 score_plugin_tuples.py --reference <vep_output.txt[.gz]> \
        --engine <id>=<path> [--engine <id>=<path> ...] \
        [--plugins <plugin>,<plugin>] [--run-command="<the run's command line>"] \
        [--custom-source <name>=<field prefix> ...] \
        [--engine-binary <id>=<sha256>] [--engine-commit <id>=<commit>] \
        [--class-table <classes.tsv>] [--classes <class>,<class>] \
        [--gwas-catalog <associations.tsv>] [--input <the run's input.vcf[.gz]>] \
        [--revel-table <revel.tsv.gz>] [--revel-no-match] \
        [--spliceai-table <scores.vcf.gz> ...] \
        [--json-cache <the candidate JSON cache>] \
        [--utrannotator-sweep <output or dir> ...] \
        [--class-population <class>=<n> ...] \
        [--pairing-tool <vep-plugin-pairs>] [--pairing-threads <n>] \
        [--work-dir <dir>] --report <report.json> [--residuals <dir>]

Inputs. The reference is an Ensembl VEP default-format output, plain or gzip. A scored
engine's path is a default-format output or a tuple file, told apart by the column count of
the first data line (14 or 5). A tuple file is what a comparator's normaliser writes:
tab-separated, no header, five columns: Location, Allele and Feature as Ensembl VEP writes
them, the Ensembl VEP field name, the decoded value, one line per row that writes it. A
default-format value is read with the percent-encoding of ``=`` and ``;`` undone, the form a
tuple file carries, and every Location with each coordinate written as its number
(``canonical_location``). Every input is opened and its column header and first data line
checked before any is scored.

Streaming. Each input is written as tuple lines, once per projection, and, when a class is
declared, a default-format input as its rows under the same key as well. Every file passes
through external ``sort`` under C collation with its scratch in ``--work-dir``, and the
sorted files are merged line by line: no input is held in memory. The named lines of
``--input`` are sorted by identifier there too and searched in place. Rows that share a
(Location, Allele, Feature) are never collapsed: a tuple counts as often as the rows of its
side write it and pairs with as many of the other side's. The scratch holds each input's
tuples twice and its rows once, so it takes several times an input's uncompressed size. The
merge and the readers of the rows walk each sorted file forward alone, so a file with a line
that sorts below the one before it, and a value of a default-format input or of a tuple
file with no row or tuple under its key among the sorted ones, would leave values without
a pair or a class: either ends the run in the error line, on the pairing tool's files as
on this script's own. An interrupt or a termination signal stops the sorts, empties the
scratch, removes the residual files the run wrote and ends in the error line; a repeat of
either is ignored while the run stops.

Pairing tool. ``--pairing-tool BINARY`` names the built ``tools/vep-plugin-pairs``, which
then reads the inputs, extracts their tuples by the rules above, sorts them and pairs the
values equal as text, on ``--pairing-threads`` threads or on every core, in place of this
script's own extraction, sort and merge. It asks this script for the source of every field
name it meets, so the registry is read here alone, and hands over its counts and the
remainder: the tuples of either side that it did not pair as equal text, in the order and
under the keys the merge hands them on, with the rows a class can read for them: the
reference's rows and a default-format engine's rows at every Location of the remainder,
those whose values all paired among them. The value rule, the classes, the residual files
and the report are this script's on either path, and both paths write the same report and
the same residual files.

Value rule, the same for every engine. Two values are equal when their text is equal, or
when both are decimal numbers that differ, in decimal arithmetic, by at most ``TOLERANCE``
of the larger of the two and by at most ``TOLERANCE`` outright: both, so that two small
numbers are never equal by their size alone and two whole numbers are equal only when they
are the same number. A pair equal as numbers alone is counted under
``text_only_differences``. The values the two sides write under one (Location, Allele,
Feature, field) are paired one to one, a value written more than once as often as each side
writes it.

Classes. ``--class-table TSV`` declares the documented classes, one per line, tab-separated:
the class (``order_only`` or a key of ``FIELD_CLASSES`` of ``compare_vep_formats.py``), the
plugins it belongs to joined by ``;``, and for each plugin, joined by ``;`` likewise, the
field names it reads there joined by ``,``, a ``*`` standing for any run of characters.
Further columns are carried into the report and onto the class's lines of the summary,
unread but for the fourth, and ``#`` lines are skipped. A run that names no table reads
``plugin_class_table.tsv`` beside this script; under a table with no row no class is declared
and the adjusted counts equal the raw ones. ``--classes A,B`` names the declared classes the
run applies: none for an empty value, and where a run names none every class but those
whose fourth column reads ``not met``. ``--gwas-catalog FILE`` is the associations
file the scored run's GWAS plugin read: ``GwasCatalogPairs.from_catalog`` of
``compare_vep_formats.py`` reads it once, and the classes the table declares for the plugin
``GWAS`` are handed it as ``gwas_pairs``. Without it a class that reads the catalog accepts
nothing, and a line of the summary says so where a class of ``GWAS`` is applied on a run
that loaded the plugin.
``--revel-table FILE`` is the table the scored run's REVEL plugin read, with its tabix index
beside it: ``RevelTable`` of ``compare_vep_formats.py`` reads its rows at a position as a
class asks, and the classes the table declares for the plugin ``REVEL`` are handed it as
``revel_table``. ``--revel-no-match`` says the run set the plugin's ``no_match``, under which
the table's transcript ids are not read. Without the table a class that reads it accepts
nothing, and a line of the summary says so likewise.
``--utrannotator-sweep PATH``, repeatable, names Ensembl VEP's default-format outputs of the
input under fixed Perl hash seeds, one run a file, or a directory, of which every file that
opens as such an output is taken: ``UtrannotatorSeedSweep.from_outputs`` reads them once,
and the classes the table declares for the plugin ``UTRAnnotator`` are handed the result as
``utrannotator_sweep``. Without it a class accepts no row it reads the seeded runs for, and a
line of the summary says so likewise.
``--json-cache DIR`` names the candidate JSON cache of the scored run:
``DownstreamSequences`` of ``compare_vep_formats.py`` reads its coding sequences, 3' UTRs
and codon tables as a class asks, and the classes the table declares for the plugin
``Downstream`` are handed it as ``downstream_sequences``, with the applied classes as
``declared``: a value two causes change is set aside only where the second cause's class is
applied too. Without it a class that recomputes the plugin's values from the cache accepts
nothing, and a line of the summary says so likewise. The classes the table declares for the
plugin ``Blosum62`` are handed the same ``downstream_sequences``. The classes the table
declares for the plugin ``MaxEntScan`` are handed the cache's introns, ``SpliceWindows`` of
``compare_vep_formats.py``, as ``splice_windows``: without the cache a class that reads a
row's place among its transcript's introns accepts nothing, and a line of the summary says
so likewise.
``--spliceai-table FILE``, repeatable, names a score file the scored run's SpliceAI plugin
read, its ``snv=`` file and its ``indel=`` file, each with its tabix index beside it:
``SpliceAiTable`` of ``compare_vep_formats.py`` reads their records through the index at the
variant of a row a class asks about, never a file whole, and the classes the table declares
for the plugin ``SpliceAI`` are handed it as ``spliceai_table``, with the cache of
``--json-cache``, whose transcripts name the gene of a row that writes no ``SYMBOL``.
Without the files ``spliceai_multi_entry_cell`` accepts nothing, and a line of the summary
says so where the class is applied on a run that loaded the plugin; a row whose gene neither
its ``SYMBOL`` nor the cache names is not set aside. Under the files the last lines of the
summary count, per engine, the rows the class was asked about and the rows it refused, and
the refused rows by the ground of the refusal (``SPLICEAI_GROUNDS`` of
``compare_vep_formats.py``): counts alone, a row once under the first ground met for it.

Adjusted counts. A tuple left unpaired is set aside on both sides when an applied class of
its own plugin reads its field and the class's predicate accepts it: ``same_items`` for
``order_only``, ``field_class_accepts`` for any other, the first to accept in table order
credited under ``masked``. A predicate reads the row, so it is consulted only with a
reference row under the (Location, Allele, Feature) that wrote the reference's value, or
wrote no value of the field where the reference has none to pair, and an engine row that did
the same on its side; a tuple file states no rows, so there the engine's side of the row is
the tuples it wrote under that key. A side that wrote no value of the field is handed to a
predicate as absent, with one exception: the classes of ``TRANSCRIPT_READING_CLASSES`` of
``compare_vep_formats.py`` (the two Downstream classes and
``plugin_hgvs_shifted_cds_coordinates``) read an empty value as a value, the plugin's own
where the first codon it reads is a stop, so they are handed the empty text where the
side's row writes the field with an empty value. Such a value is still no tuple and moves
no raw count. Rows that
share a key each stand for one unpaired tuple of the field, in the order their input writes
them: the nth unpaired copy of a value reads the nth row that wrote it, and the nth value
left with none to pair reads the nth row of the other side that wrote no value of the field.
A value left over where the other side has no such row left, as under a key two input
records share, is therefore put to no class and stays a residual. One class is decided
without a reference row, since its shape is that the reference wrote none:
``plugin_value_on_reference_skipped_record`` takes a value a default-format engine alone
wrote at a Location of a record the reference dropped, which ``ReferenceRecordScope`` and
``reference_named_locations`` of ``compare_vep_outputs.py`` tell from the two sides' rows
at the Location and the reference's first column (``dropped_by_the_reference``), with the
gVCF variant sites its rows' names and the lines of ``--input`` show left to their own class. Under a key that rows of
the reference or of a default-format engine share, the values left pair by input record,
which a row names in its ``Uploaded_variation``, before they pair in their order: a value
pairs with the other side's value whose row carries the same ``Uploaded_variation``, else
stands against that record's row of the other side that wrote no value of the field, and
the values no record places this way pair in their order with the rows still free. The
pairs equal as text or as numbers are counted before any of this and move with none of it.
The predicate is handed, as ``record_rows``, the reference's rows of the
same input record on the Feature: those at the Location under the same
``Uploaded_variation``, where it spells the record's alleles. Under an identifier the rows
state neither that they are one record's nor its alleles in order. ``--input FILE`` is the
scored run's input VCF: ``RecordRows`` of ``compare_vep_formats.py`` takes the reference's
rows at the Location and the input's lines that carry the rows' identifiers, and returns the
record's rows, with an entry of the input line, where the identifier names one line: the one
line that carries it, or of several the one that writes the row at its Location with its
allele. A line carries the first ``;``-joined element of its ID column, the name both engines
print for its rows, and none where that element is a dot, empty or ``0``. Without
``--input``, where no line carries the identifier, where none or more than
one of several lines writes the row, and where an allele repeats among the rows, which rows
of more than one record do, nothing is handed, and a class that needs the record's rows sets
nothing aside. A class of ``SIBLING_ALLELE_CLASSES`` of
``compare_vep_formats.py`` is handed as well, as ``candidate_rows``, the engine's rows of the
same input record on the Feature, found among a default-format engine's rows at the Location
by the same rule; a tuple file states no rows, so for it nothing is handed and such a class
sets nothing aside. A class of ``RECORD_REFERENCE_CLASSES`` of ``compare_vep_formats.py``
reads the reference bases a row's allele stands in for (``Entry.reference``): a row whose
``Uploaded_variation`` spells the record's alleles carries them from the parser, and a row
under an identifier takes them from the lines of ``--input`` that carry the identifier, read
by position; without ``--input``, and where no such line or two that differ give the row at
its Location, the class reads the row without its record. A difference no applied class
accepts stays in
the raw and the adjusted counts alike. A declared class that is not applied is consulted
after the applied ones: an entry in which it would set a tuple aside names it under
``not_applied`` and carries no ``adj_`` count.

Projections, both reported for every engine. ``feature`` is the tuple above. ``allele``
drops the Feature on both sides and de-duplicates: (Location, Allele, field name, value). An
allele tuple left unpaired is set aside when every unpaired feature tuple it projects from
was, at least one, and is credited to the first of their classes in table order. The
projected values pair anew, so an allele tuple can stay unpaired although its feature tuples
all paired; such a tuple stays in the raw and the adjusted counts alike.

Report. ``scores`` holds one entry per source, engine and projection, integers only, under
the count key names of the concordance provenance records: ``perl`` and ``rust`` (the
reference's and the engine's tuples), ``intersection``, ``adj_perl``, ``adj_rust``,
``adj_intersection``, ``masked`` (per applied class of the source, the keys under which it
set a field's tuple aside, zero included), ``not_applied``, ``text_only_differences``,
``residual_rows`` and ``shared_keys``, the keys of the projection, (Location, Allele,
Feature) or (Location, Allele), under which the reference or the engine writes one field of
the source more than once: in ``feature`` the keys that rows share, in ``allele`` those with
more than one value of a field. No F1 is stored: raw F1 is 2 x intersection / (perl + rust),
adjusted F1 the same over the ``adj_`` counts. ``class_table`` names the table and its
``sha256``, ``gwas_catalog`` the catalog and its ``sha256``, ``revel_table`` the REVEL table,
its ``sha256`` and whether the run set ``no_match``, ``utrannotator_sweep`` each
output read with its ``sha256``, ``run_input`` the input of
``--input`` with its ``sha256`` and its number of records, and ``classes`` lists the
declared classes with their plugins, field names, further columns and whether the run
applied them. ``not_populated`` lists each source the
reference wrote no field of. Each input is listed with its ``sha256`` and its
``unclaimed_keys``, a tuple file with ``tuples_outside_every_source``, its tuples whose
field no source owns, and an engine with the ``binary_sha256`` and the ``commit`` that
``--engine-binary`` and ``--engine-commit`` give it; ``run_command`` holds the words
``--run-command`` gave. The summary on standard output prints each entry's counts and one
line per class of it with the class's further columns, marking a class that set nothing
aside and one whose count is above the population ``--class-population CLASS=N`` gives it,
the keys the class can reach in the scored input; it names the unclaimed keys and the
sources not populated.
``json_cache`` names the directory ``--json-cache`` gave, ``spliceai_table`` the files
``--spliceai-table`` gave, by path alone. ``spliceai_refusals`` lists, under those files, the
engines whose rows ``spliceai_multi_entry_cell`` was asked about, each with the rows
``asked`` about, the rows ``refused`` and the refused rows of each of the ``grounds``.

Residuals. ``residual_rows`` counts the tuples of either side that are outside the
intersection and that no applied class set aside. ``--residuals DIR`` writes them, never
truncated, to ``DIR/<source>.<engine>.<projection>.tsv``: one line per tuple, the side
(``reference`` or ``engine``) and then the tuple's columns, five under ``feature`` and four
under ``allele``, tab-separated, in the sort order of the tuples. Every scored source has
its file, empty when nothing is left, and ``residual_rows`` is the file's line count.
``DIR`` holds no ``.tsv`` file of an earlier run, and no two scored pairs of a source and an
engine share a file name. A run that ends in the error line, its report unwritten, keeps
none of its residual files."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import mmap
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import zlib
from collections import Counter, defaultdict
from collections.abc import Callable, Iterator
from contextlib import ExitStack, closing, contextmanager, suppress
from dataclasses import dataclass
from decimal import Context, Decimal
from functools import cache
from itertools import zip_longest
from pathlib import Path

import compare_plugin_outputs
import compare_vep_formats as formats
from compare_vep_outputs import (
    ReferenceRecordScope,
    canonical_location,
    gvcf_site_locations,
    gvcf_site_locations_from_vcf,
    perl_split,
    reference_named_locations,
    replaced_reference,
    vep_rs_pos,
)

COMPONENT = "score_plugin_tuples"
DEFAULT_HEADER = (
    "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position"
    "\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra"
)
DEFAULT_COLUMNS = 14
TUPLE_COLUMNS = 5
BLANK = (b"", b"-")
TOLERANCE = Decimal("0.000001")
DECIMAL_RE = re.compile(rb"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d{1,3})?")
DECIMAL_MAX_TEXT = 64
# Wide enough for the exact difference of any two numbers DECIMAL_RE admits within
# DECIMAL_MAX_TEXT bytes, whose digits span at most 2 x (999 + 64) places.
DECIMAL_CONTEXT = Context(prec=2200)
PROJECTIONS = {
    "feature": ["Location", "Allele", "Feature", "field", "value"],
    "allele": ["Location", "Allele", "field", "value"],
}
REFERENCE_SIDE, ENGINE_SIDE = b"R", b"E"
SIDE_NAMES = {REFERENCE_SIDE: b"reference", ENGINE_SIDE: b"engine"}
IDENT_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*")
SHA256_RE = re.compile(r"[0-9a-f]{64}")
COMMIT_RE = re.compile(r"[0-9a-f]{7,64}")
RESOLVED_LIMIT = 100_000
# The registry's name of the plugin whose associations file ``--gwas-catalog`` is.
GWAS_PLUGIN = "GWAS"
# The registry's name of the plugin whose table ``--revel-table`` is.
REVEL_PLUGIN = "REVEL"
# The registry's name of the plugin whose classes read the seeded runs of ``--utrannotator-sweep``.
UTRANNOTATOR_PLUGIN = "UTRAnnotator"
# The registry's name of the plugin whose classes read the candidate cache of ``--json-cache``.
DOWNSTREAM_PLUGIN = "Downstream"
# The registry's name of the plugin whose class reads that cache's transcripts as the classes
# of ``Downstream`` do.
BLOSUM62_PLUGIN = "Blosum62"
# The registry's name of the plugin whose class reads that cache's introns.
MAXENTSCAN_PLUGIN = "MaxEntScan"
# The registry's name of the plugin whose score files ``--spliceai-table`` names.
SPLICEAI_PLUGIN = "SpliceAI"
# What ``Sources.of`` returns for a key it leaves unclaimed: falsy as ``None`` is, and told
# apart from it by identity.
UNCLAIMED = ""
# What opens the one line the pairing tool ends a refusal in.
PAIRING_ERROR = "ERROR: [vep-plugin-pairs] "
# What the fourth column of a class table's row reads where a command line applies the
# class by its name in ``--classes`` alone.
NOT_MET = "not met"
# The class table of a command line that names none.
DEFAULT_CLASS_TABLE = Path(__file__).with_name("plugin_class_table.tsv")


class ScorerError(Exception):
    """A refusal, printed as one ``ERROR: [score_plugin_tuples]`` line."""


def text(value: bytes) -> str:
    return value.decode("utf-8", "replace")


def without_chr(chromosome: str) -> str:
    """A chromosome's name without a ``chr`` prefix of any case, the form in which an input
    record's chromosome is held to a row's."""
    return chromosome[3:] if chromosome.lower().startswith("chr") else chromosome


class Sources:
    """The sources of one scored run: the plugins it loaded and its custom sources. A field
    belongs to the plugin that ``plugin_of_field`` of ``compare_plugin_outputs.py`` returns
    among the run's plugins, else to the first custom source whose prefix it starts with:
    one source at most, and never a plugin the run did not load. ``run_words``, the words of
    the run's command line, are handed to ``plugin_of_field`` as its ``flags``."""

    def __init__(self, plugins: list[str], custom: dict[str, str], run_words: list[str] | None = None) -> None:
        self.registry = list(compare_plugin_outputs.ALL_PLUGIN_NAMES)
        plugins = [compare_plugin_outputs.PLUGIN_ALIASES.get(name, name) for name in plugins]
        unknown = [name for name in plugins if name not in self.registry]
        if unknown:
            raise ScorerError(f"--plugins names {', '.join(unknown)}, no plugin of the registry")
        for name in custom:
            if name in self.registry:
                raise ScorerError(f"--custom-source {name} names a plugin of the registry")
        if not plugins and not custom:
            raise ScorerError("no source to score: name the run's plugins (--plugins) or a --custom-source")
        self.plugins = tuple(dict.fromkeys(plugins))
        self.custom = dict(custom)
        self.names = [*self.plugins, *self.custom]
        self._flags = {} if run_words is None else {"flags": tuple(run_words)}
        self._resolved: dict[bytes, str | None] = {}

    def of(self, field: bytes) -> str | None:
        """The field's source; ``None`` for a key of Ensembl VEP's own (``CORE_COLUMNS``) and
        ``UNCLAIMED`` for any other key no source claims. A field the registry refuses for
        the run's command line is a refusal."""
        try:
            return self._resolved[field]
        except KeyError:
            name = text(field)
            try:
                source = (
                    compare_plugin_outputs.plugin_of_field(name, self.plugins, **self._flags) if self.plugins else None
                )
            except ValueError as error:
                raise ScorerError(str(error)) from error
            if source is None:
                source = next((custom for custom, prefix in self.custom.items() if name.startswith(prefix)), None)
            if source is None and name.lower() not in formats.CORE_COLUMNS:
                source = UNCLAIMED
            if len(self._resolved) >= RESOLVED_LIMIT:
                self._resolved.clear()
            self._resolved[field] = source
            return source


@cache
def field_shape(pattern: str) -> re.Pattern[str]:
    """A class table's field pattern as a regular expression, ``*`` its only wildcard."""
    return re.compile(".*".join(map(re.escape, pattern.split("*"))))


class Classes:
    """The documented classes a run declares, in the order they are consulted, each with the
    plugins it belongs to and the field patterns it reads on each, and the ones the run applies:
    those ``applied`` names, every one where it is ``None``, and with ``met_alone`` every one
    but the rows whose fourth column reads ``not met``."""

    def __init__(
        self,
        table: Path | None,
        applied: list[str] | None,
        populations: dict[str, int],
        sources: Sources,
        met_alone: bool = False,
    ) -> None:
        self.declared: dict[str, dict[str, tuple[str, ...]]] = {}
        self.further: dict[str, list[str]] = {}
        self.table: dict[str, str] | None = None
        self._reading: dict[tuple[str, str], list[str]] = {}
        known = {*sources.registry, *sources.custom}
        lines: list[str] = []
        if table is not None:
            written = table.read_bytes()
            self.table = {"path": str(table), "sha256": hashlib.sha256(written).hexdigest()}
            try:
                lines = written.decode("utf-8").splitlines()
            except UnicodeDecodeError as error:
                raise ScorerError(f"{table}: the class table is not UTF-8 text") from error
        for number, line in enumerate(lines, 1):
            if not line.strip() or line.startswith("#"):
                continue
            columns = [column.strip() for column in line.split("\t")]
            if len(columns) < 3:
                raise ScorerError(f"{table}: line {number} lacks a class, its plugins and their field names")
            name, plugins, groups = columns[0], columns[1].split(";"), columns[2].split(";")
            self.further[name] = columns[3:]
            if name != formats.ORDER_ONLY and name not in formats.FIELD_CLASSES:
                raise ScorerError(f"{table}: line {number} names {name}, no class of compare_vep_formats.py")
            if name in self.declared:
                raise ScorerError(f"{table}: line {number} declares {name} a second time")
            if len(groups) != len(plugins):
                raise ScorerError(f"{table}: line {number} gives {len(plugins)} plugins and {len(groups)} field lists")
            self.declared[name] = {}
            for plugin, group in zip(plugins, groups):
                patterns = tuple(pattern.strip() for pattern in group.split(",") if pattern.strip())
                if plugin.strip() not in known or not patterns:
                    raise ScorerError(f"{table}: line {number} names {plugin.strip()!r}, no source, or no field of it")
                if plugin.strip() in self.declared[name]:
                    raise ScorerError(f"{table}: line {number} names {plugin.strip()!r} twice")
                self.declared[name][plugin.strip()] = patterns
        for option, names in (("--classes", applied or ()), ("--class-population", populations)):
            for name in names:
                if name not in self.declared:
                    raise ScorerError(f"{option} names {name}, a class the class table does not declare")
        if applied is None:
            applied = [name for name in self.declared if not (met_alone and self.further[name][:1] == [NOT_MET])]
        self.applied = frozenset(applied)
        self.populations = populations
        self.order = {name: position for position, name in enumerate(self.declared)}

    def of(self, source: str) -> list[str]:
        """The applied classes of ``source``, in table order."""
        return [name for name, plugins in self.declared.items() if source in plugins and name in self.applied]

    def of_plugin(self, plugin: str) -> frozenset[str]:
        """The declared classes the table gives ``plugin``."""
        return frozenset(name for name, plugins in self.declared.items() if plugin in plugins)

    def reading(self, source: str, field: str) -> list[str]:
        """The classes of ``source`` that read ``field``: the applied ones in table order, then
        the others."""
        found = self._reading.get((source, field))
        if found is None:
            names = [
                name
                for name, plugins in self.declared.items()
                if any(field_shape(pattern).fullmatch(field) for pattern in plugins.get(source, ()))
            ]
            found = [name for name in names if name in self.applied] + [
                name for name in names if name not in self.applied
            ]
            if len(self._reading) >= RESOLVED_LIMIT:
                self._reading.clear()
            self._reading[source, field] = found
        return found


class _Digesting(io.RawIOBase):
    """A raw reader that digests every byte read through it."""

    def __init__(self, raw: io.FileIO) -> None:
        self.raw = raw
        self.digest = hashlib.sha256()

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int | None:
        count = self.raw.readinto(buffer)
        if count:
            self.digest.update(memoryview(buffer)[:count])
        return count


class Reader:
    """The lines of a plain or gzip file as bytes, and the sha256 of the file as stored."""

    def __init__(self, path: Path) -> None:
        self._file = path.open("rb", buffering=0)
        self._digesting = _Digesting(self._file)
        self._buffered = io.BufferedReader(self._digesting, buffer_size=1 << 20)
        gzipped = self._buffered.peek(2)[:2] == b"\x1f\x8b"
        self.lines = gzip.GzipFile(fileobj=self._buffered) if gzipped else self._buffered

    def sha256(self) -> str:
        while self._buffered.read(1 << 20):
            pass
        return self._digesting.digest.hexdigest()

    def close(self) -> None:
        self._file.close()


class Sink:
    """One ``sort`` process: lines written here reach ``path`` sorted."""

    def __init__(self, command: list[str], env: dict[str, str], path: Path) -> None:
        self.path = path
        self._out = path.open("wb")
        self._process = subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=self._out, stderr=subprocess.PIPE, env=env, bufsize=1 << 20
        )
        assert self._process.stdin is not None
        self.write = self._process.stdin.write

    def close(self) -> None:
        """Ends the input and waits for the sorted file; a sort that failed is a refusal."""
        assert self._process.stdin is not None and self._process.stderr is not None
        with suppress(BrokenPipeError):
            self._process.stdin.close()
        message = " ".join(text(self._process.stderr.read()).split())
        self._process.stderr.close()
        code = self._process.wait()
        self._out.close()
        if code != 0:
            raise ScorerError(f"sort exited {code} writing {self.path.name}: {message}")

    def kill(self) -> None:
        self._process.kill()

    def reap(self) -> None:
        """Waits for a killed sort and closes its streams, whatever ``close`` did before."""
        assert self._process.stdin is not None and self._process.stderr is not None
        with suppress(OSError):
            self._process.stdin.close()
        self._process.stderr.close()
        self._process.wait()
        self._out.close()


class Sorter:
    """Starts ``sort`` under C collation, the order the merge compares bytes in."""

    def __init__(self, work: Path, memory: str, threads: int | None) -> None:
        executable = shutil.which("sort")
        if executable is None:
            raise ScorerError("no `sort` on PATH")
        self._command = [executable, "-S", memory, "-T", str(work)]
        if threads is not None:
            self._command.append(f"--parallel={threads}")
        self._env = {**os.environ, "LC_ALL": "C"}

    @contextmanager
    def sinks(self, *paths: Path, repeats: Path | None = None) -> Iterator[list[Sink]]:
        """Sinks for ``paths``, closed on leaving; a sort that fails is reported in place of
        the broken pipe its failure raises in the writer. Anything else that ends the block
        or the wait for a sort stops every sort. The file at ``repeats`` keeps a line as
        often as it is written; every other holds it once."""
        opened: list[Sink] = []
        failure: ScorerError | None = None
        try:
            for path in paths:
                opened.append(Sink([*self._command, *([] if path == repeats else ["-u"])], self._env, path))
            try:
                yield opened
            except BrokenPipeError:
                failure = ScorerError("sort stopped reading its input")
            for sink in opened:
                try:
                    sink.close()
                except ScorerError as error:
                    failure = error
        except BaseException:
            # Every sort is killed before any is waited for, so a second interrupt during
            # the wait leaves none running.
            for step in (Sink.kill, Sink.reap):
                for sink in opened:
                    step(sink)
            raise
        if failure is not None:
            raise failure


@dataclass
class Extracted:
    """One input after extraction: its sorted tuple files and, for a default-format input
    of a run that keeps rows, its sorted rows (``rows`` is the feature tuple file of a
    tuple-file input)."""

    ident: str
    path: Path
    kind: str
    sha256: str
    lines: int
    outside: int
    unclaimed: dict[str, int]
    header: str
    feature: Path
    allele: Path
    rows: Path


def decoded(value: bytes) -> bytes:
    return formats.decode_extra(value.decode("utf-8", "surrogateescape")).encode("utf-8", "surrogateescape")


def canonical(location: bytes) -> bytes:
    """The Location with each coordinate as its number; one whose digits ``int`` refuses
    stays as written."""
    written = location.strip().decode("utf-8", "surrogateescape")
    try:
        return canonical_location(written).encode("utf-8", "surrogateescape")
    except ValueError:
        return location.strip()


def listed_keys(unclaimed: dict[bytes, int]) -> dict[str, int]:
    """The keys no source claims as text, in their bytes' order, each with the times its
    input writes it. Keys of different bytes can decode to the same text, where their
    counts add."""
    listed: dict[str, int] = {}
    for name, times in sorted(unclaimed.items()):
        shown = text(name)
        listed[shown] = listed.get(shown, 0) + times
    return listed


def first_line_kind(columns: int, headed: bool, tuple_file_allowed: bool) -> str:
    """What the column count of a file's first data line makes of the file."""
    return "tuples" if columns == TUPLE_COLUMNS and tuple_file_allowed and not headed else "default"


def wrong_width(path: Path, number: int, columns: int, tuple_file_allowed: bool) -> ScorerError:
    expected = "a default-format row has 14" + (", a tuple file line 5" if tuple_file_allowed else "")
    return ScorerError(f"{path}: line {number} has {columns} columns; {expected}")


def check_header(path: Path, number: int, row: bytes) -> None:
    """Refuses a column header, a line opening with one ``#``, that names other than
    ``DEFAULT_COLUMNS`` columns."""
    named = row.count(b"\t") + 1
    if not row.startswith(b"##") and named != DEFAULT_COLUMNS:
        raise ScorerError(f"{path}: line {number}, the column header, names {named} columns, not 14")


def check_first_line(path: Path, tuple_file_allowed: bool) -> None:
    """Refuses an input that cannot be opened, whose column header names other than
    ``DEFAULT_COLUMNS`` columns, or whose first data line is neither a default-format row nor,
    where one is allowed, a tuple file line."""
    reader = Reader(path)
    try:
        headed = False
        for number, line in enumerate(reader.lines, 1):
            row = line.rstrip(b"\r\n")
            if row.startswith(b"#"):
                headed = True
                check_header(path, number, row)
            elif row:
                columns = row.count(b"\t") + 1
                kind = first_line_kind(columns, headed, tuple_file_allowed)
                if columns != (DEFAULT_COLUMNS if kind == "default" else TUPLE_COLUMNS):
                    raise wrong_width(path, number, columns, tuple_file_allowed)
                return
    except (EOFError, zlib.error, gzip.BadGzipFile) as error:
        raise ScorerError(f"{path}: cannot read the gzip stream: {error}") from error
    finally:
        reader.close()


def extract(
    ident: str,
    path: Path,
    stem: str,
    sources: Sources,
    sorter: Sorter,
    work: Path,
    tuple_file_allowed: bool,
    keep_rows: bool,
) -> Extracted:
    """Writes the input's tuples through ``sort`` and, with ``keep_rows``, its rows. Feature
    lines are ``Location, Allele, Feature, field, value``, allele lines the same without the
    Feature, row lines the key, the row's ordinal and the row as written, all tab-separated."""
    feature_path, allele_path, rows_path = (work / f"{stem}.{part}" for part in ("feature", "allele", "rows"))
    kind: str | None = None
    headed = False
    header = DEFAULT_HEADER.encode()
    count = outside = 0
    unclaimed: dict[bytes, int] = {}
    last_location = location = b""
    sources_of = sources.of

    def unclaimed_key(name: bytes) -> None:
        unclaimed[name] = unclaimed.get(name, 0) + 1
        if len(unclaimed) > RESOLVED_LIMIT:
            raise ScorerError(f"{path}: more than {RESOLVED_LIMIT} distinct keys that no source claims")

    reader = Reader(path)
    try:
        kept = [rows_path] if keep_rows else []
        with sorter.sinks(feature_path, allele_path, *kept, repeats=feature_path) as (feature, allele, *rows):
            for number, line in enumerate(reader.lines, 1):
                row = line.rstrip(b"\r\n")
                if not row:
                    continue
                if kind != "tuples" and row.startswith(b"#"):
                    headed = True
                    check_header(path, number, row)
                    if not row.startswith(b"##"):
                        header = row
                    continue
                columns = row.split(b"\t")
                if kind is None:
                    kind = first_line_kind(len(columns), headed, tuple_file_allowed)
                if len(columns) != (DEFAULT_COLUMNS if kind == "default" else TUPLE_COLUMNS):
                    raise wrong_width(path, number, len(columns), tuple_file_allowed)
                if columns[1 if kind == "default" else 0] != last_location:
                    last_location = columns[1 if kind == "default" else 0]
                    location = canonical(last_location)
                if kind == "tuples":
                    _, allele_text, feature_text, name, value = columns
                    count += 1
                    if value in BLANK:
                        continue
                    source = sources_of(name)
                    if not source:
                        outside += 1
                        if source is not None:
                            unclaimed_key(name)
                        continue
                    allele_key = b"%b\t%b\t" % (location, allele_text.strip())
                    feature.write(b"%b%b\t%b\t%b\n" % (allele_key, feature_text.strip(), name, value))
                    allele.write(b"%b%b\t%b\n" % (allele_key, name, value))
                    continue
                allele_key = b"%b\t%b\t" % (location, columns[2].strip())
                key = b"%b%b\t" % (allele_key, columns[4].strip())
                if rows:
                    rows[0].write(b"%b%012d\t%b\n" % (key, count, row))
                count += 1
                extra = columns[13]
                if extra == b"-":
                    continue
                for pair in extra.split(b";"):
                    name, separator, value = pair.partition(b"=")
                    if not separator or value in BLANK:
                        continue
                    source = sources_of(name)
                    if not source:
                        if source is not None:
                            unclaimed_key(name)
                        continue
                    if b"%3" in value:
                        value = decoded(value)
                    feature.write(b"%b%b\t%b\n" % (key, name, value))
                    allele.write(b"%b%b\t%b\n" % (allele_key, name, value))
        digest = reader.sha256()
    except (EOFError, zlib.error, gzip.BadGzipFile) as error:
        raise ScorerError(f"{path}: cannot read the gzip stream: {error}") from error
    finally:
        reader.close()
    if kind is None:
        kind = "default" if headed or not tuple_file_allowed else "tuples"
    return Extracted(
        ident=ident,
        path=path,
        kind=kind,
        sha256=digest,
        lines=count,
        outside=outside,
        unclaimed=listed_keys(unclaimed),
        header=text(header),
        feature=feature_path,
        allele=allele_path,
        rows=rows_path if kind == "default" else feature_path,
    )


class Pairing:
    """One run of the pairing tool over the reference and every engine: what it hands over
    in place of ``extract`` and of the pairs ``merge`` counts. The tool asks, a line on its
    standard output each, for the source of a field name (``field <hex>``) and for the
    canonical form of a Location it does not settle itself (``location <hex>``), and is
    answered on its standard input from ``sources`` and ``canonical``. It leaves
    ``pairs.json`` and, per engine and projection, each side's tuples that did not pair as
    equal text, sorted as ``sort`` sorts them, with each side's rows at their Locations where
    ``keep_rows`` holds. A refusal of the tool or of the registry is this run's."""

    def __init__(
        self,
        tool: Path,
        reference: Path,
        engines: list[tuple[str, Path]],
        sources: Sources,
        work: Path,
        keep_rows: bool,
        threads: int | None,
    ) -> None:
        self._directory = work / "pairs"
        self._inputs = [("reference", reference), *engines]
        command = [str(tool), "--reference", str(reference)]
        for _, path in engines:
            command += ["--engine", str(path)]
        command += ["--out", str(self._directory), "--distinct-keys", str(RESOLVED_LIMIT)]
        if keep_rows:
            command.append("--rows")
        if threads is not None:
            command += ["--threads", str(threads)]
        refused: dict[bytes, ScorerError] = {}
        with (work / "pairs.stderr").open("w+b") as errors:
            try:
                # Run by its path: a file named without a directory would be looked up on PATH.
                tool_run = subprocess.Popen(
                    command, executable=tool.absolute(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=errors
                )
            except OSError as error:
                raise ScorerError(f"cannot run the pairing tool {tool}: {error.strerror or error}") from error
            assert tool_run.stdin is not None and tool_run.stdout is not None
            try:
                for question in tool_run.stdout:
                    asked, _, subject = question.strip().partition(b" ")
                    subject = bytes.fromhex(subject.decode("ascii"))
                    if asked == b"location":
                        answer = b"location " + canonical(subject).hex().encode()
                    else:
                        try:
                            source = sources.of(subject)
                        except ScorerError as error:
                            refused[subject] = error
                            answer = b"refused"
                        else:
                            claimed = b"source " + source.encode() if source else b"unclaimed"
                            answer = b"core" if source is None else claimed
                    try:
                        tool_run.stdin.write(answer + b"\n")
                        tool_run.stdin.flush()
                    except BrokenPipeError:
                        break
                code = tool_run.wait()
            except BaseException:
                tool_run.kill()
                tool_run.wait()
                raise
            finally:
                with suppress(OSError):
                    tool_run.stdin.close()
                tool_run.stdout.close()
            errors.seek(0)
            written = text(errors.read())
        lines = written.splitlines()
        if code != 0:
            # The error line is the tool's last, and the path it names can carry a line boundary.
            _, found, told = f"\n{written}".rpartition(f"\n{PAIRING_ERROR}")
            if code != 2 or not found:
                raise ScorerError(f"the pairing tool exited {code}: {' '.join(' '.join(lines).split())}")
            told = told.removesuffix("\n")
            for name, refusal in refused.items():
                if told == f"the registry refuses the field {name.hex()}":
                    raise refusal
            # The tool writes a path's bytes that are no text as replacement characters.
            for _, path in self._inputs:
                shown = os.fsencode(path).decode("utf-8", "replace")
                if told.startswith(f"{shown}: "):
                    told = f"{path}{told[len(shown) :]}"
                    break
            raise ScorerError(told)
        # The tool is silent unless its environment asks for the time of each pass.
        for line in lines:
            print(line, file=sys.stderr)
        self._stated = json.loads((self._directory / "pairs.json").read_text(encoding="utf-8"))
        self.reference = self._input(0, 1, "reference")

    def _input(self, at: int, engine: int, side: str) -> Extracted:
        """The run's input ``at``, the reference being 0, with the files of its ``side`` of
        the remainder against the engine that is input ``engine``."""
        stated = self._stated["inputs"][at]
        ident, path = self._inputs[at]
        files = {part: self._directory / f"input{engine}.{part}.{side}" for part in ("feature", "allele", "rows")}
        return Extracted(
            ident=ident,
            path=path,
            kind=stated["kind"],
            sha256=stated["sha256"],
            lines=stated["lines"],
            outside=stated["outside"],
            unclaimed=listed_keys({bytes.fromhex(name): times for name, times in stated["unclaimed"].items()}),
            header=text(bytes.fromhex(stated["header"])),
            **files,
        )

    def engine(self, index: int) -> tuple[Extracted, Extracted, dict]:
        """The reference and the engine that is input ``index``, each as its side of their
        remainder, and the tool's counts of the pairing per projection."""
        return self._input(0, index, "reference"), self._input(index, index, "engine"), self._stated["engines"][index - 1]


def as_decimal(value: bytes) -> Decimal | None:
    """The value as a decimal number: an optional sign, digits with an optional fraction and
    an optional exponent of up to three digits, within ``DECIMAL_MAX_TEXT`` bytes."""
    if len(value) > DECIMAL_MAX_TEXT or DECIMAL_RE.fullmatch(value) is None:
        return None
    return Decimal(value.decode("ascii"))


def numbers_agree(one: bytes, other: bytes) -> bool:
    """Whether both are decimal numbers that differ by at most ``TOLERANCE`` of the larger of
    the two in size and by at most ``TOLERANCE`` outright."""
    first, second = as_decimal(one), as_decimal(other)
    if first is None or second is None:
        return False
    apart = DECIMAL_CONTEXT.subtract(first, second).copy_abs()
    larger = max(first.copy_abs(), second.copy_abs())
    return apart <= TOLERANCE and apart <= DECIMAL_CONTEXT.multiply(TOLERANCE, larger)


def left_over(values: list[bytes], paired: Counter[bytes]) -> list[bytes]:
    """``values`` in their order without the ones ``paired`` counts, each as often as counted."""
    if not paired:
        return values
    owed = dict(paired)
    left = []
    for value in values:
        if owed.get(value):
            owed[value] -= 1
        else:
            left.append(value)
    return left


def pair_values(ref_values: list[bytes], eng_values: list[bytes]) -> tuple[int, int, list[bytes], list[bytes]]:
    """Pairs two sides' values one to one, equal text first, then as many equal numbers as
    the two sides allow, a value as often as its side holds it: the pairs equal as text, the
    pairs equal as numbers alone, and the values each side has left."""
    as_text = Counter(ref_values) & Counter(eng_values)
    ref_left, eng_left = left_over(ref_values, as_text), left_over(eng_values, as_text)
    ref_numeric: Counter[bytes] = Counter()
    eng_numeric: Counter[bytes] = Counter()
    as_numbers = 0
    if ref_left and eng_left:
        # The numbers one number agrees with are a run of the other side's order that moves up
        # with it, so walking both orders and dropping the smaller of two that do not agree
        # pairs the most numbers the rule allows.
        ref_numbers = sorted((number, value) for value in ref_left if (number := as_decimal(value)) is not None)
        eng_numbers = sorted((number, value) for value in eng_left if (number := as_decimal(value)) is not None)
        i = j = 0
        while i < len(ref_numbers) and j < len(eng_numbers):
            if numbers_agree(ref_numbers[i][1], eng_numbers[j][1]):
                ref_numeric[ref_numbers[i][1]] += 1
                eng_numeric[eng_numbers[j][1]] += 1
                as_numbers += 1
                i += 1
                j += 1
            elif ref_numbers[i][0] < eng_numbers[j][0]:
                i += 1
            else:
                j += 1
        ref_left, eng_left = left_over(ref_left, ref_numeric), left_over(eng_left, eng_numeric)
    return as_text.total(), as_numbers, ref_left, eng_left


class Tally:
    """Per field name: the tuples paired, those paired as numbers alone and those each side
    has left. Per source: the tuples of each side set aside, per class the keys under which
    it set a tuple aside (``masked``) or would if it were applied (``not_applied``), and the
    keys under which a side writes one field more than once (``shared``)."""

    def __init__(self) -> None:
        self.paired: defaultdict[bytes, int] = defaultdict(int)
        self.text_only: defaultdict[bytes, int] = defaultdict(int)
        self.ref_left: defaultdict[bytes, int] = defaultdict(int)
        self.eng_left: defaultdict[bytes, int] = defaultdict(int)
        self.ref_aside: defaultdict[str, int] = defaultdict(int)
        self.eng_aside: defaultdict[str, int] = defaultdict(int)
        self.masked: defaultdict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        self.not_applied: defaultdict[str, defaultdict[str, int]] = defaultdict(lambda: defaultdict(int))
        self.shared: defaultdict[str, int] = defaultdict(int)

    def fields(self) -> set[bytes]:
        return set(self.paired) | set(self.ref_left) | set(self.eng_left)

    def pair(self, name: bytes, ref_values: list[bytes], eng_values: list[bytes]) -> tuple[list[bytes], list[bytes]]:
        as_text, as_numbers, ref_left, eng_left = pair_values(ref_values, eng_values)
        self.paired[name] += as_text + as_numbers
        self.text_only[name] += as_numbers
        self.ref_left[name] += len(ref_left)
        self.eng_left[name] += len(eng_left)
        return ref_left, eng_left

    def take(self, counted: dict) -> None:
        """Adds what the pairing tool counted of one projection: per field the pairs equal
        as text, and per source the keys under which a side writes one field more than
        once. The tool states each side's values of a field as well: its pairs and the
        values ``pair`` was handed of the field, or its counts and its remainder disagree."""
        stated = {bytes.fromhex(name): counts for name, counts in counted["values"].items()}
        for field in self.fields() | set(stated):
            reference, engine, paired = stated.get(field, (0, 0, 0))
            settled = self.paired.get(field, 0)
            handed = (settled + self.ref_left.get(field, 0), settled + self.eng_left.get(field, 0))
            if handed != (reference - paired, engine - paired):
                raise ScorerError(f"the pairing tool's counts of {text(field)} and its remainder disagree")
            self.paired[field] += paired
        for source, keys in counted["shared_keys"].items():
            self.shared[source] += keys


def ascending(path: Path) -> Iterator[bytes]:
    """The lines of a sorted file of the scratch, each with its newline. A line that sorts
    below the one before it is a refusal: the merge and the cursors read forward alone, so
    past such a line they would miss the lines it belongs with and leave their values
    without a pair or a class. ``sort`` under C collation compares lines without their
    newline, so a line that opens the next one sorts first whatever byte follows it there."""
    previous = b""
    with path.open("rb") as handle:
        for number, line in enumerate(handle, 1):
            # The bytes compare as ``sort`` compares them unless one line opens the other.
            suspect = line < previous or (len(line) < len(previous) and previous[len(line) - 1] < 0x0A)
            if suspect and line.removesuffix(b"\n") < previous.removesuffix(b"\n"):
                raise ScorerError(
                    f"{path.name}: line {number} sorts below line {number - 1}: the file is not in ascending order"
                )
            yield line
            previous = line


def merge(
    reference: Path,
    engine: Path,
    paired: defaultdict[bytes, int],
    settle: Callable[[bytes, list[bytes], list[bytes]], None],
    shared: Callable[[bytes], None],
) -> None:
    """Walks two sorted tuple files together. A line both carry is a pair, counted under its
    field name, the column before the value; any other run of lines sharing the columns
    before the value goes to ``settle`` as that prefix and the two sides' values. ``shared``
    is called with those columns wherever a side carries a second line under them. A file
    that is not in ascending order is a refusal (``ascending``)."""
    with closing(ascending(reference)) as ref, closing(ascending(engine)) as eng:
        r, e = next(ref, b""), next(eng, b"")
        last_head = None
        while r or e:
            if r == e:
                head = r[: r.rindex(b"\t")]
                paired[head[head.rindex(b"\t") + 1 :]] += 1
                if head == last_head:
                    shared(head)
                last_head = head
                r, e = next(ref, b""), next(eng, b"")
                continue
            ref_prefix = r[: r.rindex(b"\t") + 1] if r else None
            eng_prefix = e[: e.rindex(b"\t") + 1] if e else None
            ref_values: list[bytes] = []
            eng_values: list[bytes] = []
            if eng_prefix is None or (ref_prefix is not None and ref_prefix <= eng_prefix):
                prefix = ref_prefix
                while r.startswith(prefix):
                    ref_values.append(r[len(prefix) : -1])
                    r = next(ref, b"")
            if ref_prefix is None or (eng_prefix is not None and eng_prefix <= ref_prefix):
                prefix = eng_prefix
                while e.startswith(prefix):
                    eng_values.append(e[len(prefix) : -1])
                    e = next(eng, b"")
            # A run under the columns of the line both sides carried just before it is a second line there.
            if len(ref_values) > 1 or len(eng_values) > 1 or prefix[:-1] == last_head:
                shared(prefix[:-1])
            settle(prefix, ref_values, eng_values)


class Cursor:
    """A forward-only reader of one sorted file; a file that is not in ascending order is a
    refusal (``ascending``)."""

    def __init__(self, path: Path) -> None:
        self._lines = ascending(path)
        self._line = next(self._lines, b"")

    def starting(self, prefix: bytes) -> list[bytes]:
        """The lines starting with ``prefix``, a tab-terminated run of whole columns, without
        their newline; successive calls take prefixes in ascending byte order."""
        line = self._line
        lines = self._lines
        while line and line < prefix:
            line = next(lines, b"")
        found = []
        while line.startswith(prefix):
            found.append(line[:-1])
            line = next(lines, b"")
        self._line = line
        return found

    def close(self) -> None:
        self._lines.close()


class InputIndex:
    """The scored run's input VCF by identifier: the first five columns of its named records,
    sorted in the scratch and searched there, none held in memory. A record carries the first
    ``;``-joined element of its ID column, the name both engines print for its rows, and no
    name where that element is a dot, empty or ``0``."""

    def __init__(self, path: Path, sorter: Sorter, work: Path) -> None:
        index = work / "input.names"
        self.records = 0
        reader = Reader(path)
        try:
            with sorter.sinks(index, repeats=index) as (sink,):
                for line in reader.lines:
                    columns = line.rstrip(b"\r\n").split(b"\t", 5)
                    if line.startswith(b"#") or len(columns) < 5:
                        continue
                    self.records += 1
                    chromosome, position, name, reference, alternates = columns[:5]
                    name = name.split(b";", 1)[0]
                    if name not in (b".", b"", b"0"):
                        sink.write(b"%b\t%b\t%b\t%b\t%b\n" % (name, chromosome, position, reference, alternates))
            self.sha256 = reader.sha256()
        except (EOFError, zlib.error, gzip.BadGzipFile) as error:
            raise ScorerError(f"{path}: cannot read the gzip stream: {error}") from error
        finally:
            reader.close()
        if not self.records:
            raise ScorerError(f"{path}: no VCF record, a data line of five or more columns")
        self._file = index.open("rb")
        # An empty file cannot be mapped: an input of unnamed records names none.
        self._lines = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ) if index.stat().st_size else None

    def named(self, name: bytes) -> list[str]:
        """The input's lines whose ID column has the identifier ``name`` as its first
        ``;``-joined element, byte for byte, each as its first five columns in the VCF's
        order with ``name`` as the third: identifiers of different bytes can decode to one
        text."""
        lines, prefix = self._lines, name + b"\t"
        if lines is None:
            return []
        # Every line that starts before ``low`` sorts below the prefix, none from ``high`` on.
        low, high = 0, len(lines)
        while low < high:
            start = lines.rfind(b"\n", 0, (low + high) // 2) + 1
            end = lines.find(b"\n", start)
            if lines[start:end] < prefix:
                low = end + 1
            else:
                high = start
        found = []
        while lines[low : low + len(prefix)] == prefix:
            end = lines.find(b"\n", low)
            chromosome, position, reference, alternates = text(lines[low + len(prefix) : end]).split("\t")
            found.append(f"{chromosome}\t{position}\t{text(name)}\t{reference}\t{alternates}")
            low = end + 1
        return found

    def close(self) -> None:
        if self._lines is not None:
            self._lines.close()
        self._file.close()


class OutputRows:
    """Rows of an output, or a stored output plain or gzip, as the file a reader of
    ``compare_vep_outputs.py`` opens: those readers take a path and call its ``exists`` and
    ``open`` alone."""

    def __init__(self, rows: list[str] | None = None, path: Path | None = None) -> None:
        self._rows = rows
        self._path = path

    def exists(self) -> bool:
        return True

    def open(self, *_mode: object, **_how: object) -> io.TextIOBase:
        if self._path is None:
            return io.StringIO("".join(row.rstrip("\n") + "\n" for row in self._rows or ()))
        with self._path.open("rb") as stored:
            gzipped = stored.read(2) == b"\x1f\x8b"
        opener = gzip.open if gzipped else open
        return opener(self._path, "rt", encoding="utf-8", errors="replace")


def dropped_by_the_reference(
    reference_rows: list[str], engine_rows: list[str], reference: Path, input_lines: list[str]
) -> bool:
    """Whether the rows the two sides wrote at one Location, each a default-format row as
    written, are those of a record the reference dropped, as the class
    ``reference_skipped_record`` of ``compare_vep_outputs.py`` reads it, by that program's own
    routines in its own order: the gVCF variant sites the engine's generated names and the
    record lines ``input_lines`` show (``gvcf_site_locations``,
    ``gvcf_site_locations_from_vcf``) are another class's; ``ReferenceRecordScope`` takes a
    Location at which the engine wrote a structural row and the reference none; and
    ``reference_named_locations`` gives it back where an engine row there names no record or
    one the reference's output ``reference`` carries on a row, which is one read of that
    output, made only for a Location the scope took."""
    group = "output"
    engine = OutputRows(engine_rows)
    sites = gvcf_site_locations(engine)
    if input_lines:
        with tempfile.NamedTemporaryFile("w", suffix=".vcf", encoding="utf-8") as records:
            records.writelines(f"{line}\n" for line in input_lines)
            records.flush()
            sites |= gvcf_site_locations_from_vcf(Path(records.name))
    scope = ReferenceRecordScope()
    scope.record_gvcf_sites(group, sites)
    locations = set()
    for rows, by_reference in ((reference_rows, True), (engine_rows, False)):
        for row in rows:
            columns = row.split("\t", 3)
            location = canonical_location(columns[1].strip())
            locations.add(location)
            scope.observe(group, f"{location}\t{columns[2].strip()}\t", by_reference, not by_reference)
    scope.close()
    absent = scope.absent_locations(group)
    scope.record_reference_named(group, reference_named_locations(OutputRows(path=reference), engine, absent))
    return bool(locations) and locations == scope.absent_locations(group)


class LocationRows:
    """A default-format input's rows at one Location: the row lines under each key and on
    each Feature, their entries, parsed as they are asked for, and the rows of each input
    record, found once per row's Feature, identifier and allele."""

    def __init__(self, header: str) -> None:
        self._header = header
        self.by_key: dict[bytes, list[bytes]] = {}
        self.by_feature: dict[bytes, list[bytes]] = {}
        self.records: object | None = None
        self.of_record: dict[tuple[bytes, bytes, bytes], list[formats.Entry] | None] = {}
        self._parsed: dict[bytes, formats.Entry] = {}

    def take(self, location: bytes, lines: list[bytes]) -> None:
        """Replaces the rows held with ``lines``, the input's rows at ``location``."""
        self.by_key, self.by_feature, self._parsed = {}, {}, {}
        self.records, self.of_record = None, {}
        for line in lines:
            _, allele, feature, _ = line.split(b"\t", 3)
            self.by_key.setdefault(b"%b%b\t%b\t" % (location, allele, feature), []).append(line)
            self.by_feature.setdefault(feature, []).append(line)

    def entry(self, line: bytes) -> formats.Entry:
        entry = self._parsed.get(line)
        if entry is None:
            row = text(line.split(b"\t", 4)[4])
            entry = self._parsed[line] = formats.parse_default(f"{self._header}\n{row}")[0]
        return entry


class RowContext:
    """The rows behind an unpaired pair of values, as the entries the predicates of
    ``compare_vep_formats.py`` read, and the classes that accept the pair."""

    def __init__(
        self,
        reference: Extracted,
        engine: Extracted,
        sources: Sources,
        classes: Classes,
        plugin_reads: dict[str, dict[str, object]],
        run_input: InputIndex | None,
    ) -> None:
        self._reference = reference
        self._engine = engine
        self._sources = sources
        self._classes = classes
        self._input = run_input
        self._handed: dict[str, dict[str, object]] = {}
        for plugin, keywords in plugin_reads.items():
            for klass in classes.of_plugin(plugin):
                self._handed.setdefault(klass, {}).update(keywords)
        self._ref_rows: Cursor | None = None
        self._eng_rows: Cursor | None = None
        self._location: bytes | None = None
        self._ref = LocationRows(reference.header)
        self._eng = LocationRows(engine.header)
        self._key: bytes | None = None
        self._eng_lines: list[bytes] = []
        self._eng_entries: list[formats.Entry] = []
        self._order_only = formats.order_only_fields("default")
        self._dropped: tuple[bytes, bool] | None = None

    def close(self) -> None:
        for cursor in (self._ref_rows, self._eng_rows):
            if cursor is not None:
                cursor.close()

    def _load(self, key: bytes) -> None:
        """Reads the reference's rows and a default-format engine's rows of the key's
        Location, once per Location, and a tuple file's tuples under the key."""
        if self._ref_rows is None or self._eng_rows is None:
            self._ref_rows, self._eng_rows = Cursor(self._reference.rows), Cursor(self._engine.rows)
        default = self._engine.kind == "default"
        location = key[: key.index(b"\t") + 1]
        if location != self._location:
            self._location = location
            self._ref.take(location, self._ref_rows.starting(location))
            if default:
                self._eng.take(location, self._eng_rows.starting(location))
        self._key = key
        if default:
            self._eng_lines = self._eng.by_key.get(key, [])
            self._eng_entries = [self._eng.entry(line) for line in self._eng_lines]
        else:
            self._eng_lines = self._eng_rows.starting(key)

    def _under(self, key: bytes, field: str, ref_wrote: bool, eng_wrote: bool) -> list[bytes]:
        """The reference's row lines under ``key``, the rows of the key loaded. A side that
        writes a value of ``field`` under the key takes it from a row there, so a side with
        a value and no row under the key is a refusal: the reader of its sorted rows misses
        the key's rows, and its values would be left without a class."""
        if key != self._key:
            self._load(key)
        lines = self._ref.by_key.get(key, [])
        for wrote, found, side in ((ref_wrote, lines, self._reference), (eng_wrote, self._eng_lines, self._engine)):
            if wrote and not found:
                location, allele, feature = (text(part) for part in key[:-1].split(b"\t"))
                raise ScorerError(
                    f"{side.path}: no row under Location {location}, Allele {allele}, Feature {feature} was read"
                    f" from the sorted rows, where the input writes a value of {field}"
                )
        return lines

    @staticmethod
    def _writing(entries: list[formats.Entry], field: str, value: str, nth: int) -> int | None:
        """The index of the row that writes ``field`` as ``value``, or for an absent value
        writes none, after ``nth`` rows that do the same; ``None`` without such a row."""
        for index, entry in enumerate(entries):
            written = entry.fields.get(field, formats.ABSENT)
            if formats.blank(written) if value == formats.ABSENT else formats.decode_extra(written) == value:
                if not nth:
                    return index
                nth -= 1
        return None

    def _record_rows(self, line: bytes, rows: LocationRows) -> list[formats.Entry] | None:
        """The rows of one row's input record among ``rows``, its side's rows at the
        Location. Under an ``Uploaded_variation`` that spells the record's alleles, those on
        the row's Feature. Under an identifier, what ``RecordRows`` of
        ``compare_vep_formats.py`` returns for the row, given the Location's rows and the
        input's lines that carry their identifiers: the rows of the one line that carries
        the identifier, or of several the one that writes the row at its Location with its
        allele, so two rows of one identifier can read two records. ``None`` where neither
        names one record, and where an allele repeats among the rows on the Feature, which
        rows of more than one record do."""
        _, allele, feature, _, uploaded, _ = line.split(b"\t", 5)
        try:
            return rows.of_record[feature, uploaded, allele]
        except KeyError:
            found = rows.of_record[feature, uploaded, allele] = self._rows_of_record(line, feature, uploaded, rows)
            return found

    def _rows_of_record(
        self, line: bytes, feature: bytes, uploaded: bytes, rows: LocationRows
    ) -> list[formats.Entry] | None:
        group = [other for other in rows.by_feature[feature] if other.split(b"\t", 5)[4] == uploaded]
        alleles = [other.split(b"\t", 2)[1] for other in group]
        if len(set(alleles)) != len(alleles):
            return None
        if formats.uploaded_alleles(text(uploaded)) is not None:
            return [rows.entry(other) for other in group]
        if self._input is None:
            return None
        if rows.records is None:
            held = [other for lines in rows.by_feature.values() for other in lines]
            names = sorted({other.split(b"\t", 5)[4] for other in held})
            carried = [carrier for name in names for carrier in self._input.named(name)]
            rows.records = formats.RecordRows([rows.entry(other) for other in held], carried)
        return rows.records.of(rows.entry(line))

    def _tuple_entry(self, field: str, value: str) -> formats.Entry | None:
        """The tuples a tuple-file engine wrote under the key, as an entry; ``None`` where it
        wrote the field and ``value`` is absent, so that no row of it wrote none."""
        assert self._key is not None
        location, allele, feature = (text(part) for part in self._key[:-1].split(b"\t"))
        fields = {"Location": location, "Allele": allele, "Feature": feature}
        for line in self._eng_lines:
            name, _, written = line[len(self._key) :].partition(b"\t")
            fields.setdefault(text(name), text(written))
        if value != formats.ABSENT:
            fields[field] = value
        elif field in fields:
            return None
        return formats.Entry(
            key=f"{location}|{allele}|{feature}|",
            record=None,
            location=location,
            allele=allele,
            feature=feature,
            consequence="",
            fields=fields,
        )

    def _rows(
        self, key: bytes, field: str, expected: str, actual: str, ref_nth: int, eng_nth: int
    ) -> tuple[formats.Entry, formats.Entry, list[formats.Entry] | None, bytes | None, bytes] | None:
        """The reference's row of the pair, the engine's, the reference's rows of the same
        record, the row line of a default-format engine's row and that of the reference's
        row; ``None`` where a side has no row under ``key`` that wrote its value of the pair,
        or wrote none where its value is absent, after the ``ref_nth`` and ``eng_nth`` rows
        that did and that belong to other pairs."""
        lines = self._under(key, field, expected != formats.ABSENT, actual != formats.ABSENT)
        entries = [self._ref.entry(line) for line in lines]
        index = self._writing(entries, field, expected, ref_nth)
        if index is None:
            return None
        engine_line = None
        if self._engine.kind == "default":
            at = self._writing(self._eng_entries, field, actual, eng_nth)
            a = None if at is None else self._eng_entries[at]
            engine_line = None if at is None else self._eng_lines[at]
        else:
            a = self._tuple_entry(field, actual)
        if a is None:
            return None
        return entries[index], a, self._record_rows(lines[index], self._ref), engine_line, lines[index]

    def _candidate_rows(self, engine_line: bytes | None) -> list[formats.Entry] | None:
        """The engine's rows of the input record of its row ``engine_line``, found as the
        reference's are (``_record_rows``); ``None`` for a tuple file, which states no rows."""
        return None if engine_line is None else self._record_rows(engine_line, self._eng)

    def _read_record(self, entry: formats.Entry, line: bytes | None) -> None:
        """Gives the row ``line`` under an identifier, as its entry, the reference bases its
        allele stands in for (``Entry.reference``): read by position from the lines of
        ``--input`` that carry the identifier on the row's chromosome, a ``chr`` prefix set
        aside, as ``read_record_references`` of ``compare_vep_outputs.py`` reads a record. A
        row whose ``Uploaded_variation`` spells the alleles carries its own from the parser.
        A row no such line states, or two of them read differently, keeps ``None``, as every
        row under an identifier does without ``--input`` and a tuple file's entry does."""
        if line is None or self._input is None or entry.reference is not None or entry.location is None:
            return
        uploaded = line.split(b"\t", 5)[4]
        if formats.uploaded_alleles(text(uploaded)) is not None:
            return
        chromosome = without_chr(entry.location.partition(":")[0])
        read: set[str] = set()
        for carrier in self._input.named(uploaded):
            carried_on, position, _, reference, alternates = carrier.split("\t")
            first = vep_rs_pos(position)
            if first is None or without_chr(carried_on) != chromosome:
                continue
            found = replaced_reference(first, reference, perl_split(alternates, ","), entry.location, entry.allele)
            if found is not None:
                read.add(found)
        if len(read) == 1:
            entry.reference = read.pop()

    def _dropped_record(self, key: bytes, field: str) -> bool:
        """Whether ``key``'s Location holds the rows of a record the reference dropped
        (``dropped_by_the_reference``), asked once per Location. A tuple file states no rows
        and its Locations are none; the record lines handed are those of ``--input`` that
        carry the first column of an engine row at the Location."""
        if self._engine.kind != "default":
            return False
        self._under(key, field, False, True)
        location = key[: key.index(b"\t") + 1]
        if self._dropped is None or self._dropped[0] != location:
            sides = [
                [text(line.split(b"\t", 4)[4]) for lines in rows.by_feature.values() for line in lines]
                for rows in (self._ref, self._eng)
            ]
            names = sorted({line.split(b"\t", 5)[4] for lines in self._eng.by_feature.values() for line in lines})
            records = [] if self._input is None else [carrier for name in names for carrier in self._input.named(name)]
            self._dropped = (location, dropped_by_the_reference(*sides, self._reference.path, records))
        return self._dropped[1]

    @staticmethod
    def _as_written(entry: formats.Entry, field: str, value: str) -> str:
        """``value`` as a class of ``TRANSCRIPT_READING_CLASSES`` is handed it: the empty text
        for a side that left no tuple of the field because its row writes the field with an
        empty value, which those classes read as a value; any other value as it is."""
        return "" if value == formats.ABSENT and entry.fields.get(field) == "" else value

    @staticmethod
    def _numbered(values: list[bytes]) -> list[tuple[bytes, int]]:
        """Each value with the number of its copies before it; the copies are adjacent."""
        numbered: list[tuple[bytes, int]] = []
        for value in values:
            nth = numbered[-1][1] + 1 if numbered and numbered[-1][0] == value else 0
            numbered.append((value, nth))
        return numbered

    def pairs(
        self, key: bytes, name: bytes, ref_left: list[bytes], eng_left: list[bytes]
    ) -> list[tuple[bytes | None, bytes | None, int, int]]:
        """The values the two sides have left of a field under ``key`` as the pairs ``judge``
        takes: a value of each side, ``None`` for a side with none to pair, and for each the
        count ``judge`` reads its row by. Where no class reads the field, the engine is a
        tuple file or no rows share the key, the values pair in their order. Under a key rows
        share they pair by input record first, the ``Uploaded_variation`` of the row each
        value reads: a value with the other side's value of the same record, else with that
        record's row of the other side that wrote no value of the field; the values no record
        places pair in their order, a side with none to pair reading its rows that wrote no
        value and are not taken. Keys come in ascending byte order."""
        ref_items, eng_items = self._numbered(ref_left), self._numbered(eng_left)
        field = text(name)
        if self._engine.kind == "default" and self._classes.reading(self._sources.of(name), field):
            ref_lines = self._under(key, field, bool(ref_items), bool(eng_items))
            if len(ref_lines) > 1 or len(self._eng_lines) > 1:
                ref_entries = [self._ref.entry(line) for line in ref_lines]
                engine = (self._eng_entries, self._eng_lines)
                return self._by_record(field, ref_items, eng_items, (ref_entries, ref_lines), engine)
        in_order: list[tuple[bytes | None, bytes | None, int, int]] = []
        for at, (ref_item, eng_item) in enumerate(zip_longest(ref_items, eng_items)):
            ref_value, ref_nth = (None, at - len(ref_items)) if ref_item is None else ref_item
            eng_value, eng_nth = (None, at - len(eng_items)) if eng_item is None else eng_item
            in_order.append((ref_value, eng_value, ref_nth, eng_nth))
        return in_order

    def _placed(
        self, field: str, items: list[tuple[bytes, int]], entries: list[formats.Entry], lines: list[bytes]
    ) -> tuple[list[bytes | None], list[bytes]]:
        """Of one side's rows under a key, as entries and as row lines: the record of the row
        each value of ``items`` reads, ``None`` without such a row, and the records of the rows
        that wrote no value of the field, in the rows' order."""
        records = [line.split(b"\t", 5)[4] for line in lines]
        reads = (self._writing(entries, field, text(value), nth) for value, nth in items)
        wrote_none = (formats.blank(entry.fields.get(field, formats.ABSENT)) for entry in entries)
        none = [record for record, blank in zip(records, wrote_none) if blank]
        return [None if at is None else records[at] for at in reads], none

    @staticmethod
    def _stood(free: list[int], records: list[bytes | None], none: list[bytes]) -> tuple[list[int], dict[int, int]]:
        """Stands each value of ``free`` against its record's first row of the other side that
        wrote no value and is not taken, ``none`` being those rows' records: the values still
        free, and per value stood the row's count among those rows."""
        still: list[int] = []
        stood: dict[int, int] = {}
        for at in free:
            taken = set(stood.values())
            row = next((row for row, record in enumerate(none) if record == records[at] and row not in taken), None)
            if records[at] is None or row is None:
                still.append(at)
            else:
                stood[at] = row
        return still, stood

    def _by_record(
        self,
        field: str,
        ref_items: list[tuple[bytes, int]],
        eng_items: list[tuple[bytes, int]],
        reference: tuple[list[formats.Entry], list[bytes]],
        engine: tuple[list[formats.Entry], list[bytes]],
    ) -> list[tuple[bytes | None, bytes | None, int, int]]:
        """The pairs of ``pairs`` under a key rows share; ``reference`` and ``engine`` are each
        side's rows under the key, as entries and as row lines."""
        ref_records, ref_none = self._placed(field, ref_items, *reference)
        eng_records, eng_none = self._placed(field, eng_items, *engine)
        found: list[tuple[bytes | None, bytes | None, int, int]] = []
        ref_free: list[int] = []
        eng_free = list(range(len(eng_items)))
        for at, (value, nth) in enumerate(ref_items):
            record = ref_records[at]
            same = next((other for other in eng_free if record is not None and eng_records[other] == record), None)
            if same is None:
                ref_free.append(at)
            else:
                eng_free.remove(same)
                found.append((value, eng_items[same][0], nth, eng_items[same][1]))
        ref_free, ref_stood = self._stood(ref_free, ref_records, eng_none)
        eng_free, eng_stood = self._stood(eng_free, eng_records, ref_none)
        found += [(ref_items[at][0], None, ref_items[at][1], row) for at, row in ref_stood.items()]
        found += [(None, eng_items[at][0], row, eng_items[at][1]) for at, row in eng_stood.items()]
        # A side with none to pair reads its next row that wrote no value and is not taken; past
        # the last of them the count names no row, and the value is put to no class.
        ref_rows = [row for row in range(len(ref_none)) if row not in eng_stood.values()]
        eng_rows = [row for row in range(len(eng_none)) if row not in ref_stood.values()]
        ref_rows += range(len(ref_none), len(ref_none) + len(eng_free))
        eng_rows += range(len(eng_none), len(eng_none) + len(ref_free))
        for at, (ref_at, eng_at) in enumerate(zip_longest(ref_free, eng_free)):
            ref_value, ref_nth = (None, ref_rows[at - len(ref_free)]) if ref_at is None else ref_items[ref_at]
            eng_value, eng_nth = (None, eng_rows[at - len(eng_free)]) if eng_at is None else eng_items[eng_at]
            found.append((ref_value, eng_value, ref_nth, eng_nth))
        return found

    def judge(
        self, key: bytes, name: bytes, ref_value: bytes | None, eng_value: bytes | None, ref_nth: int, eng_nth: int
    ) -> str | None:
        """The class that accepts the field's two unpaired values under ``key`` (a side that
        wrote none is ``None``): the first class of the field's source reading the field to
        accept, the applied ones consulted first. ``ref_nth`` and ``eng_nth`` count the other
        unpaired pairs of the key and field that carry the same value on that side, so that
        each pair reads rows of its own. Keys come in ascending byte order."""
        field = text(name)
        expected = formats.ABSENT if ref_value is None else text(ref_value)
        actual = formats.ABSENT if eng_value is None else text(eng_value)
        rows: tuple | None = None
        rows_read = False
        for klass in self._classes.reading(self._sources.of(name), field):
            if klass == formats.ORDER_ONLY:
                separator = self._order_only.get(field)
                paired = bool(separator) and ref_value is not None and eng_value is not None
                accepted = paired and formats.same_items(expected, actual, separator)
            elif klass == formats.PLUGIN_VALUE_ON_REFERENCE_SKIPPED_RECORD:
                accepted = ref_value is None and eng_value is not None and self._dropped_record(key, field)
            else:
                if not rows_read:
                    rows, rows_read = self._rows(key, field, expected, actual, ref_nth, eng_nth), True
                handed = (expected, actual)
                if rows is not None and klass in formats.TRANSCRIPT_READING_CLASSES:
                    handed = (self._as_written(rows[0], field, expected), self._as_written(rows[1], field, actual))
                keywords = self._handed.get(klass, {})
                if rows is not None and klass in formats.SIBLING_ALLELE_CLASSES:
                    keywords = {**keywords, "candidate_rows": self._candidate_rows(rows[3])}
                if rows is not None and klass in formats.RECORD_REFERENCE_CLASSES:
                    self._read_record(rows[0], rows[4])
                    self._read_record(rows[1], rows[3])
                accepted = rows is not None and formats.field_class_accepts(
                    klass, rows[0], field, *handed, rows[1], record_rows=rows[2], **keywords
                )
            if accepted:
                return klass
        return None


def residual_path(directory: Path, source: str, engine: str, projection: str) -> Path:
    return directory / f"{source}.{engine}.{projection}.tsv"


class Residuals:
    """The residual files of one engine: per source and projection, every tuple outside the
    intersection that no applied class set aside, one line of its side and its columns."""

    def __init__(self, directory: Path, engine: str, written: list[Path]) -> None:
        self._directory = directory
        self._engine = engine
        self._written = written
        self._files: dict[tuple[str, str], io.BufferedWriter] = {}

    def _path(self, source: str, projection: str) -> Path:
        """The path of one residual file, noted among the files the run wrote."""
        path = residual_path(self._directory, source, self._engine, projection)
        self._written.append(path)
        return path

    def write(self, projection: str, source: str, prefix: bytes, unexplained: list[tuple[bytes, bytes]]) -> None:
        """Writes the ``(value, side)`` pairs left unexplained under ``prefix``, the tuple's
        columns before its value, to the source's file."""
        file = self._files.get((source, projection))
        if file is None:
            file = self._files[source, projection] = self._path(source, projection).open("wb")
        file.write(b"".join(b"%b\t%b%b\n" % (SIDE_NAMES[side], prefix, value) for value, side in sorted(unexplained)))

    def close(self) -> None:
        for file in self._files.values():
            file.close()

    def complete(self, scored: list[tuple[str, str]]) -> None:
        """Leaves an empty file for each scored ``(source, projection)`` without a residual."""
        for source, projection in scored:
            if (source, projection) not in self._files:
                self._path(source, projection).write_bytes(b"")


class Projection:
    """One projection of one engine: its tally, and the book-keeping after the unpaired
    tuples of one key are judged."""

    def __init__(self, name: str, sources: Sources, classes: Classes, residuals: Residuals | None) -> None:
        self.tally = Tally()
        self._name = name
        self._sources = sources
        self._applied = classes.applied
        self._order = classes.order
        self._residuals = residuals
        self._shared_key: dict[str, bytes] = {}

    def shared(self, head: bytes) -> None:
        """Counts the key of ``head``, a tuple's columns before its value, once for the
        source of its field; heads come in ascending byte order."""
        cut = head.rindex(b"\t")
        source = self._sources.of(head[cut + 1 :])
        if self._shared_key.get(source) != head[:cut]:
            self._shared_key[source] = head[:cut]
            self.tally.shared[source] += 1

    def _record(self, source: str, prefix: bytes, accepted: set[str], unexplained: list[tuple[bytes, bytes]]) -> None:
        """Counts the key once under each class that accepted a tuple of it and writes the
        ``(value, side)`` pairs left unexplained."""
        for klass in accepted:
            (self.tally.masked if klass in self._applied else self.tally.not_applied)[source][klass] += 1
        if unexplained and self._residuals is not None:
            self._residuals.write(self._name, source, prefix, unexplained)


class FeatureScore(Projection):
    """The feature projection of one engine. Each tuple left unpaired is written to
    ``verdicts`` without its Feature, with its side and the class that accepted it (empty
    for none), the lines the allele projection reads."""

    def __init__(
        self, context: RowContext, verdicts: Sink, sources: Sources, classes: Classes, residuals: Residuals | None
    ) -> None:
        super().__init__("feature", sources, classes, residuals)
        self._context = context
        self._verdicts = verdicts

    def settle(self, prefix: bytes, ref_values: list[bytes], eng_values: list[bytes]) -> None:
        tally = self.tally
        body = prefix[:-1]
        cut = body.rindex(b"\t") + 1
        key, name = body[:cut], body[cut:]
        ref_left, eng_left = tally.pair(name, ref_values, eng_values)
        allele_key = key[: key.rindex(b"\t", 0, -1) + 1]
        source = self._sources.of(name)
        accepted: set[str] = set()
        unexplained: list[tuple[bytes, bytes]] = []
        for ref_value, eng_value, ref_nth, eng_nth in self._context.pairs(key, name, ref_left, eng_left):
            klass = self._context.judge(key, name, ref_value, eng_value, ref_nth, eng_nth)
            tag = b"" if klass is None else klass.encode()
            for value, side, aside in (
                (ref_value, REFERENCE_SIDE, tally.ref_aside),
                (eng_value, ENGINE_SIDE, tally.eng_aside),
            ):
                if value is None:
                    continue
                self._verdicts.write(b"%b%b\t%b\t%b\t%b\n" % (allele_key, name, value, side, tag))
                if klass in self._applied:
                    aside[source] += 1
                else:
                    unexplained.append((value, side))
            if klass is not None:
                accepted.add(klass)
        self._record(source, prefix, accepted, unexplained)


class AlleleScore(Projection):
    """The allele projection of one engine, the classes of its unpaired tuples read from the
    feature projection's sorted verdict lines."""

    def __init__(self, verdicts: Path, sources: Sources, classes: Classes, residuals: Residuals | None) -> None:
        super().__init__("allele", sources, classes, residuals)
        self._verdicts = {REFERENCE_SIDE: Cursor(verdicts), ENGINE_SIDE: Cursor(verdicts)}

    def close(self) -> None:
        for cursor in self._verdicts.values():
            cursor.close()

    def settle(self, prefix: bytes, ref_values: list[bytes], eng_values: list[bytes]) -> None:
        tally = self.tally
        name = prefix[:-1].rsplit(b"\t", 1)[1]
        ref_left, eng_left = tally.pair(name, ref_values, eng_values)
        source = self._sources.of(name)
        accepted: set[str] = set()
        unexplained: list[tuple[bytes, bytes]] = []
        for values, side, aside in (
            (ref_left, REFERENCE_SIDE, tally.ref_aside),
            (eng_left, ENGINE_SIDE, tally.eng_aside),
        ):
            for lookup, value in sorted((b"%b%b\t%b\t" % (prefix, value, side), value) for value in values):
                # The verdict without a class sorts first, so a class there means no feature
                # tuple of this allele tuple was left unexplained.
                tags = [text(line[len(lookup) :]) for line in self._verdicts[side].starting(lookup)]
                if tags and tags[0] and all(tag in self._applied for tag in tags):
                    aside[source] += 1
                    accepted.add(min(tags, key=self._order.__getitem__))
                else:
                    unexplained.append((value, side))
                    if tags and tags[0]:
                        held = [tag for tag in tags if tag not in self._applied]
                        accepted.add(min(held, key=self._order.__getitem__))
        self._record(source, prefix, accepted, unexplained)


def unshared(head: bytes) -> None:
    """Stands in for ``Projection.shared`` over a remainder, whose keys a side writes a field
    under more than once the pairing tool has counted over every tuple."""


def score_engine(
    reference: Extracted,
    engine: Extracted,
    sorter: Sorter,
    verdicts: Path,
    sources: Sources,
    classes: Classes,
    residuals: Residuals | None,
    plugin_reads: dict[str, dict[str, object]],
    run_input: InputIndex | None,
    counted: dict | None = None,
) -> dict[str, Tally]:
    """The tallies of one engine against the reference, per projection. With ``counted``,
    what the pairing tool counted of the two, their tuple files are its remainder: the
    merge hands every line of them to the settling step, and the tool's pairs and its
    shared keys are added to what that step tallies."""
    context = RowContext(reference, engine, sources, classes, plugin_reads, run_input)
    try:
        with sorter.sinks(verdicts) as (sink,):
            feature = FeatureScore(context, sink, sources, classes, residuals)
            shared = feature.shared if counted is None else unshared
            merge(reference.feature, engine.feature, feature.tally.paired, feature.settle, shared)
    finally:
        context.close()
    allele = AlleleScore(verdicts, sources, classes, residuals)
    try:
        shared = allele.shared if counted is None else unshared
        merge(reference.allele, engine.allele, allele.tally.paired, allele.settle, shared)
    finally:
        allele.close()
    tallies = {"feature": feature.tally, "allele": allele.tally}
    if counted is not None:
        for projection, tally in tallies.items():
            tally.take(counted[projection])
    return tallies


def source_entries(tally: Tally, sources: Sources, classes: Classes) -> dict[str, dict]:
    """The tally summed per source under the report's count key names. ``masked`` lists every
    applied class of the source, and an entry in which a class that is not applied would set
    a tuple aside carries no adjusted count."""
    sums: dict[str, list[int]] = {}
    for name in tally.fields():
        counts = sums.setdefault(sources.of(name), [0, 0, 0, 0])
        counts[0] += tally.paired.get(name, 0)
        counts[1] += tally.ref_left.get(name, 0)
        counts[2] += tally.eng_left.get(name, 0)
        counts[3] += tally.text_only.get(name, 0)
    entries: dict[str, dict] = {}
    for source, (paired, ref_left, eng_left, text_only) in sums.items():
        perl, rust = paired + ref_left, paired + eng_left
        adj_perl = perl - tally.ref_aside.get(source, 0)
        adj_rust = rust - tally.eng_aside.get(source, 0)
        not_applied = dict(sorted(tally.not_applied.get(source, {}).items()))
        entry = {"perl": perl, "rust": rust, "intersection": paired}
        if not not_applied:
            entry |= {"adj_perl": adj_perl, "adj_rust": adj_rust, "adj_intersection": paired}
        masked = tally.masked.get(source, {})
        entries[source] = entry | {
            "masked": {klass: masked.get(klass, 0) for klass in classes.of(source)},
            "not_applied": not_applied,
            "text_only_differences": text_only,
            "residual_rows": adj_perl + adj_rust - 2 * paired,
            "shared_keys": tally.shared.get(source, 0),
        }
    return entries


def score(
    reference: Path,
    engines: list[tuple[str, Path]],
    plugins: list[str],
    custom: dict[str, str],
    work_dir: Path | None,
    *,
    sort_memory: str = "2G",
    sort_threads: int | None = None,
    residuals_dir: Path | None = None,
    class_table: Path | None = None,
    applied: list[str] | None = None,
    populations: dict[str, int] | None = None,
    engine_facts: dict[str, dict[str, str]] | None = None,
    run_command: list[str] | None = None,
    gwas_catalog: Path | None = None,
    revel_table: Path | None = None,
    revel_no_match: bool = False,
    run_input: Path | None = None,
    utrannotator_sweep: list[Path] | None = None,
    json_cache: Path | None = None,
    pairing_tool: Path | None = None,
    pairing_threads: int | None = None,
    met_alone: bool = False,
    spliceai_table: list[Path] | None = None,
) -> dict:
    """The report for ``reference`` against each ``(id, path)`` of ``engines``, over the
    plugins the scored run loaded and its custom sources. With ``residuals_dir``, each
    engine's residual files are written there. ``class_table`` declares the classes,
    ``applied`` names those applied (all of them when ``None``, less with ``met_alone`` the rows
    whose fourth column reads ``not met``), ``populations`` gives a class
    the keys it can reach, ``engine_facts`` gives an engine id the facts its report entry
    carries, ``run_command`` holds the words of the scored run's command line,
    ``json_cache`` names the candidate JSON cache of the scored run,
    ``spliceai_table`` names the score files its SpliceAI plugin read,
    ``utrannotator_sweep`` names the outputs of Ensembl VEP's seeded runs of the input, or
    directories that hold them,
    ``revel_table`` is the table its REVEL plugin read and ``revel_no_match`` whether the run
    set that plugin's ``no_match``,
    ``gwas_catalog`` is the associations file its GWAS plugin read, and ``run_input`` is the
    VCF the run annotated. ``pairing_tool`` names the built pairing tool, which then extracts,
    sorts and pairs in place of ``extract`` and of the merge's pairs, on ``pairing_threads``
    threads or on every core."""
    if not engines:
        raise ScorerError("no engine to score")
    for ident in engine_facts or {}:
        if all(ident != other for other, _ in engines):
            raise ScorerError(f"{ident} is no engine of this run")
    sources = Sources(plugins, custom, run_command)
    classes = Classes(class_table, applied, populations or {}, sources, met_alone)
    check_first_line(reference, tuple_file_allowed=False)
    for _, path in engines:
        check_first_line(path, tuple_file_allowed=True)
    # What the classes of a plugin are handed beside the rows, by the registry's plugin name:
    # the keywords of ``field_class_accepts`` that plugin's classes read.
    plugin_reads: dict[str, dict[str, object]] = {}
    catalog = None
    if gwas_catalog is not None:
        catalog = digested(gwas_catalog)
        plugin_reads[GWAS_PLUGIN] = {"gwas_pairs": formats.GwasCatalogPairs.from_catalog(gwas_catalog)}
    table = None
    if revel_table is not None:
        table = {**digested(revel_table), "no_match": revel_no_match}
        plugin_reads[REVEL_PLUGIN] = {"revel_table": formats.RevelTable(revel_table, revel_no_match)}
    sweep = None
    if utrannotator_sweep:
        outputs = default_format_outputs(utrannotator_sweep)
        sweep = [digested(path) for path in outputs]
        plugin_reads[UTRANNOTATOR_PLUGIN] = {"utrannotator_sweep": formats.UtrannotatorSeedSweep.from_outputs(outputs)}
    cache = None
    if json_cache is not None:
        cache = {"path": str(json_cache)}
        plugin_reads[DOWNSTREAM_PLUGIN] = {
            "downstream_sequences": formats.DownstreamSequences(json_cache),
            "declared": classes.applied,
        }
        if classes.of_plugin(BLOSUM62_PLUGIN):
            plugin_reads[BLOSUM62_PLUGIN] = {"downstream_sequences": plugin_reads[DOWNSTREAM_PLUGIN]["downstream_sequences"]}
        if classes.of_plugin(MAXENTSCAN_PLUGIN):
            plugin_reads[MAXENTSCAN_PLUGIN] = {"splice_windows": formats.SpliceWindows(json_cache)}
    score_files = None
    score_table = None
    spliceai_refusals: list[dict] | None = None
    if spliceai_table:
        score_files = [{"path": str(path)} for path in spliceai_table]
        score_table = formats.SpliceAiTable(spliceai_table, json_cache)
        spliceai_refusals = []
        plugin_reads[SPLICEAI_PLUGIN] = {"spliceai_table": score_table}
    if residuals_dir is not None:
        for name in (*sources.names, *(ident for ident, _ in engines)):
            if IDENT_RE.fullmatch(name) is None:
                raise ScorerError(f"{name!r} cannot be part of a residual file name")
        names = [f"{source}.{ident}" for ident, _ in engines for source in sources.names]
        if len(set(names)) != len(names):
            raise ScorerError("two scored pairs of a source and an engine share a residual file name")
        if residuals_dir.is_dir() and any(residuals_dir.glob("*.tsv")):
            raise ScorerError(f"--residuals {residuals_dir} holds .tsv files of an earlier run")
    for directory in (work_dir, residuals_dir):
        if directory is not None:
            directory.mkdir(parents=True, exist_ok=True)
    if residuals_dir is not None and not os.access(residuals_dir, os.W_OK):
        raise ScorerError(f"cannot write in --residuals {residuals_dir}")
    seen: dict[str, set[str]] = defaultdict(set)
    scores: dict[tuple[str, int, str], dict] = {}
    scored: list[Extracted] = []
    residual_files: list[Path] = []
    # A row is read only by a class, so a run that declares none keeps no rows.
    keep_rows = bool(classes.declared)
    scratch: tempfile.TemporaryDirectory | None = None
    named_input: InputIndex | None = None
    read_input: dict | None = None
    try:
        scratch = tempfile.TemporaryDirectory(prefix=f"{COMPONENT}_", dir=work_dir)
        with scratch, ExitStack() as opened:
            work = Path(scratch.name)
            sorter = Sorter(work, sort_memory, sort_threads)
            if run_input is not None:
                named_input = InputIndex(run_input, sorter, work)
                opened.callback(named_input.close)
                read_input = {"path": str(run_input), "sha256": named_input.sha256, "records": named_input.records}
            pairing = None
            if pairing_tool is None:
                ref = extract(
                    "reference", reference, "input0", sources, sorter, work, tuple_file_allowed=False, keep_rows=keep_rows
                )
            else:
                pairing = Pairing(pairing_tool, reference, engines, sources, work, keep_rows, pairing_threads)
                ref = pairing.reference
            if ref.lines == 0:
                raise ScorerError(f"{reference}: the reference has no data row")
            for index, (ident, path) in enumerate(engines, 1):
                if pairing is None:
                    engine = extract(
                        ident, path, f"input{index}", sources, sorter, work, tuple_file_allowed=True, keep_rows=keep_rows
                    )
                    against, counted = ref, None
                else:
                    against, engine, counted = pairing.engine(index)
                residuals = None if residuals_dir is None else Residuals(residuals_dir, ident, residual_files)
                if score_table is not None:
                    score_table.refusals.clear()
                try:
                    verdicts = work / f"input{index}.verdicts"
                    tallies = score_engine(
                        against, engine, sorter, verdicts, sources, classes, residuals, plugin_reads, named_input, counted
                    )
                finally:
                    if residuals is not None:
                        residuals.close()
                if score_table is not None and spliceai_refusals is not None:
                    refused = score_table.refusals.counts()
                    if refused["asked"]:
                        spliceai_refusals.append({"engine": ident, **refused})
                for part in ("feature", "allele", "rows", "verdicts"):
                    (work / f"input{index}.{part}").unlink(missing_ok=True)
                scored.append(engine)
                for projection, tally in tallies.items():
                    for name in tally.fields():
                        seen[sources.of(name)].add(text(name))
                    for source, entry in source_entries(tally, sources, classes).items():
                        scores[source, index, projection] = {
                            "source": source,
                            "engine": ident,
                            "projection": projection,
                            **entry,
                        }
                if residuals is not None:
                    residuals.complete([(source, projection) for source, at, projection in scores if at == index])
        # A signal can still end the run here, so the report is put together where a failure
        # removes the residual files.
        populated = {
            source for (source, _, projection), entry in scores.items() if projection == "feature" and entry["perl"]
        }
        return {
            "tuple": PROJECTIONS,
            "value_rule": {
                "numeric_tolerance": str(TOLERANCE),
                "numeric_tolerance_applies": ["to the difference as a share of the larger number", "to the difference"],
            },
            "plugins": list(sources.plugins),
            "run_command": run_command,
            "class_table": classes.table,
            "gwas_catalog": catalog,
            "revel_table": table,
            "utrannotator_sweep": sweep,
            "json_cache": cache,
            "spliceai_table": score_files,
            "spliceai_refusals": spliceai_refusals,
            "run_input": read_input,
            "classes": [
                {
                    "class": name,
                    "plugins": {plugin: list(patterns) for plugin, patterns in plugins.items()},
                    "applied": name in classes.applied,
                    "further_columns": classes.further[name],
                }
                | ({"population": classes.populations[name]} if name in classes.populations else {})
                for name, plugins in classes.declared.items()
            ],
            "reference": input_entry(ref),
            "engines": [
                {"id": engine.ident, **(engine_facts or {}).get(engine.ident, {}), **input_entry(engine)}
                for engine in scored
            ],
            "sources": {source: sorted(seen[source]) for source in sorted(seen)},
            "scores": [scores[key] for key in sorted(scores, key=lambda key: (key[0], key[1], key[2] != "feature"))],
            "not_populated": sorted(set(sources.names) - populated),
        }
    except BaseException:
        # A stop signal can end the removal of the scratch on leaving the block above; taken
        # once, it cannot end this one.
        if scratch is not None:
            with suppress(OSError):
                scratch.cleanup()
        for path in residual_files:
            path.unlink(missing_ok=True)
        raise


def digested(path: Path) -> dict[str, str]:
    with path.open("rb") as handle:
        return {"path": str(path), "sha256": hashlib.file_digest(handle, "sha256").hexdigest()}


def names_default_columns(path: Path) -> bool:
    """Whether the file opens as a default-format output: its first line that is no ``##``
    comment is the column header."""
    try:
        reader = Reader(path)
    except OSError:
        return False
    try:
        head = reader.lines.read(1 << 16)
    except (OSError, EOFError, zlib.error, gzip.BadGzipFile):
        return False
    finally:
        reader.close()
    first = next((line for line in head.split(b"\n") if not line.startswith(b"##")), b"")
    return first.startswith(b"#Uploaded_variation\t")


def default_format_outputs(paths: list[Path]) -> list[Path]:
    """The outputs ``paths`` name: a file as it is, and of a directory every file in it that
    opens as a default-format output, in name order."""
    outputs: list[Path] = []
    for path in paths:
        if not path.is_dir():
            outputs.append(path)
            continue
        inside = sorted(entry for entry in path.iterdir() if entry.is_file() and names_default_columns(entry))
        if not inside:
            raise ScorerError(f"{path}: no default-format output in the directory")
        outputs += inside
    return outputs


def input_entry(extracted: Extracted) -> dict:
    entry = {"path": str(extracted.path), "sha256": extracted.sha256, "format": extracted.kind}
    if extracted.kind == "default":
        entry["rows"] = extracted.lines
    else:
        entry["tuples"] = extracted.lines
        entry["tuples_outside_every_source"] = extracted.outside
    entry["unclaimed_keys"] = extracted.unclaimed
    return entry


def replaced_whole(path: Path) -> bool:
    """Whether a report bound for ``path`` is written beside it and moved onto it, so that the
    path holds a whole report or what it held before: a path that is free or a regular file.
    A link, a device or a pipe is written as it is, since the move would replace it."""
    return not os.path.lexists(path) or (path.is_file() and not path.is_symlink())


def write_report(report: dict, path: Path, residuals_dir: Path | None) -> None:
    """Writes the report, whole or not at all where ``replaced_whole`` holds, a report that
    replaces another taking its mode; a run whose report is not written keeps none of its
    residual files."""
    staged = path
    try:
        if replaced_whole(path):
            # The process id keeps two runs bound for one path from writing one staged file.
            staged = path.with_name(f"{path.name}.{os.getpid()}.tmp")
        staged.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        if staged != path:
            # A free path has no mode to take, and a file system may refuse to set one.
            with suppress(OSError):
                shutil.copymode(path, staged)
            os.replace(staged, path)
    except BaseException:
        if staged != path:
            with suppress(OSError):
                staged.unlink()
        if residuals_dir is not None:
            for entry in report["scores"]:
                written = residual_path(residuals_dir, entry["source"], entry["engine"], entry["projection"])
                written.unlink(missing_ok=True)
        raise


@contextmanager
def one_interrupt() -> Iterator[None]:
    """Takes an interrupt or a termination signal as one ``KeyboardInterrupt`` meanwhile, so
    that either stops the sorts and empties the scratch, and ignores a repeat of either while
    the run stops. A signal the process ignores stays ignored. A thread that cannot set a
    signal handler runs as it is, and so does a signal whose handler was not set from Python:
    it reads as ``None``, which cannot be set again afterwards."""
    previous: dict[int, object] = {}

    def interrupt(number: int, frame: object) -> None:
        for taken in previous:
            signal.signal(taken, signal.SIG_IGN)
        raise KeyboardInterrupt

    for number in (signal.SIGINT, signal.SIGTERM):
        handler = signal.getsignal(number)
        if handler is None or handler is signal.SIG_IGN:
            continue
        try:
            signal.signal(number, interrupt)
        except ValueError:
            break
        previous[number] = handler
    try:
        yield
    finally:
        for number, handler in previous.items():
            signal.signal(number, handler)


def named_pairs(option: str, values: list[str], what: str) -> list[tuple[str, str]]:
    """``NAME=VALUE`` arguments as pairs, each name once."""
    pairs: list[tuple[str, str]] = []
    for value in values:
        name, separator, rest = value.partition("=")
        if not separator or not rest or IDENT_RE.fullmatch(name) is None:
            raise ScorerError(f"{option} takes {what}, not {value!r}")
        if any(name == other for other, _ in pairs):
            raise ScorerError(f"{option} names {name} twice")
        pairs.append((name, rest))
    return pairs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True, type=Path, help="Ensembl VEP default-format output")
    parser.add_argument(
        "--engine",
        required=True,
        action="append",
        metavar="ID=PATH",
        help="a scored engine's default-format output or tuple file; repeat per engine",
    )
    parser.add_argument(
        "--engine-binary",
        action="append",
        default=[],
        metavar="ID=SHA256",
        help="the sha256 of the binary that wrote an engine's output, carried into the report",
    )
    parser.add_argument(
        "--engine-commit",
        action="append",
        default=[],
        metavar="ID=COMMIT",
        help="the source commit of the binary that wrote an engine's output, carried into the report",
    )
    parser.add_argument(
        "--plugins",
        default="",
        metavar="A,B",
        help=(
            "the plugins the scored run loaded, by their registry names or aliases;"
            " a field is resolved among them alone"
        ),
    )
    parser.add_argument(
        "--run-command",
        default=None,
        metavar="COMMAND",
        help="the scored run's command line, whose flags the registry reads beside the run's plugins",
    )
    parser.add_argument(
        "--custom-source",
        action="append",
        default=[],
        metavar="NAME=PREFIX",
        help="group the fields starting with PREFIX as the source NAME; repeat per source",
    )
    parser.add_argument("--work-dir", type=Path, default=None, help="directory for the sort scratch")
    parser.add_argument("--report", required=True, type=Path, help="JSON report path")
    parser.add_argument(
        "--residuals",
        type=Path,
        default=None,
        metavar="DIR",
        help="directory for the residual files, <source>.<engine>.<projection>.tsv",
    )
    parser.add_argument(
        "--class-table",
        type=Path,
        default=DEFAULT_CLASS_TABLE,
        metavar="TSV",
        help=(
            "the documented classes: per line a class, its plugins and their field names"
            " (default: plugin_class_table.tsv beside this script)"
        ),
    )
    parser.add_argument(
        "--gwas-catalog",
        type=Path,
        default=None,
        metavar="FILE",
        help="the associations file the run's GWAS plugin read, handed to the classes of GWAS",
    )
    parser.add_argument(
        "--revel-table",
        type=Path,
        default=None,
        metavar="FILE",
        help="the table the run's REVEL plugin read, its tabix index beside it, handed to the classes of REVEL",
    )
    parser.add_argument(
        "--revel-no-match",
        action="store_true",
        help="the run set the REVEL plugin's no_match: a class of REVEL reads no row's transcript ids",
    )
    parser.add_argument(
        "--utrannotator-sweep",
        action="append",
        type=Path,
        default=[],
        metavar="PATH",
        help=(
            "an Ensembl VEP default-format output of the input under one Perl hash seed, or a directory"
            " of them, handed to the classes of UTRAnnotator; repeat per run"
        ),
    )
    parser.add_argument(
        "--json-cache",
        type=Path,
        default=None,
        metavar="DIR",
        help=(
            "the candidate JSON cache of the scored run: the classes of Downstream and of Blosum62 are handed its "
            "sequences, the class of MaxEntScan its introns, and spliceai_multi_entry_cell its transcripts' symbols, "
            "which name the gene of a row that writes no SYMBOL"
        ),
    )
    parser.add_argument(
        "--spliceai-table",
        action="append",
        type=Path,
        default=[],
        metavar="FILE",
        help=(
            "a score file the run's SpliceAI plugin read, its tabix index beside it, handed to the classes of"
            " SpliceAI with the cache of --json-cache; repeat for the snv= file and the indel= file"
        ),
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        metavar="FILE",
        help="the VCF the scored run annotated, read for the records its rows name by an identifier, the first `;`-joined element of a record's ID column",
    )
    parser.add_argument(
        "--classes",
        default=None,
        metavar="A,B",
        help=(
            "the declared classes to apply (default: every one whose fourth column does not read"
            " 'not met'; an empty value: none)"
        ),
    )
    parser.add_argument(
        "--class-population",
        action="append",
        default=[],
        metavar="CLASS=N",
        help="the keys a class can reach in this input; a count above it is marked",
    )
    parser.add_argument(
        "--pairing-tool",
        type=Path,
        default=None,
        metavar="BINARY",
        help=(
            "the built tools/vep-plugin-pairs: it extracts, sorts and pairs on every core in place of this"
            " script's own extraction, sort and merge, and the report and the residual files are the same"
        ),
    )
    parser.add_argument(
        "--pairing-threads", type=int, default=None, metavar="N", help="threads of the pairing tool (default: every core)"
    )
    parser.add_argument("--sort-memory", default="2G", help="buffer of each sort process, as `sort -S` takes it")
    parser.add_argument("--sort-threads", type=int, default=None, help="threads of each sort process")
    args = parser.parse_args(argv)
    try:
        engines = [(ident, Path(path)) for ident, path in named_pairs("--engine", args.engine, "ID=PATH")]
        custom = dict(named_pairs("--custom-source", args.custom_source, "NAME=PREFIX"))
        populations = {}
        for name, count in named_pairs("--class-population", args.class_population, "CLASS=N"):
            if not count.isdecimal():
                raise ScorerError(f"--class-population takes CLASS=N, not {name}={count}")
            populations[name] = int(count)
        applied = None if args.classes is None else [name.strip() for name in args.classes.split(",") if name.strip()]
        plugins = [name.strip() for name in args.plugins.split(",") if name.strip()]
        try:
            run_command = None if args.run_command is None else shlex.split(args.run_command)
        except ValueError as error:
            raise ScorerError(f"--run-command takes the run's command line, not {args.run_command!r}") from error
        facts: dict[str, dict[str, str]] = {}
        for option, what, key, values, shape in (
            ("--engine-binary", "ID=SHA256", "binary_sha256", args.engine_binary, SHA256_RE),
            ("--engine-commit", "ID=COMMIT", "commit", args.engine_commit, COMMIT_RE),
        ):
            for ident, value in named_pairs(option, values, what):
                if shape.fullmatch(value) is None:
                    raise ScorerError(f"{option} takes {what}, not {ident}={value}")
                facts.setdefault(ident, {})[key] = value
        inputs = [args.reference, *(path for _, path in engines)]
        inputs += [
            path
            for path in (args.class_table, args.gwas_catalog, args.revel_table, args.input, args.pairing_tool)
            if path is not None
        ]
        inputs += [path for path in args.utrannotator_sweep if not path.is_dir()]
        inputs += args.spliceai_table
        for path in inputs:
            if not path.is_file():
                raise ScorerError(f"not a file: {path}")
        if args.json_cache is not None and not args.json_cache.is_dir():
            raise ScorerError(f"not a directory: {args.json_cache}")
        args.report.parent.mkdir(parents=True, exist_ok=True)
        if (
            args.report.is_dir()
            or (args.report.exists() and not os.access(args.report, os.W_OK))
            or (replaced_whole(args.report) and not os.access(args.report.parent, os.W_OK))
        ):
            raise ScorerError(f"cannot write the report {args.report}")
        with one_interrupt():
            report = score(
                args.reference,
                engines,
                plugins,
                custom,
                args.work_dir,
                sort_memory=args.sort_memory,
                sort_threads=args.sort_threads,
                residuals_dir=args.residuals,
                class_table=args.class_table,
                gwas_catalog=args.gwas_catalog,
                revel_table=args.revel_table,
                revel_no_match=args.revel_no_match,
                utrannotator_sweep=args.utrannotator_sweep,
                json_cache=args.json_cache,
                spliceai_table=args.spliceai_table,
                run_input=args.input,
                applied=applied,
                populations=populations,
                engine_facts=facts,
                run_command=run_command,
                pairing_tool=args.pairing_tool,
                pairing_threads=args.pairing_threads,
                met_alone=True,
            )
            try:
                write_report(report, args.report, args.residuals)
            except OSError as error:
                raise ScorerError(f"cannot write the report {args.report}: {error.strerror or error}") from error
    except (ScorerError, formats.ComparatorError) as error:
        print(f"ERROR: [{COMPONENT}] {error}", file=sys.stderr)
        return 2
    except OSError as error:
        print(f"ERROR: [{COMPONENT}] {' '.join(str(error).split())}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print(f"ERROR: [{COMPONENT}] interrupted before the report was written", file=sys.stderr)
        return 130
    # A key of an input can carry bytes the stream's encoding lacks.
    reconfigure = getattr(sys.stdout, "reconfigure", None)
    if reconfigure is not None:
        reconfigure(errors="backslashreplace")
    further = {
        item["class"]: f" [{', '.join(item['further_columns'])}]" if item["further_columns"] else ""
        for item in report["classes"]
    }
    for entry in report["scores"]:
        row = f"{entry['source']}\t{entry['engine']}\t{entry['projection']}"
        counts = " ".join(f"{key}={value}" for key, value in entry.items() if isinstance(value, int))
        print(f"{row}\t{counts}")
        for klass, count in entry["masked"].items():
            population = populations.get(klass)
            note = " (none)" if count == 0 else ""
            if population is not None and count > population:
                note = f" (above its reachable population of {population})"
            print(f"{row}\tclass {klass} set aside {count}{note}{further[klass]}")
        for klass, count in entry["not_applied"].items():
            print(f"{row}\tclass {klass} is not applied and would set aside {count}: raw counts only{further[klass]}")
    for whose, described in (
        ("the reference", report["reference"]),
        *((f"the engine {item['id']}", item) for item in report["engines"]),
    ):
        if described["unclaimed_keys"]:
            keys = ", ".join(f"{name} ({times})" for name, times in described["unclaimed_keys"].items())
            print(f"keys of {whose} that no source claims: {keys}")
    if report["not_populated"]:
        print(f"not populated in the reference: {', '.join(report['not_populated'])}")
    for plugin, key, option, clause in (
        (GWAS_PLUGIN, "gwas_catalog", "--gwas-catalog", "that reads the catalog sets nothing aside"),
        (REVEL_PLUGIN, "revel_table", "--revel-table", "that reads the table sets nothing aside"),
        (UTRANNOTATOR_PLUGIN, "utrannotator_sweep", "--utrannotator-sweep", "sets aside no row it reads the seeded runs for"),
        (DOWNSTREAM_PLUGIN, "json_cache", "--json-cache", "that reads the candidate cache sets nothing aside"),
        (MAXENTSCAN_PLUGIN, "json_cache", "--json-cache", "that reads the candidate cache sets nothing aside"),
    ):
        applied = any(item["applied"] and plugin in item["plugins"] for item in report["classes"])
        if report[key] is None and applied and plugin in report["plugins"]:
            print(f"no {option}: a class of {plugin} {clause}")
    reads_score_files = any(item["applied"] and item["class"] == formats.SPLICEAI_MULTI_ENTRY_CELL for item in report["classes"])
    if report["spliceai_table"] is None and reads_score_files and SPLICEAI_PLUGIN in report["plugins"]:
        print(f"no --spliceai-table: {formats.SPLICEAI_MULTI_ENTRY_CELL} sets nothing aside")
    elif report["json_cache"] is None and reads_score_files and SPLICEAI_PLUGIN in report["plugins"]:
        print(f"no --json-cache: {formats.SPLICEAI_MULTI_ENTRY_CELL} sets aside no row that writes no SYMBOL")
    for refused in report["spliceai_refusals"] or ():
        for line in formats.spliceai_refusal_lines(refused, f", engine {refused['engine']}"):
            print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
