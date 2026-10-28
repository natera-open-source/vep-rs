#!/usr/bin/env python3
"""Per-column comparison of VEP output files in every output format.

Three modes, one JSON report with a ``verdict`` each, exit 0 on ``green`` and 1 on ``red``:

``--format {default,tab,vcf,json,parquet} --expected F --actual F``
    One Ensembl VEP output (expected) against one candidate output (actual) in one format.
    Every row of the default and tab formats, every CSQ entry and record line of the VCF
    format, and every consequence object and record object of the JSON format is one
    comparable entry; every field of every entry is compared verbatim. A Parquet actual is
    read back through ``scripts/adapters/parquet_to_vep_tab.py`` and compared with a tab
    expected as a multiset. ``--input`` (the input VCF) keys the VCF and JSON entries by
    input record (a line's CHROM, POS, ID, REF and ALT, which both engines write back
    unchanged; a JSON object's ``input`` names its record the same way, since Ensembl VEP
    writes a structural variant's empty or ``0`` QUAL, FILTER or INFO back as ``.``), so a
    record one side skipped keeps every later record aligned.

``--self-parity --default F --tab F --vcf F --json F [--parquet P]``
    One engine's formats against each other, the default format as the pivot: every record,
    every (allele, feature) entry and every column the two formats both declare must carry
    the same value once the format's own encoding is undone.

``--warnings --expected F --actual F [--expected-skipped F --actual-skipped F]``
    Two engines' ``_warnings.txt`` files, and optionally their ``--skipped_variants_file``
    outputs, as sets of (input line number, description) with the quoted input line and any
    path masked.

Masks, floors and field pairs (comparison mode):

* ``--manifest`` is a golden-corpus ``manifest.json`` whose ``divergences``,
  ``vep_rs_only_tuples`` and ``field_divergences`` name keys documented to differ, each
  under an ``expected_divergence`` class; ``--masks CLASS,...`` names the classes applied.
  A documented consequence set or field value of an applied class is checked against the
  actual instead of the expected (an actual value other than the documented one is a stale
  divergence, red); a key documented as candidate-only is set aside, and a record the
  reference skipped keeps its record-level line and object uncompared; a key of a class not
  named is compared like any other. The class ``regulatory_feature_rows`` is built in rather
  than read from a manifest: it drops the expected ``RegulatoryFeature`` and ``MotifFeature``
  rows in every format for a candidate cache without regulatory data, recomputes the
  ``most_severe_consequence`` of each JSON record that lost a row from the rows that remain,
  and in ``--warnings`` mode requires the candidate's notice about missing regulatory data
  (or the notice of releases 0.1.0 to 0.5.0) and leaves it out of the comparison. The class
  ``order_only`` is built in too: the set-valued columns
  ``DOMAINS``, ``FREQS`` and ``ZYG`` (and the ``Extra`` text carrying them) agree when their
  items agree as a multiset in any order, since Ensembl VEP prints the protein features in
  the order of the cache's array and the frequencies and the per-sample zygosities in hash
  order; every such agreement is counted per column as ``order_only`` in the report, never
  passed as equal. The class ``individual_allele_order`` is built in as well: Ensembl VEP
  lists the two alternate alleles of an individual's copy of a record in hash order
  (``Parser/VCF.pm:639``), so an entry carrying ``IND`` (``ind`` on a JSON record object)
  agrees when the copy's alternate rows come in any order within each feature (a reference
  row its number or its reference allele identifies keeping its place), when its
  ``allele_string`` and ``UPLOADED_ALLELE`` name
  the same reference allele and the same alternate alleles in any order, and when its
  ``ALLELE_NUM`` values differ by one permutation of the alternates' numbers applied to every
  row of the copy, the reference allele keeping 0; every such agreement is counted per column
  as ``individual_allele_order``, the reordered entries under ``counts``, and an entry
  without ``IND`` is compared like any other. A consistent renumbering whose allele strings
  do not move is accepted too, since ``UPLOADED_ALLELE`` prints the record's own string
  wherever an original allele string is kept.
* ``--floors COL=MIN,...``: the populated (present and not ``-``) count of column ``COL``
  must reach ``MIN`` on both sides. ``COL`` is the default-format name in every format; in
  the JSON format the count is taken under the column's JSON key (lower-cased and renamed
  as the JSON factory renames it, SIFT and PolyPhen under either of their two keys, a
  position under either of its start and end keys), and a column with no row count there
  has no floor there: a record-level column, which the JSON format carries once per record
  or not at all, ``Feature_type``, which names the consequence list and is no key of its
  entries, and ``Extra``, the text formats' container of the keyed fields; the report lists
  such a column under ``floors_not_applicable`` and the text formats of the same case check
  its floor.
* ``--field-pairs 'COL=EXPECTED=>ACTUAL;...'`` (repeatable): a value pair accepted as
  agreeing wherever it occurs in ``COL``.
* ``--nearest-sites <json_cache>``: ``NEAREST`` values (and the JSON ``nearest``) agree
  through the candidate cache's protein-coding transcription start sites when every expected
  value is among the actual values and every actual value starts at a site of an expected
  value. Ensembl VEP names one transcript per winning site, whichever its interval tree held
  first, and the candidate names every transcript at the site, so the ids differ where the
  sites agree.
* ``--keep-only regulatory_feature_rows`` in comparison mode inverts the built-in row mask:
  both inputs retain only ``RegulatoryFeature`` and ``MotifFeature`` consequence rows,
  plus the VCF/JSON record entries of records retaining such a row. A JSON record that
  loses non-regulatory rows has its ``most_severe_consequence`` recomputed from its retained
  rows. Other record fields compare verbatim. This cannot be combined with
  ``--masks regulatory_feature_rows`` or used in warnings or self-parity mode.

A comparison is green when no expected entry is missing, no actual entry is unexpected
(the two one-sided key counts), every compared field agrees after masks and field pairs,
the shared entries come in the same order (not for Parquet, which is sorted by variant, and
not under ``--unordered``, for a run whose reference rows Ensembl VEP emits in hash order:
the allele buckets of ``--pick_allele`` and ``--pick_allele_gene`` and the gene buckets of
``--per_gene``, ``OutputFactory.pm:593``, ``:605`` and ``:855``; a golden corpus marks such
a flag set ``unordered``), no documented divergence is stale and every floor is met.

The parsers follow the comparators of ``crates/vep-cli/tests/common/mod.rs``. The output
grammar they read is Ensembl VEP's at release/116.2: ``OutputFactory/VEP_output.pm:99-136``
(the thirteen fixed columns and the ``Extra`` column of ``KEY=value`` pairs joined by ``;``,
``;`` in a value as ``%3B``), ``OutputFactory/Tab.pm:125-129`` and ``:187-198`` (every field
a column, ``-`` for an undefined value), ``OutputFactory/VCF.pm:100-119`` and ``:379-420``
(the CSQ fields in the ``##INFO`` header's ``Format:`` list, lists joined with ``&``, ``,``
and ``|`` rewritten to ``&``, ``;`` to ``%3B``, whitespace to ``_``, an undefined value
empty), ``OutputFactory/JSON.pm:93-136`` and ``:229-389`` (keys lower-cased, ``-`` values
dropped, ``YES`` as 1, ``*_position`` split into ``*_start`` and ``*_end``, SIFT and
PolyPhen split into prediction and score, ``DOMAINS`` into ``db``/``name`` objects, the
renames of ``%RENAME_KEYS``, the consequence lists named by feature type, numbers
numberified), ``OutputFactory.pm:883-931`` with ``Utils.pm:141-166`` (``Uploaded_variation``
and ``Location``), ``BaseVEP.pm:802-863`` and ``BaseRunner.pm:364-380`` (the warnings and
skipped-variants files), and ``%OVERLAP_CONSEQUENCES`` of ensembl-variation
``Utils/Constants.pm`` at release/116 (the term ranks behind ``most_severe_consequence``).
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import subprocess
import sys
import zlib
from collections import Counter, defaultdict, deque
from collections.abc import Callable, Collection, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TypeVar

COMPONENT = "compare_vep_formats"
ADAPTER = Path(__file__).resolve().parent.parent / "adapters" / "parquet_to_vep_tab.py"
FORMATS = ("default", "tab", "vcf", "json", "parquet")

# The field value of an entry that lacks the field.
ABSENT = "<absent>"
# The class a vep_rs_only_tuples entry carries when every record it names is one the
# reference dropped before annotation.
REFERENCE_SKIPPED_RECORD = "reference_skipped_record"
# The built-in mask class dropping the reference's regulatory rows (``without_regulatory_rows``).
REGULATORY_FEATURE_ROWS = "regulatory_feature_rows"
REGULATORY_KINDS = {"regulatory_feature", "motif_feature"}
# The 0.1.0-to-0.5.0 and the cache-aware warnings-file phrases for --regulatory on a cache
# without regulatory data.
REGULATORY_NOTICE = "regulatory-feature annotation is not implemented"
MISSING_REGULATORY_NOTICE = "this cache holds no regulatory data"
EXAMPLES_PER_GROUP = 4

# OutputFactory/JSON.pm:98-112 at release/116.2.
JSON_RENAME = {
    "consequence": "consequence_terms",
    "gene": "gene_id",
    "allele": "variant_allele",
    "symbol": "gene_symbol",
    "symbol_source": "gene_symbol_source",
    "overlapbp": "bp_overlap",
    "overlappc": "percentage_overlap",
    "refseq": "refseq_transcript_ids",
    "ensp": "protein_id",
    "chr": "seq_region_name",
    "variation_name": "id",
    "sv": "colocated_structural_variants",
    "clinical_impact": "clinvar_somatic_classification",
}
JSON_CONSEQUENCE_LISTS = (
    ("transcript_consequences", "transcript"),
    ("intergenic_consequences", "intergenic"),
    ("regulatory_feature_consequences", "regulatory_feature"),
    ("motif_feature_consequences", "motif_feature"),
)
# The engine identity vep-rs writes in every JSON object and Ensembl VEP never writes, set
# aside like ``input`` rather than compared; its header forms (``## vep-rs <version>``,
# ``##vep-rs=<version>``) fall with the other ``##`` lines the text parsers skip.
JSON_IDENTITY_KEY = "vep_rs_version"
# Columns the default format prints on every row of a record while the JSON format
# carries them once per record, in colocated_variants or at the top level, or not at all:
# the JSON factory's add_colocated_variant_info returns its hash unchanged (JSON.pm:493-495),
# so FREQS (OutputFactory.pm:1149-1163) has no JSON key, JSON.pm:419-426 copies the AF
# keys alone out of the per-allele hash that carries MAX_AF and MAX_AF_POPS, and ZYG
# (OutputFactory.pm:974-988) is a top-level key of the record object.
RECORD_LEVEL_COLUMNS = {
    "Uploaded_variation",
    "Location",
    "Existing_variation",
    "VARIANT_CLASS",
    "MINIMISED",
    "NEAREST",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "VAR_SYNONYMS",
    "FREQS",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLINVAR_SOMATIC_CLASSIFICATION",
    "ZYG",
}
RECORD_LEVEL_COLUMN_RE = re.compile(r"(^|_)AF$")
# The built-in class under which a set-valued column agrees when its items agree as a
# multiset: Ensembl VEP prints DOMAINS in the order of the cache's protein_features array
# (ensembl-variation BaseTranscriptVariation.pm:634-643 at release/116 copies it from the
# core API's Translation.pm:772-819, which buckets a translation's features by analysis and
# flattens the buckets in Perl hash order, the database's rows inside each), FREQS in hash
# order (OutputFactory.pm:1152-1153) and the ZYG list of --individual_zyg in the hash order
# of its samples (OutputFactory.pm:983-987), so none of the orders carries meaning.
ORDER_ONLY = "order_only"
ORDER_ONLY_COLUMNS = ("DOMAINS", "FREQS", "ZYG")
# The marker ``order_only_fields`` gives a JSON array field in place of a list separator.
JSON_ARRAY = "json"
# The JSON field of each order-only column with its separator: the ``domains`` array of a
# consequence entry, and the record entry's ``zyg`` array, whose strings ``json_scalar``
# joins with ``,``; FREQS has no JSON key.
ORDER_ONLY_JSON_FIELDS = {"domains": JSON_ARRAY, "zyg": ","}
# The built-in class under which an individual's alleles agree as a set. Under --individual,
# Parser/VCF.pm:635-639 at release/116.2 collects the alternate alleles an individual's
# genotype names in a hash and writes the copy's allele string as the reference allele
# followed by the hash's keys, so a copy carrying two alternate alleles lists them in hash
# order, and that order runs through the copy's rows within each feature, its ALLELE_NUM
# values (the allele's index in the string, BaseVariationFeatureOverlapAllele.pm:217-218 at
# release/116) and the JSON record's allele_string (JSON.pm:241). UPLOADED_ALLELE follows it
# only on a record with no original allele string (OutputFactory.pm:1360 prefers that
# string, which minimise_alleles sets on an indel record and under --minimal, Parser.pm:852,
# 881, 978-982).
INDIVIDUAL_ALLELE_ORDER = "individual_allele_order"
# The record-level token of an individual's copy: IND on a default, tab or CSQ entry
# (OutputFactory.pm:971), ind on a JSON record object.
INDIVIDUAL_FIELDS = ("IND", "ind")
# The fields carrying the copy's allele string, and the allele's index in it.
ALLELE_STRING_FIELDS = ("allele_string", "UPLOADED_ALLELE", "uploaded_allele")
ALLELE_NUMBER_FIELDS = ("ALLELE_NUM", "allele_num")
# The fields carrying the reference allele of a --show_ref_allele row.
REFERENCE_ALLELE_FIELDS = ("REF_ALLELE", "ref_allele")
# The note of a floor the JSON format cannot count.
FLOOR_NOT_APPLICABLE = (
    "a column the JSON format has no row count for (a record-level column is carried once "
    "per record or not at all, Feature_type names the consequence list, Extra is the text "
    "formats' container), so the floor is checked in the text formats of the same case"
)
# Config.pm:501-508 at release/116.2: --vcf implies --symbol, --biotype and --numbers, so a
# VCF run carries these columns where a default run without those flags does not.
VCF_IMPLIED_COLUMNS = {"SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE", "EXON", "INTRON"}

# ensembl-variation Utils/Constants.pm %OVERLAP_CONSEQUENCES at release/116, SO term to rank.
CONSEQUENCE_RANK = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "feature_elongation": 9,
    "feature_truncation": 10,
    "inframe_insertion": 11,
    "inframe_deletion": 12,
    "missense_variant": 13,
    "protein_altering_variant": 14,
    "splice_donor_5th_base_variant": 15,
    "splice_region_variant": 16,
    "splice_donor_region_variant": 17,
    "splice_polypyrimidine_tract_variant": 18,
    "incomplete_terminal_codon_variant": 19,
    "start_retained_variant": 20,
    "stop_retained_variant": 21,
    "synonymous_variant": 22,
    "coding_sequence_variant": 23,
    "mature_miRNA_variant": 24,
    "5_prime_UTR_variant": 25,
    "3_prime_UTR_variant": 26,
    "non_coding_transcript_exon_variant": 27,
    "intron_variant": 28,
    "NMD_transcript_variant": 29,
    "non_coding_transcript_variant": 30,
    "coding_transcript_variant": 31,
    "upstream_gene_variant": 32,
    "downstream_gene_variant": 33,
    "TFBS_ablation": 34,
    "TFBS_amplification": 35,
    "TF_binding_site_variant": 36,
    "regulatory_region_ablation": 37,
    "regulatory_region_amplification": 38,
    "regulatory_region_variant": 39,
    "intergenic_variant": 40,
    "sequence_variant": 41,
}


class ComparatorError(Exception):
    """A file this comparator cannot read or an option it cannot apply."""


def read_text(path: Path) -> str:
    """The text of ``path``, gunzipped when its name ends in ``.gz``."""
    try:
        if str(path).endswith(".gz"):
            with gzip.open(path, "rt", encoding="utf-8") as fh:
                return fh.read()
        return Path(path).read_text(encoding="utf-8")
    except (OSError, EOFError, UnicodeDecodeError, zlib.error) as exc:
        raise ComparatorError(f"cannot read {path}: {getattr(exc, 'strerror', None) or exc}") from exc


def text_lines(text: str) -> list[str]:
    """The lines of ``text`` without their terminators, no trailing empty line."""
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return [line.removesuffix("\r") for line in lines]


def input_records(path: Path) -> list[str]:
    """The data lines of an input VCF, in order."""
    return [line for line in text_lines(read_text(path)) if not line.startswith("#")]


@dataclass
class Entry:
    """One row of the default or tab format, one CSQ entry or record line of the VCF
    format, or one consequence object or record object of the JSON format.

    ``key`` is the entry's identity within its format: the Location (default, tab) or the
    input record ordinal (VCF, JSON), the allele, the feature and the feature type, plus a
    ``#n`` counter on exact duplicates.
    """

    key: str
    record: int | None
    location: str | None
    allele: str
    feature: str
    consequence: str
    fields: dict[str, str] = field(default_factory=dict)


def dedup_keys(entries: list[Entry]) -> None:
    """Suffixes the second and later occurrences of a key with ``#n``."""
    seen: Counter[str] = Counter()
    for e in entries:
        seen[e.key] += 1
        n = seen[e.key]
        if n > 1:
            e.key = f"{e.key}#{n}"


def canonical_json_text(text: str) -> str:
    """A JSON object's text with its members sorted (``json_compact``), or ``text`` as it is
    when it is not a JSON object. Ensembl VEP encodes the ``GA4GH_VRS`` object with its keys
    in hash order, which differs from row to row, so the text formats compare the parsed
    object, as the JSON format's nested object does."""
    try:
        value = json.loads(text)
    except ValueError:
        return text
    return json_compact(value) if isinstance(value, dict) else text


def numberified_json_text(text: str) -> str:
    """``canonical_json_text`` with every string member that is a whole number made a number:
    the ``GA4GH_VRS`` object of the text formats carries its ``start`` as a string, and the
    JSON format's nested object has both bounds as numbers (``OutputFactory/JSON.pm:288``,
    ``numberify``)."""
    try:
        value = json.loads(text)
    except ValueError:
        return text
    if not isinstance(value, dict):
        return text

    def numberify(v):
        if isinstance(v, dict):
            return {k: numberify(x) for k, x in v.items()}
        if isinstance(v, list):
            return [numberify(x) for x in v]
        return int(v) if isinstance(v, str) and v.isdecimal() else v

    return json_compact(numberify(value))


# The fields whose text is a JSON object, compared in canonical form.
JSON_OBJECT_FIELDS = ("GA4GH_VRS",)


def parse_default(text: str) -> list[Entry]:
    """The default format: the columns of the ``#`` header line, the ``Extra`` column's
    ``key=value`` pairs becoming fields of their own while ``Extra`` keeps its raw text, a
    ``JSON_OBJECT_FIELDS`` value in canonical form in both."""
    cols: list[str] = []
    entries: list[Entry] = []
    for line in text_lines(text):
        if line.startswith("#"):
            if not line.startswith("##"):
                cols = line[1:].split("\t")
            continue
        vals = line.split("\t")
        if len(vals) != len(cols):
            raise ComparatorError(f"column count on line: {line}")
        fields = dict(zip(cols, vals))
        for name in JSON_OBJECT_FIELDS:
            if name in fields:
                fields[name] = canonical_json_text(fields[name])
        extra = fields.get("Extra")
        if extra is not None and extra != "-":
            pairs: list[str] = []
            for kv in extra.split(";"):
                k, sep, v = kv.partition("=")
                if sep and k in JSON_OBJECT_FIELDS:
                    v = canonical_json_text(v)
                    kv = f"{k}={v}"
                fields[k] = v
                pairs.append(kv)
            fields["Extra"] = ";".join(pairs)
        location = fields.get("Location", "")
        allele = fields.get("Allele", "")
        feature = fields.get("Feature", "")
        entries.append(
            Entry(
                key="|".join([location, allele, feature, fields.get("Feature_type", "")]),
                record=None,
                location=location,
                allele=allele,
                feature=feature,
                consequence=fields.get("Consequence", ""),
                fields=fields,
            )
        )
    dedup_keys(entries)
    return entries


def parse_tab(text: str) -> list[Entry]:
    """The tab format: every column named in the ``#`` header line."""
    return parse_default(text)


class InputIndex:
    """Input records by identity, for mapping an output record back to its ordinal.
    Identical input lines are matched in order of appearance."""

    def __init__(self, inputs: list[str], key: Callable[[str], str]):
        self.unmatched: dict[str, deque[int]] = defaultdict(deque)
        for i, line in enumerate(inputs):
            self.unmatched[key(line)].append(i)

    def record(self, key: str) -> int | None:
        queue = self.unmatched.get(key)
        return queue.popleft() if queue else None


def vcf_line_key(line: str) -> str:
    """A VCF data line's identity: CHROM, POS, ID, REF and ALT."""
    return "\t".join(line.split("\t")[:5])


def csq_format_fields(line: str) -> list[str]:
    """The CSQ field names of a ``##INFO=<ID=CSQ`` header line's ``Format:`` list."""
    parts = line.split("Format: ", 1)
    if len(parts) < 2:
        raise ComparatorError(f"CSQ header without a Format list: {line}")
    fmt = parts[1]
    while fmt.endswith('">'):
        fmt = fmt[:-2]
    return fmt.split("|")


def header_columns(text: str) -> list[str]:
    """The columns of a default or tab output's ``#`` header line, empty without one."""
    for line in text_lines(text):
        if line.startswith("#") and not line.startswith("##"):
            return line[1:].split("\t")
    return []


def parse_vcf(text: str, inputs: list[str]) -> list[Entry]:
    """The VCF format: one entry per CSQ entry, keyed by the input record ordinal, the
    allele, the feature and the feature type, plus one record entry per line carrying the
    seven fixed columns, the INFO without CSQ and whether a CSQ was present. A line
    matching no input record is keyed by its output ordinal."""
    csq_fields: list[str] = []
    entries: list[Entry] = []
    index = InputIndex(inputs, vcf_line_key)
    lines_seen = 0
    for line in text_lines(text):
        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = csq_format_fields(line)
            continue
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 8:
            raise ComparatorError(f"VCF line has fewer than 8 columns: {line}")
        record = index.record(vcf_line_key(line))
        rec_label = f"rec{record}" if record is not None else f"line{lines_seen}"
        lines_seen += 1
        info_rest: list[str] = []
        csq: str | None = None
        for kv in cols[7].split(";"):
            if kv.startswith("CSQ="):
                csq = kv[4:]
            else:
                info_rest.append(kv)
        fields = dict(zip(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"], cols[:7]))
        fields["INFO_without_CSQ"] = ";".join(info_rest)
        fields["CSQ_present"] = "yes" if csq is not None else "no"
        entries.append(
            Entry(
                key=f"{rec_label}|record",
                record=record,
                location=None,
                allele="",
                feature="",
                consequence="",
                fields=fields,
            )
        )
        if csq is None:
            continue
        for chunk in csq.split(","):
            vals = chunk.split("|")
            if len(vals) != len(csq_fields):
                raise ComparatorError(f"CSQ chunk width on {rec_label}: {chunk}")
            cfields = dict(zip(csq_fields, vals))
            allele = cfields.get("Allele", "")
            feature = cfields.get("Feature", "")
            entries.append(
                Entry(
                    key=f"{rec_label}|{allele}|{feature}|{cfields.get('Feature_type', '')}",
                    record=record,
                    location=None,
                    allele=allele,
                    feature=feature,
                    consequence=cfields.get("Consequence", ""),
                    fields=cfields,
                )
            )
    dedup_keys(entries)
    return entries


def _json_float(f: float) -> str:
    """A float as a JSON text, shortest round trip, the exponent without sign padding."""
    if math.isnan(f) or math.isinf(f):
        return "null"
    s = repr(f)
    if "e" in s:
        mantissa, exponent = s.split("e")
        return f"{mantissa}e{int(exponent)}"
    return s


def json_compact(v) -> str:
    """A JSON value as compact text with object keys sorted."""
    if v is None or isinstance(v, bool):
        return json.dumps(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return _json_float(v)
    if isinstance(v, str):
        return json.dumps(v, ensure_ascii=False)
    if isinstance(v, list):
        return "[" + ",".join(json_compact(x) for x in v) + "]"
    members = ",".join(f"{json.dumps(k, ensure_ascii=False)}:{json_compact(x)}" for k, x in sorted(v.items()))
    return "{" + members + "}"


def json_scalar(v) -> str:
    """A JSON value as one comparable string: a string itself, an array its items joined
    with ``,``, anything else its compact JSON text."""
    if isinstance(v, str):
        return v
    if isinstance(v, list):
        return ",".join(json_scalar(x) for x in v)
    return json_compact(v)


def json_objects(text: str) -> list[tuple[int, dict]]:
    """The (line ordinal, object) of every non-empty line of a JSON Lines output."""
    out = []
    for n, line in enumerate(text_lines(text)):
        if not line.strip():
            continue
        try:
            obj = json.loads(line)
        except ValueError as exc:
            raise ComparatorError(f"JSON line {n + 1} does not parse: {exc}") from exc
        if not isinstance(obj, dict):
            raise ComparatorError(f"JSON line {n + 1} is not an object")
        out.append((n, obj))
    return out


def json_record_label(index: InputIndex, n: int, obj: dict) -> tuple[int | None, str]:
    """The input record ordinal of a JSON object and its label. An object names its record
    by the first five columns of its ``input``, as a VCF line does: Ensembl VEP stores a
    structural variant's line with an empty or ``0`` QUAL, FILTER or INFO written as ``.``
    (``Parser/VCF.pm:464-466``), so the whole ``input`` may differ from the line as read."""
    inp = obj.get("input")
    record = index.record(vcf_line_key(inp) if isinstance(inp, str) else "")
    return record, (f"rec{record}" if record is not None else f"obj{n}")


def parse_json(text: str, inputs: list[str]) -> list[Entry]:
    """The JSON format: one record entry per object (the top-level scalars other than the
    engine identity, keyed by the input record the ``input`` line names) and one entry per
    element of each consequence list, keyed by the record, ``variant_allele``, the feature id
    and the list's feature type. An object whose ``input`` matches no record is keyed by its
    output ordinal."""
    index = InputIndex(inputs, vcf_line_key)
    entries: list[Entry] = []
    set_aside = {name for name, _ in JSON_CONSEQUENCE_LISTS} | {"input", JSON_IDENTITY_KEY}
    for n, obj in json_objects(text):
        record, rec_label = json_record_label(index, n, obj)
        fields = {k: json_scalar(v) for k, v in obj.items() if k not in set_aside}
        entries.append(
            Entry(
                key=f"{rec_label}|record",
                record=record,
                location=None,
                allele="",
                feature="",
                consequence="",
                fields=fields,
            )
        )
        for list_name, kind in JSON_CONSEQUENCE_LISTS:
            items = obj.get(list_name)
            if not isinstance(items, list):
                continue
            for item in items:
                if not isinstance(item, dict):
                    raise ComparatorError(f"{rec_label}: a {list_name} element is not an object")
                cfields = {k: json_scalar(v) for k, v in item.items()}
                allele = cfields.get("variant_allele", "")
                if kind == "intergenic":
                    feature = ""
                else:
                    feature = cfields.get(f"{kind}_id", "")
                entries.append(
                    Entry(
                        key=f"{rec_label}|{allele}|{feature}|{kind}",
                        record=record,
                        location=None,
                        allele=allele,
                        feature=feature,
                        consequence=cfields.get("consequence_terms", ""),
                        fields=cfields,
                    )
                )
    dedup_keys(entries)
    return entries


def parquet_as_tab(path: Path) -> str:
    """A Parquet output read back as tab text through the adapter."""
    proc = subprocess.run(
        [sys.executable, str(ADAPTER), str(path)], capture_output=True, text=True, check=False
    )
    if proc.returncode != 0:
        raise ComparatorError(f"the Parquet adapter failed on {path}: {' '.join(proc.stderr.split())}")
    return proc.stdout


def parse_format(fmt: str, path: Path, inputs: list[str] | None) -> list[Entry]:
    """Parses one output file in ``fmt``; VCF and JSON need the input records."""
    if fmt == "parquet":
        return parse_tab(parquet_as_tab(path))
    text = read_text(path)
    if fmt == "default":
        return parse_default(text)
    if fmt == "tab":
        return parse_tab(text)
    if inputs is None:
        raise ComparatorError(f"--format {fmt} needs --input, the input VCF")
    if fmt == "vcf":
        return parse_vcf(text, inputs)
    return parse_json(text, inputs)


def dash(s: str) -> str:
    """The default format's spelling of absence for an empty VCF or JSON value."""
    return "-" if s == "" else s


def decode_extra(s: str) -> str:
    """Undoes the percent-encoding of ``=`` and ``;`` in Extra and CSQ values."""
    return s.replace("%3D", "=").replace("%3B", ";")


def documented_form(s: str) -> str:
    """A field value on the default format's spelling: ``-`` for an empty or absent
    value, percent-encoding undone."""
    return "-" if s == ABSENT else dash(decode_extra(s))


def consequence_set(s: str) -> str:
    """Consequence terms as a sorted, de-duplicated, comma-joined set."""
    terms = sorted({t for t in re.split(r"[&,]", s) if t})
    return ",".join(terms)


def _record_indices(item: dict) -> list[int]:
    out = []
    for r in item.get("record_indices") or []:
        if isinstance(r, int) and not isinstance(r, bool) and r >= 0:
            out.append(r)
    return out


def _s(item: dict, key: str) -> str:
    v = item.get(key)
    return v if isinstance(v, str) else ""


T = TypeVar("T")


class Documented:
    """The keys a corpus manifest documents as divergent, restricted to the classes applied:
    consequence sets per key, keys only the candidate emits, field values per key, and the
    records the reference skipped. Keys are addressable by Location or by record ordinal.
    ``masks`` is the set of classes requested, the built-in ones included."""

    def __init__(self) -> None:
        self.masks: set[str] = set()
        self.by_location: dict[tuple[str, str, str], list[str]] = defaultdict(list)
        self.by_record: dict[tuple[int, str, str], list[str]] = defaultdict(list)
        self.extra_by_location: set[tuple[str, str, str]] = set()
        self.extra_by_record: set[tuple[int, str, str]] = set()
        self.fields_by_location: dict[tuple[str, str, str], dict[str, str]] = defaultdict(dict)
        self.fields_by_record: dict[tuple[int, str, str], dict[str, str]] = defaultdict(dict)
        self.skipped_records: set[int] = set()
        self.available: Counter[str] = Counter()
        self.applied: Counter[str] = Counter()

    @classmethod
    def from_manifest(cls, manifest: dict | None, masks: set[str]) -> Documented:
        d = cls()
        d.masks = set(masks)
        if not manifest:
            return d
        for _, item in d._applied_items(manifest, "divergences", masks):
            sets = [x for x in item.get("vep_rs_consequence_sets") or [] if isinstance(x, str)]
            loc, allele, feature = _s(item, "location"), _s(item, "allele"), _s(item, "feature")
            d.by_location[(loc, allele, feature)].extend(sets)
            for r in _record_indices(item):
                d.by_record[(r, allele, feature)].extend(sets)
        for klass, item in d._applied_items(manifest, "vep_rs_only_tuples", masks):
            loc, allele, feature = _s(item, "location"), _s(item, "allele"), _s(item, "feature")
            d.extra_by_location.add((loc, allele, feature))
            for r in _record_indices(item):
                d.extra_by_record.add((r, allele, feature))
                if klass == REFERENCE_SKIPPED_RECORD:
                    d.skipped_records.add(r)
        for _, item in d._applied_items(manifest, "field_divergences", masks):
            fname, value = _s(item, "field"), _s(item, "vep_rs")
            loc, allele, feature = _s(item, "location"), _s(item, "allele"), _s(item, "feature")
            d.fields_by_location[(loc, allele, feature)][fname] = value
            for r in _record_indices(item):
                d.fields_by_record[(r, allele, feature)][fname] = value
        return d

    def _applied_items(self, manifest: dict, section: str, masks: set[str]) -> Iterator[tuple[str, dict]]:
        """The (class, item) pairs of one manifest section whose class ``masks`` names;
        every item's class is counted in ``available``, each yielded one's in ``applied``."""
        for item in manifest.get(section) or []:
            klass = _s(item, "expected_divergence")
            self.available[klass] += 1
            if klass in masks:
                self.applied[klass] += 1
                yield klass, item

    @staticmethod
    def _lookup(
        by_location: dict[tuple[str, str, str], T], by_record: dict[tuple[int, str, str], T], e: Entry
    ) -> T | None:
        """The value documented on the entry's key: by Location when the entry has one, then
        by record ordinal when it has one."""
        allele, feature = dash(e.allele), dash(e.feature)
        if e.location is not None:
            value = by_location.get((e.location, allele, feature))
            if value is not None:
                return value
        if e.record is not None:
            return by_record.get((e.record, allele, feature))
        return None

    def _divergent_fields(self, e: Entry) -> dict[str, str] | None:
        return self._lookup(self.fields_by_location, self.fields_by_record, e)

    def divergent_field(self, e: Entry, fname: str) -> str | None:
        """The documented candidate value of ``fname`` on the entry's key; the manifest
        spells the default-format name, the JSON format its lower-case form."""
        m = self._divergent_fields(e)
        if m is None:
            return None
        if fname in m:
            return m[fname]
        for k, v in m.items():
            if k.lower() == fname.lower():
                return v
        return None

    def has_divergent_fields(self, e: Entry) -> bool:
        return self._divergent_fields(e) is not None

    def divergent(self, e: Entry) -> list[str] | None:
        """The documented candidate consequence sets of the entry's key, if documented."""
        return self._lookup(self.by_location, self.by_record, e)

    def extra(self, e: Entry) -> bool:
        """Whether the entry's key is documented as one only the candidate emits."""
        allele, feature = dash(e.allele), dash(e.feature)
        if e.location is not None and (e.location, allele, feature) in self.extra_by_location:
            return True
        return e.record is not None and (e.record, allele, feature) in self.extra_by_record

    def record_divergent(self, record: int) -> bool:
        """Whether a documented key belongs to the record, so its summary values may differ."""
        return any(r == record for r, _, _ in self.by_record) or any(
            r == record for r, _, _ in self.extra_by_record
        )

    def reference_skipped(self, record: int) -> bool:
        return record in self.skipped_records


def populated_counts(entries: list[Entry]) -> Counter:
    """Per field, the entries carrying a value other than ``-``."""
    counts: Counter[str] = Counter()
    for e in entries:
        for name, value in e.fields.items():
            if dash(value) != "-":
                counts[name] += 1
    return counts


def extra_pairs(extra: str) -> dict[str, str]:
    """The ``key=value`` pairs of a default-format ``Extra`` value."""
    if extra in ("-", ABSENT):
        return {}
    return dict(kv.partition("=")[::2] for kv in extra.split(";"))


NEAREST_FIELDS = ("NEAREST", "nearest")


def nearest_list(value: str) -> list[str]:
    """The values of a ``NEAREST`` field as every format joins them (``,`` or ``&``)."""
    if value in ("-", ABSENT):
        return []
    return [v for v in re.split(r"[,&]", value) if v]


class NearestSites:
    """The protein-coding transcription start sites of a candidate JSON cache, by transcript
    stable id, gene stable id and gene symbol: a transcript's site is its start on the forward
    strand and its end on the reverse (``AnnotationType/Transcript.pm:749``), a gene's or
    symbol's the sites of its transcripts."""

    def __init__(self, sites: dict[str, set[int]]) -> None:
        self.sites = sites

    @classmethod
    def from_cache(cls, cache: Path) -> NearestSites:
        transcripts = cache / "transcripts"
        if not transcripts.is_dir():
            raise ComparatorError(f"--nearest-sites: {cache} has no transcripts/ directory")
        sites: dict[str, set[int]] = defaultdict(set)
        for shard in sorted(transcripts.glob("*/*.json")) + sorted(transcripts.glob("*/*.json.gz")):
            try:
                parsed = json.loads(read_text(shard))
            except ValueError as exc:
                raise ComparatorError(f"{shard} is not JSON: {exc}") from exc
            for t in cls.shard_transcripts(shard, parsed):
                if t.get("biotype") != "protein_coding":
                    continue
                forward = str(t.get("strand", "1")) == "1"
                site = cls.coordinate(shard, t, "start" if forward else "end")
                for key in ("stable_id", "gene_stable_id", "gene_symbol"):
                    value = t.get(key)
                    if value:
                        sites[value].add(site)
        return cls(dict(sites))

    @staticmethod
    def shard_transcripts(shard: Path, parsed: object) -> list[dict]:
        """The transcript objects of a parsed shard: a bare array, or a converted cache's
        range-keyed map of arrays; a ``null`` in either is a placeholder and is dropped."""
        if isinstance(parsed, dict):
            if not all(isinstance(wrapped, list) for wrapped in parsed.values()):
                raise ComparatorError(f"{shard}: a range-keyed map value is not a transcript array")
            parsed = [t for wrapped in parsed.values() for t in wrapped]
        if not isinstance(parsed, list):
            raise ComparatorError(f"{shard} is neither a transcript array nor a range-keyed map")
        transcripts = [t for t in parsed if t is not None]
        if not all(isinstance(t, dict) for t in transcripts):
            raise ComparatorError(f"{shard}: a transcript is not an object")
        return transcripts

    @staticmethod
    def coordinate(shard: Path, transcript: dict, key: str) -> int:
        """The ``key`` coordinate of a transcript, an integer or a string of digits."""
        value = transcript.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, str)) or not str(value).isdecimal():
            raise ComparatorError(f"{shard}: transcript {transcript.get('stable_id')!r} has no integer {key}")
        return int(value)

    def agree(self, expected: str, actual: str) -> bool:
        """Whether two ``NEAREST`` values name the same winning sites: every expected value is
        among the actual values and every actual value starts at a site of an expected value;
        a value the cache lacks starts nowhere, so it never agrees."""
        expected_values = nearest_list(expected)
        actual_values = nearest_list(actual)
        if not expected_values or not actual_values:
            return expected == actual
        if any(e not in actual_values for e in expected_values):
            return False
        expected_sites = set().union(*(self.sites.get(e, set()) for e in expected_values))
        return all(self.sites.get(a, set()) & expected_sites for a in actual_values)


def extra_accepted(expected: str, actual: str, accept: Callable[[str, str, str], bool]) -> bool:
    """Whether ``actual`` is ``expected`` with at least one differing ``key=value`` pair that
    ``accept(key, expected value, actual value)`` allows substituted, the keys and their
    order kept; a difference in the order of the pairs is not accepted."""
    ap = extra_pairs(actual)
    rewritten = []
    substituted = False
    for kv in expected.split(";"):
        k, sep, v = kv.partition("=")
        if sep and k in ap and v != ap[k] and accept(k, v, ap[k]):
            kv = f"{k}={ap[k]}"
            substituted = True
        rewritten.append(kv)
    return substituted and ";".join(rewritten) == actual


def extra_nearest_accepted(expected: str, actual: str, sites: NearestSites) -> bool:
    """``extra_accepted`` for an agreeing ``NEAREST`` value."""
    return extra_accepted(expected, actual, lambda k, v, a: k == "NEAREST" and sites.agree(v, a))


def extra_pairs_accepted(expected: str, actual: str, field_pairs: dict[str, set[tuple[str, str]]]) -> bool:
    """``extra_accepted`` for the accepted field pairs."""
    return extra_accepted(expected, actual, lambda k, v, a: (v, a) in field_pairs.get(k, ()))


def extra_order_only_accepted(expected: str, actual: str, order_only: dict[str, str]) -> bool:
    """``extra_accepted`` for the set-valued columns, their items agreeing as a multiset."""
    return extra_accepted(expected, actual, lambda k, v, a: k in order_only and same_items(v, a, order_only[k]))


def list_items(value: str, separator: str) -> list[str] | None:
    """The items of a list-valued field: the ``separator``-joined parts of a text format's
    value, or the elements of the JSON array ``json_scalar`` joined (``JSON_ARRAY``), each
    as its compact text; ``None`` when the value is not such a list."""
    if separator != JSON_ARRAY:
        return value.split(separator)
    try:
        items = json.loads(f"[{value}]")
    except ValueError:
        return None
    return [json_compact(item) for item in items]


def same_items(expected: str, actual: str, separator: str) -> bool:
    """Whether two list-valued field texts carry the same items the same number of times,
    in any order; an empty or absent value agrees with itself alone."""
    if expected in ("", "-", ABSENT) or actual in ("", "-", ABSENT):
        return expected == actual
    expected_items = list_items(expected, separator)
    actual_items = list_items(actual, separator)
    if expected_items is None or actual_items is None:
        return False
    return Counter(expected_items) == Counter(actual_items)


def order_only_fields(fmt: str) -> dict[str, str]:
    """The fields the ``order_only`` class compares as multisets in ``fmt``, each with its
    list separator: ``,`` in the default and tab formats (``convert_arrayref``,
    Utils.pm:185-198), ``&`` in VCF (VCF.pm:391), and in JSON the arrays of
    ``ORDER_ONLY_JSON_FIELDS``."""
    if fmt == "vcf":
        return dict.fromkeys(ORDER_ONLY_COLUMNS, "&")
    if fmt == "json":
        return dict(ORDER_ONLY_JSON_FIELDS)
    return dict.fromkeys(ORDER_ONLY_COLUMNS, ",")


def json_floor_fields(column: str) -> tuple[str, ...] | None:
    """The JSON fields a floor on the default-format ``column`` is counted under
    (JSON.pm:311-385): the name lower-cased and renamed as ``%RENAME_KEYS`` renames it,
    ``Feature`` under the per-kind id, SIFT and PolyPhen under either of their two keys, a
    position under either of its start and end keys (JSON.pm:340-341 delete the start when
    it is ``?`` and keep the end, which the default format prints as ``?-N``); ``None`` for
    a column with no row count in the JSON format: a record-level column, ``Feature_type``,
    which chooses the consequence list an entry joins and is deleted from the entry
    (JSON.pm:348-350), and ``Extra``, the text formats' container of the keyed fields
    (VEP_output.pm:114-133)."""
    key = column.lower()
    if column in RECORD_LEVEL_COLUMNS or RECORD_LEVEL_COLUMN_RE.search(column) or key in ("feature_type", "extra"):
        return None
    if key == "feature":
        return tuple(f"{kind}_id" for _, kind in JSON_CONSEQUENCE_LISTS if kind != "intergenic")
    pm = POSITION_RE.match(column)
    if pm:
        coord = pm.group(1).lower()
        return (f"{coord}_start", f"{coord}_end")
    if key in ("sift", "polyphen"):
        return (f"{key}_prediction", f"{key}_score")
    return (JSON_RENAME.get(key, key),)


def populated_under(entries: list[Entry], fields: Collection[str]) -> int:
    """The consequence entries (never a record's own entry, which repeats ``strand`` and
    the coordinates) carrying any of ``fields`` with a value other than ``-``: the JSON
    factory drops a key whose default-format value is ``-`` and keeps the allele keys
    (JSON.pm:327), so a present key is a populated column, an empty ``domains`` array
    included (JSON.pm:363-371 keep the key when no item re-splits)."""
    return sum(
        1
        for e in entries
        if not e.key.endswith("|record") and any(f in e.fields and e.fields[f] != "-" for f in fields)
    )


def compare_entries(
    expected: list[Entry],
    actual: list[Entry],
    documented: Documented,
    floors: dict[str, int] | None = None,
    field_pairs: dict[str, set[tuple[str, str]]] | None = None,
    check_order: bool = True,
    examples: int = EXAMPLES_PER_GROUP,
    nearest_sites: NearestSites | None = None,
    order_only: dict[str, str] | None = None,
    floor_fields: dict[str, tuple[str, ...] | None] | None = None,
) -> dict:
    """Compares expected against actual entries and returns the format's report: the
    one-sided key counts, every field with its differing count and examples grouped by the
    expected consequence set, the stale documented divergences, the order of the shared
    keys, the floors, and the ``verdict``. ``order_only`` names the fields compared as
    multisets with their list separators (``order_only_fields``), each agreement counted in
    ``documented.applied`` under the class; ``floor_fields`` maps every floor's column to the
    fields ``populated_under`` counts it by in a format whose keys differ from the column
    names, ``None`` for a column the format has no row count for; a format without the map
    counts the column itself over every entry. ``INDIVIDUAL_ALLELE_ORDER`` in
    ``documented.masks`` compares an individual's alleles as a set: the alternate rows of one
    copy may come in any order within a feature, its allele strings agree on the reference
    allele and the set of alternate alleles, and its ``ALLELE_NUM`` values differ by one
    permutation of the alternates' numbers applied to every row of the copy; each agreement
    is counted under the class, the reordered entries under ``counts``."""
    floors = floors or {}
    field_pairs = field_pairs or {}
    order_only = order_only or {}
    floor_fields = floor_fields or {}
    actual_by_key = {e.key: e for e in actual}
    expected_by_key = {e.key: e for e in expected}

    missing: list[Entry] = []
    unexpected: list[Entry] = []
    mismatches: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    mismatch_totals: Counter[str] = Counter()
    pairs_accepted: Counter[str] = Counter()
    nearest_accepted: Counter[str] = Counter()
    order_only_accepted: Counter[str] = Counter()
    order_only_examples: dict[str, list[str]] = defaultdict(list)
    allele_order = INDIVIDUAL_ALLELE_ORDER in documented.masks
    allele_order_accepted: Counter[str] = Counter()
    allele_order_examples: dict[str, list[str]] = defaultdict(list)
    stale: list[str] = []
    compared = skipped_documented = annotated_without_reference = 0

    def carries_individual(e: Entry) -> bool:
        """Whether the entry itself carries one of ``INDIVIDUAL_FIELDS`` with a value."""
        return any(dash(e.fields.get(f, "")) != "-" for f in INDIVIDUAL_FIELDS)

    individual_records = {
        e.key.split("|", 1)[0] for e in expected if e.key.endswith("|record") and carries_individual(e)
    }

    def is_individual(e: Entry) -> bool:
        """An entry of an individual's copy of a record: it carries IND itself (a default, tab
        or CSQ entry) or its record object carries ind (a JSON consequence entry)."""
        return carries_individual(e) or e.key.split("|", 1)[0] in individual_records

    def copy_of(e: Entry) -> tuple[str, ...]:
        """The copy an entry belongs to: its record label, name and individual."""
        individual = tuple(e.fields.get(f, "") for f in INDIVIDUAL_FIELDS)
        return (e.key.split("|", 1)[0], e.fields.get("Uploaded_variation", ""), *individual)

    def copy_feature(e: Entry) -> tuple[str, ...]:
        """The copy and feature an entry belongs to: ``copy_of`` with its feature and its kind."""
        return (*copy_of(e), e.feature, e.key.rsplit("|", 1)[-1].partition("#")[0])

    def allele_number(e: Entry) -> str | None:
        """The entry's allele number, under either field name, or ``None``."""
        return next((e.fields[name] for name in ALLELE_NUMBER_FIELDS if name in e.fields), None)

    def allele_numbers(entries: list[Entry]) -> dict[tuple[str, ...], dict[str, set[str]]]:
        """Per copy, the allele numbers each of its alleles carries across the copy's rows."""
        numbers: dict[tuple[str, ...], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
        for e in entries:
            if not is_individual(e):
                continue
            number = allele_number(e)
            if number is not None:
                numbers[copy_of(e)][e.allele].add(number)
        return numbers

    expected_numbers = allele_numbers(expected) if allele_order else {}
    actual_numbers = allele_numbers(actual) if allele_order else {}

    def renumbered_alternates(e: Entry, a: Entry, expected_value: str, actual_value: str) -> bool:
        """Whether the row's allele numbers, ``expected_value`` and ``actual_value``, differ by
        one permutation of the alternates' numbers applied to every row of the copy: each
        allele carries one number on each side, the row's own among them, the two sides number
        the same alleles with the same numbers, and the allele numbered 0, the reference, is
        the same allele on both (Parser/VCF.pm:639 keeps the reference first)."""
        expected_map = expected_numbers.get(copy_of(e), {})
        actual_map = actual_numbers.get(copy_of(a), {})
        if not expected_map or expected_map.keys() != actual_map.keys():
            return False
        if any(len(numbers) != 1 for numbers in (*expected_map.values(), *actual_map.values())):
            return False
        expected_number_of = {allele: next(iter(numbers)) for allele, numbers in expected_map.items()}
        actual_number_of = {allele: next(iter(numbers)) for allele, numbers in actual_map.items()}
        if expected_number_of.get(e.allele) != expected_value or actual_number_of.get(a.allele) != actual_value:
            return False
        if sorted(expected_number_of.values()) != sorted(actual_number_of.values()):
            return False
        expected_reference = [allele for allele, number in expected_number_of.items() if number == "0"]
        actual_reference = [allele for allele, number in actual_number_of.items() if number == "0"]
        return expected_reference == actual_reference

    def reference_row(e: Entry) -> bool:
        """Whether the entry is identifiably the copy's reference-allele row: numbered 0, or
        its allele is the reference allele it prints (``-`` for an insertion's); a row
        carrying neither field is not."""
        number = allele_number(e)
        if number is not None:
            return number == "0"
        reference = next((dash(e.fields[name]) for name in REFERENCE_ALLELE_FIELDS if name in e.fields), None)
        return reference is not None and reference == e.allele

    def same_alleles(expected_value: str, actual_value: str) -> bool:
        """The same reference allele followed by the same alternate alleles in any order."""
        if expected_value in ("", "-", ABSENT) or actual_value in ("", "-", ABSENT):
            return expected_value == actual_value
        ref, _, alts = expected_value.partition("/")
        actual_ref, _, actual_alts = actual_value.partition("/")
        return ref == actual_ref and same_items(alts, actual_alts, "/")

    def allele_order_accepts(e: Entry, a: Entry, name: str, ev: str, av: str) -> bool:
        """Whether a field of an individual's entry differs in the order of the copy's
        alleles alone; in the Extra text every differing pair must."""
        if not (allele_order and is_individual(e)):
            return False

        def reordered_only(key: str, expected_value: str, actual_value: str) -> bool:
            if key in ALLELE_STRING_FIELDS:
                return same_alleles(expected_value, actual_value)
            if key in ALLELE_NUMBER_FIELDS:
                return renumbered_alternates(e, a, expected_value, actual_value)
            return key in order_only and same_items(expected_value, actual_value, order_only[key])

        if name == "Extra":
            return extra_accepted(ev, av, reordered_only)
        return name in ALLELE_STRING_FIELDS + ALLELE_NUMBER_FIELDS and reordered_only(name, ev, av)

    def order_key(e: Entry) -> str:
        """The key the order comparison reads: an individual's consequence entry by its copy
        and feature (``copy_feature``) and whether it is the reference row, without its
        allele and duplicate counter, since the copy's alternate rows come in the order of
        its alleles while the reference row keeps its place."""
        if not allele_order or e.key.endswith("|record") or not is_individual(e):
            return e.key
        return "|".join([*copy_feature(e), "reference" if reference_row(e) else "alternate"])

    for e in expected:
        a = actual_by_key.get(e.key)
        sets = documented.divergent(e)
        if sets is not None:
            skipped_documented += 1
            if a is None:
                if sets:
                    stale.append(f"{e.key}: documented candidate sets {sets}, but the candidate emits no row")
            else:
                got = consequence_set(a.consequence)
                if not any(consequence_set(s) == got for s in sets):
                    stale.append(
                        f"{e.key}: documented candidate sets {sets}, the candidate now emits "
                        f"{a.consequence!r} (expected {e.consequence!r})"
                    )
            continue
        if a is None:
            missing.append(e)
            continue
        compared += 1
        record_summary_differs = (
            e.key.endswith("|record") and e.record is not None and documented.record_divergent(e.record)
        )
        # The default format's Extra is the raw text of the fields compared one by one
        # below, so it differs wherever one of them is documented to.
        extra_differs = documented.has_divergent_fields(e)
        for name in sorted(set(e.fields) | set(a.fields)):
            if record_summary_differs and name == "most_severe_consequence":
                continue
            if extra_differs and name == "Extra":
                continue
            ev = e.fields.get(name, ABSENT)
            av = a.fields.get(name, ABSENT)
            documented_value = documented.divergent_field(e, name)
            if documented_value is not None:
                if documented_form(av) != documented_form(documented_value):
                    stale.append(
                        f"{e.key}: field {name} documented as {documented_value!r}, the candidate "
                        f"now prints {av!r} (expected {ev!r})"
                    )
                continue
            if ev == av:
                continue
            if (ev, av) in field_pairs.get(name, ()) or (
                name == "Extra" and extra_pairs_accepted(ev, av, field_pairs)
            ):
                pairs_accepted[name] += 1
                continue
            if nearest_sites is not None and (
                (name in NEAREST_FIELDS and nearest_sites.agree(ev, av))
                or (name == "Extra" and extra_nearest_accepted(ev, av, nearest_sites))
            ):
                nearest_accepted[name] += 1
                continue
            if (name in order_only and same_items(ev, av, order_only[name])) or (
                name == "Extra" and order_only and extra_order_only_accepted(ev, av, order_only)
            ):
                order_only_accepted[name] += 1
                documented.applied[ORDER_ONLY] += 1
                if len(order_only_examples[name]) < examples:
                    order_only_examples[name].append(f"{e.key}: expected {ev!r}, got {av!r}")
                continue
            if allele_order_accepts(e, a, name, ev, av):
                allele_order_accepted[name] += 1
                documented.applied[INDIVIDUAL_ALLELE_ORDER] += 1
                if len(allele_order_examples[name]) < examples:
                    allele_order_examples[name].append(f"{e.key}: expected {ev!r}, got {av!r}")
                continue
            mismatch_totals[name] += 1
            group = mismatches[name][e.consequence]
            if len(group) < examples:
                group.append(f"{e.key}: expected {ev!r}, got {av!r}")
    for a in actual:
        if a.key in expected_by_key or documented.extra(a):
            continue
        # Only the record-level line or object is without a counterpart; a consequence
        # entry on a skipped record is a candidate-only key, documented or unexpected.
        if a.key.endswith("|record") and a.record is not None and documented.reference_skipped(a.record):
            annotated_without_reference += 1
            continue
        unexpected.append(a)

    shared_expected = [e for e in expected if e.key in actual_by_key and documented.divergent(e) is None]
    shared_actual = [a for a in actual if a.key in expected_by_key and documented.divergent(a) is None]
    order = None
    reordered = 0
    if check_order:
        expected_order = [order_key(e) for e in shared_expected]
        actual_order = [order_key(a) for a in shared_actual]
        if expected_order != actual_order:
            first = next((i for i, (x, y) in enumerate(zip(expected_order, actual_order)) if x != y), 0)
            order = {
                "first_divergence": first,
                "expected": shared_expected[first].key if first < len(shared_expected) else None,
                "actual": shared_actual[first].key if first < len(shared_actual) else None,
            }
        else:
            reordered = sum(1 for e, a in zip(shared_expected, shared_actual) if e.key != a.key)
            if reordered:
                documented.applied[INDIVIDUAL_ALLELE_ORDER] += reordered

    expected_populated = populated_counts(expected)
    actual_populated = populated_counts(actual)
    columns: dict[str, dict] = {}
    floors_not_applicable: list[str] = []
    for name in sorted(set(expected_populated) | set(actual_populated) | set(mismatch_totals) | set(floors)):
        expected_count = expected_populated.get(name, 0)
        actual_count = actual_populated.get(name, 0)
        floor = floors.get(name)
        floor_met = None
        floor_detail: dict[str, str | list[str]] = {}
        if floor is not None:
            counted = floor_fields.get(name)
            if name not in floor_fields:
                floor_met = expected_count >= floor and actual_count >= floor
            elif counted is None:
                floors_not_applicable.append(name)
                floor_detail["floor_note"] = FLOOR_NOT_APPLICABLE
            else:
                expected_count = populated_under(expected, counted)
                actual_count = populated_under(actual, counted)
                floor_detail["floor_fields"] = list(counted)
                floor_met = expected_count >= floor and actual_count >= floor
        col = {
            "expected_populated": expected_count,
            "actual_populated": actual_count,
            "differing": mismatch_totals.get(name, 0),
            "pairs_accepted": pairs_accepted.get(name, 0),
            "nearest_accepted": nearest_accepted.get(name, 0),
            "order_only": order_only_accepted.get(name, 0),
            "individual_allele_order": allele_order_accepted.get(name, 0),
            "floor": floor,
            "floor_met": floor_met,
            **floor_detail,
        }
        if name in mismatches:
            col["examples"] = dict(mismatches[name])
        if name in order_only_examples:
            col["order_only_examples"] = order_only_examples[name]
        if name in allele_order_examples:
            col["individual_allele_order_examples"] = allele_order_examples[name]
        columns[name] = col

    reasons: list[str] = []
    if missing:
        reasons.append(f"missing in the candidate: {len(missing)} entries")
    if unexpected:
        reasons.append(f"unexpected in the candidate: {len(unexpected)} entries")
    for name in sorted(mismatch_totals):
        reasons.append(f"field {name}: {mismatch_totals[name]} differing entries")
    if stale:
        reasons.append(f"documented divergences not observed: {len(stale)}")
    if order is not None:
        reasons.append(f"order differs at shared position {order['first_divergence']}")
    for name in sorted(floors):
        col = columns[name]
        if col["floor_met"] is False:
            reasons.append(
                f"floor {name}={floors[name]}: expected populated {col['expected_populated']}, "
                f"actual populated {col['actual_populated']}"
            )

    def grouped(entries: list[Entry]) -> dict[str, dict]:
        by_cons: dict[str, list[Entry]] = defaultdict(list)
        for e in entries:
            by_cons[e.consequence].append(e)
        return {
            cons: {"count": len(items), "examples": [e.key for e in items[:examples]]}
            for cons, items in sorted(by_cons.items())
        }

    return {
        "verdict": "red" if reasons else "green",
        "reasons": reasons,
        "counts": {
            "expected_entries": len(expected),
            "actual_entries": len(actual),
            "compared": compared,
            "documented_divergences_skipped": skipped_documented,
            "record_entries_on_reference_skipped_records": annotated_without_reference,
            "individual_entries_reordered": reordered,
        },
        "one_sided": {"expected_only": len(missing), "actual_only": len(unexpected)},
        "missing": grouped(missing),
        "unexpected": grouped(unexpected),
        "columns": columns,
        "floors_not_applicable": floors_not_applicable,
        "stale_divergences": stale[:20],
        "order": order,
    }


def is_regulatory_row(e: Entry) -> bool:
    """A RegulatoryFeature or MotifFeature consequence entry in any format: the key's last
    segment is the Feature_type (default, tab, VCF) or the JSON list's kind."""
    last_segment = e.key.rsplit("|", 1)[-1].partition("#")[0]
    return feature_kind(last_segment) in REGULATORY_KINDS


def filter_consequence_rows(
    entries: list[Entry], keep: Callable[[Entry], bool], keep_bare_records: bool
) -> list[Entry]:
    """The entries whose consequence rows ``keep`` accepts, in their order, with the VCF/JSON
    record entries: every one when ``keep_bare_records``, else those of a record that kept a
    row. A record that lost a row has its ``most_severe_consequence`` recomputed from the
    rows it kept, as the candidate computes it from the rows it writes. Mutates the kept
    record entries; raises ``ComparatorError`` on a term ``CONSEQUENCE_RANK`` lacks."""
    terms_by_record: dict[str, set[str]] = defaultdict(set)
    lost_rows: set[str] = set()
    for e in entries:
        if e.key.endswith("|record"):
            continue
        record = e.key.split("|", 1)[0]
        if keep(e):
            terms_by_record[record].update(t for t in re.split(r"[&,]", e.consequence) if t)
        else:
            lost_rows.add(record)
    kept: list[Entry] = []
    for e in entries:
        if not e.key.endswith("|record"):
            if keep(e):
                kept.append(e)
            continue
        record = e.key.split("|", 1)[0]
        if not keep_bare_records and record not in terms_by_record:
            continue
        if record in lost_rows and "most_severe_consequence" in e.fields:
            terms = terms_by_record.get(record, set())
            unranked = sorted(terms.difference(CONSEQUENCE_RANK))
            if unranked:
                raise ComparatorError(f"{e.key}: no rank for the consequence terms {unranked}")
            e.fields = {**e.fields, "most_severe_consequence": most_severe(terms)}
        kept.append(e)
    return kept


def without_regulatory_rows(expected: list[Entry]) -> list[Entry]:
    """The reference entries without their regulatory rows (``regulatory_feature_rows``): a
    record the reference annotates against a regulatory feature and no transcript compares
    as its intergenic row, every record entry stays (``filter_consequence_rows``)."""
    return filter_consequence_rows(expected, lambda e: not is_regulatory_row(e), True)


def only_regulatory_rows(entries: list[Entry]) -> list[Entry]:
    """The regulatory rows of either side and the record entries of the records that have
    one (``--keep-only regulatory_feature_rows``; ``filter_consequence_rows``)."""
    return filter_consequence_rows(entries, is_regulatory_row, False)


def compare_files(
    fmt: str,
    expected_path: Path,
    actual_path: Path,
    inputs: list[str] | None,
    documented: Documented,
    masks: set[str],
    floors: dict[str, int],
    field_pairs: dict[str, set[tuple[str, str]]],
    examples: int,
    nearest_sites: NearestSites | None = None,
    unordered: bool = False,
    keep_only: str | None = None,
) -> dict:
    """The comparison-mode report for one format. ``unordered`` compares the shared entries
    without their order (a run whose reference rows Ensembl VEP emits in hash order); a
    Parquet actual is compared without order in any case."""
    if keep_only is not None and keep_only != REGULATORY_FEATURE_ROWS:
        raise ComparatorError(f"unknown --keep-only class: {keep_only}")
    if keep_only == REGULATORY_FEATURE_ROWS and REGULATORY_FEATURE_ROWS in masks:
        raise ComparatorError("--keep-only regulatory_feature_rows conflicts with --masks regulatory_feature_rows")
    documented.masks |= masks
    expected_format = "tab" if fmt == "parquet" else fmt
    expected = parse_format(expected_format, expected_path, inputs)
    actual = parse_format(fmt, actual_path, inputs)
    if REGULATORY_FEATURE_ROWS in masks:
        kept = without_regulatory_rows(expected)
        # A class that dropped nothing is unused, as a manifest class matching no item is.
        if len(kept) < len(expected):
            documented.applied[REGULATORY_FEATURE_ROWS] += len(expected) - len(kept)
        expected = kept
    if keep_only == REGULATORY_FEATURE_ROWS:
        expected = only_regulatory_rows(expected)
        actual = only_regulatory_rows(actual)
    report = compare_entries(
        expected,
        actual,
        documented,
        floors=floors,
        field_pairs=field_pairs,
        check_order=fmt != "parquet" and not unordered,
        examples=examples,
        nearest_sites=nearest_sites,
        order_only=order_only_fields(fmt) if ORDER_ONLY in masks else None,
        floor_fields={name: json_floor_fields(name) for name in floors} if fmt == "json" else None,
    )
    report.update(
        {
            "mode": "compare",
            "format": fmt,
            "expected": str(expected_path),
            "actual": str(actual_path),
            "masks": {
                "requested": sorted(masks),
                "applied": dict(sorted(documented.applied.items())),
                "unused": sorted(masks - set(documented.applied)),
                "available": dict(sorted(documented.available.items())),
            },
            "field_pairs": {name: sorted(pairs) for name, pairs in sorted(field_pairs.items())},
            "nearest_sites": None if nearest_sites is None else len(nearest_sites.sites),
            "unordered": unordered,
            "keep_only": keep_only,
        }
    )
    return report


POSITION_RE = re.compile(r"(\w+?)_position$", re.IGNORECASE)
POSITION_PART_RE = re.compile(r"^\d+(/\d)?\d*$")
TOOL_SCORE_RE = re.compile(r"([a-z_]+)?\(?([\d.]+)?\)?", re.IGNORECASE)
DOMAIN_RE = re.compile(r"(\w+):(\w+)")
PERL_NUMBER_RE = re.compile(r"^\s*[+-]?(\d+\.?\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?)")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def perl_numify(s: str) -> str:
    """``0 + $s``: the leading number of a string, 0 when it has none."""
    m = PERL_NUMBER_RE.match(s)
    if not m:
        return "0"
    text = m.group(1)
    try:
        value = float(text)
    except ValueError:
        return "0"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    if value == int(value) and "e" not in text.lower():
        return str(int(value))
    return _json_float(value)


def perl_numeric_or_text(word: str) -> int | str:
    """A ``\\w+`` capture as ``numberify`` leaves it: ASCII digits alone become an integer (Perl's
    ``+ 0`` on a string ``looks_like_number`` accepts); anything else stays text."""
    return int(word) if word.isascii() and word.isdigit() else word


def format_coords(start, end) -> str:
    """``Utils.pm:141-166``: ``start``, ``start-end``, or the pair low to high."""
    if start is not None and end is not None:
        try:
            s, e = float(start), float(end)
        except (TypeError, ValueError):
            return f"{start}-{end}"
        if s > e:
            return f"{end}-{start}"
        if s == e:
            return f"{start}"
        return f"{start}-{end}"
    if start is not None:
        return f"{start}-?"
    if end is not None:
        return f"?-{end}"
    return "-"


def uploaded_variation(obj: dict) -> str:
    """``Uploaded_variation`` for a record of VCF input as ``OutputFactory.pm:883-930``
    derives it from the JSON object's ``id``, ``seq_region_name``, ``start``,
    ``allele_string`` and ``input``: the ID when it is one, else ``chr_pos_ref/alt``."""
    chrom = json_scalar(obj.get("seq_region_name", ""))
    start = obj.get("start")
    allele_string = json_scalar(obj.get("allele_string", ""))
    try:
        var_name = f"{chrom}_{int(float(start))}_{allele_string}"
    except (TypeError, ValueError):
        var_name = f"{chrom}_{start}_{allele_string}"
    vn = obj.get("id")
    vn = json_scalar(vn) if vn is not None else None
    if vn not in (None, "", ".") and vn != var_name:
        return vn
    line = obj.get("input")
    cols = line.split("\t") if isinstance(line, str) else []
    if len(cols) >= 5:
        chrom_in, pos, ref, alt = cols[0], cols[1], cols[3], cols[4]
        alt = alt.replace(",", "/")
        return "_".join([chrom_in, perl_numify(pos), f"{ref}/{alt}" if alt != "" else ref])
    if vn not in (None, ".") and vn != var_name:
        return vn
    return f"{chrom}_{start}_{allele_string}"


def record_location(obj: dict) -> str:
    """``Location`` as ``OutputFactory.pm:931`` prints it from the JSON object's
    ``seq_region_name``, ``start`` and ``end``."""
    return f"{json_scalar(obj.get('seq_region_name', ''))}:{format_coords(obj.get('start'), obj.get('end'))}"


def feature_kind(feature_type: str) -> str:
    """``JSON.pm:348-349``: the consequence list a feature type's entries go to
    (``Transcript`` to ``transcript``, ``RegulatoryFeature`` to ``regulatory_feature``, an
    undefined type to ``intergenic``); a list name is returned unchanged."""
    kind = (feature_type if feature_type not in ("", "-") else "intergenic").lower()
    if kind == "intergenic" or "_feature" in kind:
        return kind
    return kind.replace("feature", "_feature")


def vcf_escape(s: str) -> str:
    """``VCF.pm:401-404`` on a value already joined for the default format."""
    s = s.replace(",", "&").replace(";", "%3B")
    s = re.sub(r"\s+", "_", s)
    return s.replace("|", "&")


def reverse_complement(allele: str) -> str:
    """The reverse complement of a sequence allele; any other allele unchanged."""
    if allele and re.fullmatch(r"[ACGTNacgtn-]+", allele):
        return allele.translate(COMPLEMENT)[::-1]
    return allele


def values_equal(a: str, b: str) -> bool:
    """Equality of two scalar texts, numerically when both are numbers."""
    if a == b:
        return True
    try:
        return float(a) == float(b)
    except ValueError:
        return False


def lists_equal(a: list[str], b: list[str]) -> bool:
    return len(a) == len(b) and all(values_equal(x, y) for x, y in zip(a, b))


def default_list(value: str) -> list[str] | None:
    """A default-format value as the list ``convert_arrayref`` joined, percent-encoding
    undone; ``None`` for ``-``."""
    if value in ("-", ABSENT):
        return None
    return decode_extra(value).split(",")


def json_list(value) -> list[str] | None:
    """A JSON value as a list of scalar texts, a string split on ``,`` as the default
    format's list would be; ``None`` when absent."""
    if value is None:
        return None
    if isinstance(value, list):
        return [json_scalar(x) for x in value]
    return json_scalar(value).split(",")


def most_severe(terms: set[str]) -> str:
    """``JSON.pm:387-388``: the lowest-ranked term, ``?`` when there is none."""
    ranked = sorted((t for t in terms if t != "?"), key=lambda t: (CONSEQUENCE_RANK.get(t, 10**6), t))
    return ranked[0] if ranked else "?"


class SelfParity:
    """One engine's outputs on one input, every consequence entry keyed by the record's
    ``Uploaded_variation`` and ``Location``, its allele, its feature and its feature kind,
    and compared column by column against the default format."""

    def __init__(self, default_text: str, tab_text: str, vcf_text: str, json_text: str,
                 parquet_text: str | None = None, inputs: list[str] | None = None,
                 examples: int = EXAMPLES_PER_GROUP) -> None:
        self.examples = examples
        objects = json_objects(json_text)
        if inputs is None:
            inputs = [obj["input"] for _, obj in objects if isinstance(obj.get("input"), str)]
        index = InputIndex(inputs, vcf_line_key)
        self.records: dict[str, dict] = {}
        self.objects_by_name: dict[tuple[str, str], list[dict]] = defaultdict(list)
        self.parts: dict[str, tuple[str, str, str]] = {}
        for n, obj in objects:
            _, label = json_record_label(index, n, obj)
            self.records[label] = obj
            self.objects_by_name[(uploaded_variation(obj), record_location(obj))].append(obj)
        self.default = self._rekey_rows(parse_default(default_text))
        self.tab = self._rekey_rows(parse_tab(tab_text))
        self.parquet = self._rekey_rows(parse_tab(parquet_text)) if parquet_text is not None else None
        vcf_entries = parse_vcf(vcf_text, inputs)
        self.vcf_records = [e for e in vcf_entries if e.key.endswith("|record")]
        self.vcf = self._rekey_by_record(vcf_entries)
        self.json = self._rekey_by_record(parse_json(json_text, inputs))
        self.tab_columns = header_columns(tab_text)
        self.parquet_columns = header_columns(parquet_text) if parquet_text is not None else []
        self.csq_columns = self._csq_columns(vcf_text)
        self.default_columns: set[str] = set()
        for e in self.default:
            self.default_columns.update(e.fields)

    @staticmethod
    def _csq_columns(text: str) -> list[str]:
        for line in text_lines(text):
            if line.startswith("##INFO=<ID=CSQ"):
                return csq_format_fields(line)
        return []

    def _rekeyed(self, keyed: list[tuple[Entry, str, str, str]]) -> list[Entry]:
        """Entries keyed by (name, Location, allele, feature, kind), duplicates numbered,
        each key's parts kept in ``self.parts`` (a name may itself carry ``|``)."""
        out = [
            Entry(
                key=f"{name}|{location}|{e.allele}|{dash(e.feature)}|{kind}",
                record=e.record,
                location=e.location,
                allele=e.allele,
                feature=e.feature,
                consequence=e.consequence,
                fields=e.fields,
            )
            for e, name, location, kind in keyed
        ]
        dedup_keys(out)
        for entry, (_, name, location, kind) in zip(out, keyed):
            self.parts[entry.key] = (name, location, kind)
        return out

    def _rekey_rows(self, rows: list[Entry]) -> list[Entry]:
        """Default-format rows keyed by their own ``Uploaded_variation`` and ``Location``."""
        return self._rekeyed(
            [
                (
                    row,
                    row.fields.get("Uploaded_variation", ""),
                    row.fields.get("Location", ""),
                    feature_kind(row.fields.get("Feature_type", "")),
                )
                for row in rows
            ]
        )

    def _rekey_by_record(self, entries: list[Entry]) -> list[Entry]:
        """VCF and JSON consequence entries keyed by the name and Location of the JSON
        object of their record; an entry whose record has no object keeps its label."""
        keyed = []
        for e in entries:
            if e.key.endswith("|record"):
                continue
            label, _, _, kind = e.key.split("|", 3)
            kind = feature_kind(kind.partition("#")[0])
            obj = self.records.get(label)
            if obj is None:
                keyed.append((e, label, "", kind))
            else:
                keyed.append((e, uploaded_variation(obj), record_location(obj), kind))
        return self._rekeyed(keyed)

    @staticmethod
    def _by_key(entries: list[Entry]) -> dict[str, Entry]:
        return {e.key: e for e in entries}

    def _pair(self, name: str, pivot: dict[str, Entry], other: dict[str, Entry], compare_fields) -> dict:
        pivot_only = sorted(set(pivot) - set(other))
        other_only = sorted(set(other) - set(pivot))
        shared = sorted(set(pivot) & set(other))
        differing: dict[str, list[str]] = defaultdict(list)
        totals: Counter[str] = Counter()
        compared_columns: set[str] = set()
        for key in shared:
            for column, detail in compare_fields(pivot[key], other[key], compared_columns):
                totals[column] += 1
                if len(differing[column]) < self.examples:
                    differing[column].append(f"{key}: {detail}")
        reasons = []
        if pivot_only:
            reasons.append(f"{name}: {len(pivot_only)} entries of the default format are absent")
        if other_only:
            reasons.append(f"{name}: {len(other_only)} entries have no default-format row")
        for column in sorted(totals):
            reasons.append(f"{name}: column {column} differs on {totals[column]} entries")
        return {
            "shared_entries": len(shared),
            "one_sided": {"default_only": len(pivot_only), f"{name}_only": len(other_only)},
            "default_only_examples": pivot_only[: self.examples],
            f"{name}_only_examples": other_only[: self.examples],
            "columns_compared": sorted(compared_columns),
            "differing": {c: {"count": totals[c], "examples": ex} for c, ex in sorted(differing.items())},
            "reasons": reasons,
        }

    def _compare_tab_like(self, columns: list[str]):
        declared = [c for c in columns if c != "Extra"]

        def compare(p: Entry, o: Entry, compared: set[str]):
            for column in declared:
                compared.add(column)
                raw = p.fields.get(column, ABSENT)
                pv = raw.replace("%3B", ";")
                ov = o.fields.get(column, ABSENT)
                if pv == ABSENT and ov == "-":
                    continue
                if pv != ov:
                    yield column, f"default {raw!r}, other {ov!r}"

        return compare

    def _compare_vcf(self):
        columns = [c for c in self.csq_columns if c not in VCF_IMPLIED_COLUMNS or c in self.default_columns]

        def compare(p: Entry, o: Entry, compared: set[str]):
            strand = self._object_for(p).get("strand")
            for column in columns:
                compared.add(column)
                pv = p.fields.get(column, ABSENT)
                ov = o.fields.get(column, "")
                if column == "Allele":
                    expected = pv if pv != ABSENT else "-"
                    if strand == -1:
                        expected = reverse_complement(expected)
                    expected = vcf_escape(expected.replace("%3B", ";"))
                elif pv in (ABSENT, "-"):
                    expected = ""
                else:
                    expected = vcf_escape(pv.replace("%3B", ";"))
                if expected != ov:
                    yield column, f"default {pv!r} (as CSQ {expected!r}), vcf {ov!r}"

        return compare

    def _compare_json(self):
        def scalar_differs(want: str | None, got: str | None) -> bool:
            if want is None or got is None:
                return (want is None) != (got is None)
            return not values_equal(want, got)

        def compare(p: Entry, o: Entry, compared: set[str]):
            kind = self.parts[p.key][2]
            remaining = set(o.fields)
            for column, raw in p.fields.items():
                if column in ("Extra", "Feature_type") or column in RECORD_LEVEL_COLUMNS:
                    continue
                if RECORD_LEVEL_COLUMN_RE.search(column):
                    continue
                compared.add(column)
                key = column.lower()
                pm = POSITION_RE.match(column)
                if pm:
                    coord = pm.group(1).lower()
                    start_key, end_key = f"{coord}_start", f"{coord}_end"
                    remaining.discard(start_key)
                    remaining.discard(end_key)
                    want_start = want_end = None
                    if raw not in ("-", ""):
                        parts = raw.split("-")
                        start_text = parts[0]
                        end_text = parts[1] if len(parts) > 1 and parts[1] != "" else None
                        want_start = start_text if POSITION_PART_RE.match(start_text) else None
                        if end_text is not None and (end_text.isdigit() or "?" in end_text):
                            want_end = end_text
                        else:
                            want_end = start_text
                        if end_text is not None and not POSITION_PART_RE.match(want_end):
                            want_end = None
                    for jkey, want in ((start_key, want_start), (end_key, want_end)):
                        got = o.fields.get(jkey)
                        if scalar_differs(want, got):
                            yield column, f"default {raw!r}, json {jkey}={got!r}"
                    continue
                if key in ("sift", "polyphen"):
                    pred_key, score_key = f"{key}_prediction", f"{key}_score"
                    remaining.discard(pred_key)
                    remaining.discard(score_key)
                    want_pred = want_score = None
                    if raw not in ("-", ""):
                        tm = TOOL_SCORE_RE.search(raw)
                        if tm:
                            want_pred, want_score = tm.group(1) or None, tm.group(2)
                    for jkey, want in ((pred_key, want_pred), (score_key, want_score)):
                        got = o.fields.get(jkey)
                        if scalar_differs(want, got):
                            yield column, f"default {raw!r}, json {jkey}={got!r}"
                    continue
                if key == "domains":
                    remaining.discard("domains")
                    want = None
                    if raw not in ("-", ""):
                        # OutputFactory/JSON.pm:363-371 re-splits each entry with ``(\w+):(\w+)`` and
                        # pushes ``{db => $1, name => $2}`` if ``$1 && $2``. A failed match leaves ``$1``
                        # and ``$2`` at the list's last successful one, so such an entry repeats the
                        # previous object and is dropped only before the first match. A capture of ``0``
                        # is false to Perl, so its entry and every unmatched one after it push nothing.
                        # ``numberify`` then runs over the record (JSON.pm:288), so a name of digits
                        # alone is a number.
                        want = []
                        db = name = None
                        perl_false = (None, "", "0")
                        for item in decode_extra(raw).split(","):
                            dm = DOMAIN_RE.search(item)
                            if dm:
                                db, name = dm.group(1), dm.group(2)
                            if db not in perl_false and name not in perl_false:
                                domain = {"db": perl_numeric_or_text(db), "name": perl_numeric_or_text(name)}
                                want.append(json_compact(domain))
                    # json_scalar joins the list's compact objects with ",", which each object
                    # also contains, so the comparison is on the joined text.
                    got = o.fields.get("domains")
                    if (want is None) != (got is None) or (want is not None and ",".join(want) != got):
                        yield column, f"default {raw!r}, json domains={got!r}"
                    continue
                if column in JSON_OBJECT_FIELDS:
                    remaining.discard(key)
                    want = None if raw in ("-", "") else numberified_json_text(raw)
                    got = o.fields.get(key)
                    if want != got:
                        yield column, f"default {raw!r}, json {key}={got!r}"
                    continue
                jkey = f"{kind}_id" if key == "feature" else JSON_RENAME.get(key, key)
                remaining.discard(jkey)
                want = default_list(raw)
                if want is not None:
                    want = ["1" if x == "YES" else x for x in want]
                if key in ("allele", "used_ref", "given_ref") and want is None:
                    want = ["-"]
                got = o.fields.get(jkey)
                got_list = got.split(",") if got is not None else None
                if (want is None) != (got_list is None) or (want is not None and not lists_equal(want, got_list)):
                    yield column, f"default {raw!r}, json {jkey}={got!r}"
            for jkey in sorted(remaining):
                compared.add(jkey)
                yield jkey, f"json carries {jkey}={o.fields[jkey]!r} and the default format no column for it"

        return compare

    def _object_for(self, e: Entry) -> dict:
        """The JSON object whose name and Location an entry's key carries, else empty."""
        name, location, _ = self.parts[e.key]
        return self.objects_by_name.get((name, location), [{}])[0]

    def _record_checks(self) -> dict:
        """Record-level values the JSON format carries once: ``most_severe_consequence``
        against the rank order of the object's own terms (``JSON.pm:387-388``),
        ``colocated_variants`` ids against the rows' ``Existing_variation``, the ``zyg``
        value against the rows' ``ZYG`` (the list of ``--individual_zyg`` as a multiset, since
        its order carries no meaning; the one zygosity of ``--individual`` as itself), and the
        VCF ``ID`` column against the record's name."""
        problems: list[str] = []
        existing_by_name: dict[tuple[str, str], set[str]] = defaultdict(set)
        zyg_by_name: dict[tuple[str, str], set[str]] = defaultdict(set)
        for row in self.default:
            name, location, _ = self.parts[row.key]
            existing_by_name[(name, location)].add(row.fields.get("Existing_variation", "-"))
            zyg = row.fields.get("ZYG", "-")
            if zyg != "-":
                zyg_by_name[(name, location)].add(",".join(sorted(zyg.split(","))))
        for (name, location), objects in sorted(self.objects_by_name.items()):
            want_ids: set[str] = set()
            for value in existing_by_name.get((name, location), set()):
                want_ids.add(",".join(sorted(default_list(value) or [])))
            got_ids: set[str] = set()
            got_zyg: set[str] = set()
            for obj in objects:
                terms = set()
                for list_name, _ in JSON_CONSEQUENCE_LISTS:
                    for item in obj.get(list_name) or []:
                        if isinstance(item, dict):
                            terms.update(json_list(item.get("consequence_terms")) or [])
                # An object with no consequence has no row to hold its zyg to: the
                # reference-homozygous copy of --individual (OutputFactory.pm:511 takes the
                # alternate alleles alone without --process_ref_homs) and a record carried
                # with no allele to annotate.
                if "zyg" in obj and terms:
                    got_zyg.add(",".join(sorted(json_scalar(obj["zyg"]).split(","))))
                want = most_severe(terms)
                got = obj.get("most_severe_consequence")
                if got is None and not terms and "id" not in obj:
                    continue
                if got != want:
                    problems.append(
                        f"{name}|{location}: most_severe_consequence {got!r}, the object's terms give {want!r}"
                    )
                colocated = obj.get("colocated_variants")
                ids: list[str] = []
                if isinstance(colocated, list):
                    ids = [json_scalar(c.get("id", "")) for c in colocated if isinstance(c, dict)]
                got_ids.add(",".join(sorted(ids)))
            if (name, location) in existing_by_name and want_ids != got_ids:
                problems.append(
                    f"{name}|{location}: Existing_variation {sorted(want_ids)} "
                    f"against colocated_variants ids {sorted(got_ids)}"
                )
            want_zyg = zyg_by_name.get((name, location), set())
            if (name, location) in existing_by_name and want_zyg != got_zyg:
                problems.append(f"{name}|{location}: ZYG {sorted(want_zyg)} against the objects' zyg {sorted(got_zyg)}")
        for e in self.vcf_records:
            label = e.key.split("|", 1)[0]
            obj = self.records.get(label)
            vcf_id = e.fields.get("ID", ".")
            if obj is None or vcf_id == "." or "id" not in obj:
                continue
            name = uploaded_variation(obj)
            if name != vcf_id:
                problems.append(f"{label}: VCF ID {vcf_id!r} against the name {name!r}")
        return {"problems": problems}

    def run(self) -> dict:
        pivot = self._by_key(self.default)
        pairs = {
            "tab": self._pair("tab", pivot, self._by_key(self.tab), self._compare_tab_like(self.tab_columns)),
            "vcf": self._pair("vcf", pivot, self._by_key(self.vcf), self._compare_vcf()),
            "json": self._pair("json", pivot, self._by_key(self.json), self._compare_json()),
        }
        if self.parquet is not None:
            pairs["parquet"] = self._pair(
                "parquet", pivot, self._by_key(self.parquet), self._compare_tab_like(self.parquet_columns)
            )
        records = self._record_checks()
        reasons = [r for pair in pairs.values() for r in pair["reasons"]] + records["problems"]
        keyed_or_record_checked = {"Uploaded_variation", "Location", "Existing_variation", "ZYG"}
        not_compared = {
            "default_record_level": sorted(
                c for c in self.default_columns
                if (c in RECORD_LEVEL_COLUMNS and c not in keyed_or_record_checked)
                or RECORD_LEVEL_COLUMN_RE.search(c)
            ),
            "vcf_implied": sorted(
                c for c in self.csq_columns if c in VCF_IMPLIED_COLUMNS and c not in self.default_columns
            ),
            "tab_undeclared": sorted(c for c in self.default_columns if c != "Extra" and c not in self.tab_columns),
            "vcf_undeclared": sorted(
                c for c in self.default_columns
                if c != "Extra" and c not in self.csq_columns and c not in RECORD_LEVEL_COLUMNS
            ),
        }
        return {
            "mode": "self-parity",
            "verdict": "red" if reasons else "green",
            "reasons": reasons,
            "records": len(self.records),
            "default_entries": len(pivot),
            "pairs": pairs,
            "record_checks": records,
            "not_compared": not_compared,
        }


WARNING_LINE_RE = re.compile(r"^WARNING:\s*(?:line (\d+) skipped|variant skipped)\s*(?:\((.*?)\): )?(.*)$")
WARNING_OTHER_RE = re.compile(r"^WARNING:\s*(.*?)(?:\s+(?:on )?line (\d+))?$")
SKIPPED_LINE_RE = re.compile(r"^Line\s+(\S+)\s*\t(.*)\t(.*)$")
PATH_RE = re.compile(r"(?<!\S)(?:~|\.{1,2})?/[^\s]+")


def mask_description(text: str) -> str:
    """A warning description with every path token replaced by ``<path>``."""
    return PATH_RE.sub("<path>", text).strip()


def parse_warning_lines(text: str) -> set[tuple[int | None, str]]:
    """The (input line number, description) set of a ``_warnings.txt`` or a
    ``--skipped_variants_file`` output; ``BaseVEP.pm:802-863`` and
    ``BaseRunner.pm:364-380`` are the two line shapes."""
    out: set[tuple[int | None, str]] = set()
    for line in text_lines(text):
        if not line.strip() or line.startswith("[VEP skipped the following variants from"):
            continue
        m = SKIPPED_LINE_RE.match(line)
        if m:
            number = int(m.group(1)) if m.group(1).isdigit() else None
            out.add((number, mask_description(m.group(3))))
            continue
        m = WARNING_LINE_RE.match(line)
        if m and (m.group(1) is not None or line.startswith("WARNING: variant skipped")):
            number = int(m.group(1)) if m.group(1) is not None else None
            out.add((number, mask_description(m.group(3))))
            continue
        m = WARNING_OTHER_RE.match(line)
        if m:
            number = int(m.group(2)) if m.group(2) is not None else None
            out.add((number, mask_description(m.group(1))))
            continue
        out.add((None, mask_description(line)))
    return out


def compare_warning_files(expected_path: Path, actual_path: Path, masks: Collection[str] = ()) -> dict:
    """The two files as one-sided entry sets. Under ``regulatory_feature_rows`` the candidate's
    notice that it annotates no regulatory feature is required and left out of the sets."""
    expected = parse_warning_lines(read_text(expected_path))
    actual = parse_warning_lines(read_text(actual_path))
    notice_missing = False
    if REGULATORY_FEATURE_ROWS in masks:
        notice = {
            entry for entry in actual
            if REGULATORY_NOTICE in entry[1] or MISSING_REGULATORY_NOTICE in entry[1]
        }
        notice_missing = not notice
        actual -= notice

    def rows(items):
        numbered_first = sorted(items, key=lambda x: (x[0] is None, x[0] or 0, x[1]))
        return [{"line": n, "description": d} for n, d in numbered_first]

    return {
        "expected": str(expected_path),
        "actual": str(actual_path),
        "expected_entries": len(expected),
        "actual_entries": len(actual),
        "shared": len(expected & actual),
        "expected_only": rows(expected - actual),
        "actual_only": rows(actual - expected),
        "regulatory_notice_missing": notice_missing,
    }


def compare_warnings(
    expected: Path, actual: Path, expected_skipped: Path | None, actual_skipped: Path | None, masks: Collection[str] = ()
) -> dict:
    """The warnings-mode report: the ``_warnings.txt`` pair and, when given, the
    skipped-variants pair, each as one-sided sets; ``masks`` applies to the warnings pair."""
    sections = {"warnings": compare_warning_files(expected, actual, masks)}
    if expected_skipped is not None or actual_skipped is not None:
        if expected_skipped is None or actual_skipped is None:
            raise ComparatorError("--expected-skipped and --actual-skipped go together")
        sections["skipped"] = compare_warning_files(expected_skipped, actual_skipped)
    reasons = []
    for name, section in sections.items():
        if section["expected_only"]:
            reasons.append(f"{name}: {len(section['expected_only'])} entries only in the expected file")
        if section["actual_only"]:
            reasons.append(f"{name}: {len(section['actual_only'])} entries only in the actual file")
        if section["regulatory_notice_missing"]:
            reasons.append(f"{name}: the regulatory notice is missing from the actual file")
    return {
        "mode": "warnings",
        "format": "warnings",
        "verdict": "red" if reasons else "green",
        "reasons": reasons,
        **sections,
    }


def parse_floors(spec: str | None) -> dict[str, int]:
    """``COL=MIN,...`` into a column-to-minimum map."""
    floors: dict[str, int] = {}
    if not spec:
        return floors
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        name, sep, value = item.partition("=")
        if not sep or not name or not value.isdigit():
            raise ComparatorError(f"--floors expects COL=MIN with a non-negative integer, got {item!r}")
        floors[name] = int(value)
    return floors


def parse_field_pairs(specs: list[str] | None) -> dict[str, set[tuple[str, str]]]:
    """``COL=EXPECTED=>ACTUAL`` items, several per argument separated by ``;``."""
    pairs: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for spec in specs or []:
        for item in spec.split(";"):
            if not item.strip():
                continue
            name, sep, rest = item.partition("=")
            expected, arrow, actual = rest.partition("=>")
            if not sep or not name or not arrow:
                raise ComparatorError(f"--field-pairs expects COL=EXPECTED=>ACTUAL, got {item!r}")
            pairs[name].add((expected, actual))
    return dict(pairs)


def parse_masks(spec: str | None) -> set[str]:
    return {m.strip() for m in (spec or "").split(",") if m.strip()}


def load_manifest(path: Path | None) -> dict | None:
    if path is None:
        return None
    try:
        manifest = json.loads(read_text(path))
    except ValueError as exc:
        raise ComparatorError(f"{path} is not JSON: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ComparatorError(f"{path} is not a manifest object")
    return manifest


def write_report(report: dict, path: Path | None) -> None:
    if path is None:
        return
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    except OSError as exc:
        raise ComparatorError(f"cannot write {path}: {exc.strerror or exc}") from exc


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="compare_vep_formats.py",
        description="Per-column comparison of VEP output files; see the module docstring for the rules.",
    )
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--self-parity", action="store_true", help="one engine's formats against each other")
    mode.add_argument("--warnings", action="store_true", help="two engines' warnings files against each other")
    ap.add_argument("--format", choices=FORMATS, help="the format of --expected and --actual")
    ap.add_argument("--expected", type=Path, help="the Ensembl VEP output, or its warnings file")
    ap.add_argument("--actual", type=Path, help="the candidate output, or its warnings file")
    ap.add_argument("--input", type=Path, help="the input VCF (needed for --format vcf and json)")
    ap.add_argument("--manifest", type=Path, help="the corpus manifest.json documenting divergences")
    ap.add_argument(
        "--masks",
        help="comma-separated intended-divergence classes applied: the manifest's, and the built-in "
        "regulatory_feature_rows (the one class --warnings mode takes), order_only (DOMAINS, "
        "FREQS and ZYG compared as multisets) and individual_allele_order (an individual's "
        "alleles compared as a set)",
    )
    ap.add_argument(
        "--keep-only", choices=(REGULATORY_FEATURE_ROWS,),
        help="comparison mode: retain only RegulatoryFeature/MotifFeature rows and their records on both sides",
    )
    ap.add_argument(
        "--floors",
        help="COL=MIN,... minimum populated counts per default-format column, both sides; in JSON "
        "counted under the column's JSON key, a column with no row count there (record-level, "
        "Feature_type, Extra) reported as not applicable",
    )
    ap.add_argument("--field-pairs", action="append", help="COL=EXPECTED=>ACTUAL[;...] value pairs accepted as equal")
    ap.add_argument(
        "--nearest-sites", type=Path, help="the candidate JSON cache whose start sites settle NEAREST values"
    )
    ap.add_argument("--default", type=Path, help="self-parity: the default-format output")
    ap.add_argument("--tab", type=Path, help="self-parity: the --tab output")
    ap.add_argument("--vcf", type=Path, help="self-parity: the --vcf output")
    ap.add_argument("--json", type=Path, help="self-parity: the --json output")
    ap.add_argument("--parquet", type=Path, help="self-parity: the Parquet output (file or directory)")
    ap.add_argument("--expected-skipped", type=Path, help="warnings: the expected --skipped_variants_file output")
    ap.add_argument("--actual-skipped", type=Path, help="warnings: the actual --skipped_variants_file output")
    ap.add_argument(
        "--unordered",
        action="store_true",
        help="compare the shared entries without their order (a run whose reference rows Ensembl VEP "
        "emits in hash order: --pick_allele, --pick_allele_gene, --per_gene)",
    )
    ap.add_argument("--report", type=Path, help="where the JSON report is written")
    ap.add_argument("--examples", type=int, default=EXAMPLES_PER_GROUP, help="examples kept per difference group")
    return ap


def run(args: argparse.Namespace) -> dict:
    if args.keep_only is not None and (args.warnings or args.self_parity):
        raise ComparatorError("--keep-only is available only in --format comparison mode")
    if args.warnings:
        if args.expected is None or args.actual is None:
            raise ComparatorError("--warnings needs --expected and --actual")
        return compare_warnings(args.expected, args.actual, args.expected_skipped, args.actual_skipped, parse_masks(args.masks))
    if args.self_parity:
        missing = [name for name in ("default", "tab", "vcf", "json") if getattr(args, name) is None]
        if missing:
            raise ComparatorError("--self-parity needs --" + ", --".join(missing))
        inputs = input_records(args.input) if args.input is not None else None
        parquet_text = parquet_as_tab(args.parquet) if args.parquet is not None else None
        report = SelfParity(
            read_text(args.default),
            read_text(args.tab),
            read_text(args.vcf),
            read_text(args.json),
            parquet_text,
            inputs,
            args.examples,
        ).run()
        report["files"] = {
            name: str(getattr(args, name)) for name in ("default", "tab", "vcf", "json", "parquet", "input")
            if getattr(args, name) is not None
        }
        return report
    if args.format is None or args.expected is None or args.actual is None:
        raise ComparatorError("--format, --expected and --actual are required")
    masks = parse_masks(args.masks)
    documented = Documented.from_manifest(load_manifest(args.manifest), masks)
    inputs = input_records(args.input) if args.input is not None else None
    floors = parse_floors(args.floors)
    field_pairs = parse_field_pairs(args.field_pairs)
    nearest_sites = NearestSites.from_cache(args.nearest_sites) if args.nearest_sites is not None else None
    return compare_files(
        args.format,
        args.expected,
        args.actual,
        inputs,
        documented,
        masks,
        floors,
        field_pairs,
        args.examples,
        nearest_sites,
        unordered=args.unordered,
        keep_only=args.keep_only,
    )


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = run(args)
        write_report(report, args.report)
    except ComparatorError as exc:
        print(f"ERROR: [{COMPONENT}] {exc}", file=sys.stderr)
        return 2
    label = report.get("format") or report["mode"]
    print(f"{label}: {report['verdict']}")
    for reason in report["reasons"]:
        print(f"  {reason}")
    return 0 if report["verdict"] == "green" else 1


if __name__ == "__main__":
    sys.exit(main())
