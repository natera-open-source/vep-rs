"""Write the supplementary columns of fastVEP's tab output as plugin tuples.

A fastVEP run with `--sa-dir` adds one tab column per loaded source after the
built-in columns. This script writes those columns as plugin tuples, the form
in which a comparator's values are scored against Ensembl VEP's plugin and
`--custom` fields: tab-separated, five columns, no header: the row's Location,
Allele and Feature, the Ensembl VEP field name, the decoded value. It reads the
tab output as fastVEP writes it or as fastvep_normalize.py rewrites it: the two
carry the same supplementary columns, and the input is never changed.

A line ends at its line feed, with the carriage return before it where the
file has one. Any other carriage return is part of its cell and is written
`%0D`, in the Location, Allele and Feature cells too.

The column header is the last line before the first row that starts with one
`#`. A cell of `-` holds nothing. Any other cell holds one entry, or several
joined with `,`. An entry is a `|`-joined string whose first part, the
alternate allele or the gene symbol, is dropped. `fastvep_field_map.tsv` names
the Ensembl VEP field of each part after it, one row per field:

    fastvep_column  header name of the supplementary column
    part            N, counted from 1 after the dropped part; N-M for parts N
                    to M joined with `|`; N.word or N.score for one half of a
                    part written `word(score)`
    vep_field       Ensembl VEP field name the value is written under
    source          plugin or `--custom` source the field belongs to
    kind            text or number, the value rule a scorer applies
    note            free text
    list_separator  `-`, or the separator Ensembl VEP writes where the column
                    joins a list with `&`
    reference_set   `-`, or the expected set of the source's committed corpus
                    whose reference declares the field

A value is percent-decoded, except that a tab and every character that
breaks a line stay encoded, so a tuple stays one line. A number keeps the text
fastVEP printed. An empty part writes no tuple. A populated part without a
row in the map, and the populated half of a `word(score)` part without one,
is dropped and counted. Equal values from several entries of one cell are
written once. The input and the tuple file are UTF-8 whatever the locale, and
a byte that is not UTF-8, raw or percent-encoded, is written as it is.
`--fields` keeps the named Ensembl VEP fields alone, for a run whose database
or reference holds fewer fields than the map; `--field-map` reads another map.

One summary line on stderr gives the tuples written, the rows read and the
rows that end before the last mapped column, the tuples of each kept field of
the mapped columns (a field that wrote none reads 0), the mapped columns, the
dropped parts and halves by column, the values `--fields` left out, the parts
a `word(score)` row could not split, the cells with several entries and the
`FV_` columns the map has no row for.

A run ends with one `ERROR: [fastvep_plugin_tuples]` line and exit status 1
when the input cannot be opened, has no column header line, names none of the
map's columns, names one of them more than once or lacks the column of a field
`--fields` names; when `--fields` names no field, as a whole or beside a comma,
or a field the map lacks; when the map cannot be read or is malformed; and
when the tuple file is the input or cannot be opened. The tuple file is opened
only after those checks, so a refused run leaves an earlier tuple file as it
was. A malformed command line is a usage error with exit status 2. A read or
a write that fails after its file is opened ends in Python's own error, with
a partial tuple file once rows are written, and an exit status that is not 0.

Usage:
    python3 fastvep_plugin_tuples.py <fastvep.tab> <tuples.tsv> \\
        [--field-map <map.tsv>] [--fields <name,name>]
"""

import argparse
import re
import sys
from collections import Counter
from pathlib import Path
from typing import NamedTuple
from urllib.parse import quote, unquote_to_bytes

FIELD_MAP = Path(__file__).with_name("fastvep_field_map.tsv")
MAP_HEADER = [
    "fastvep_column", "part", "vep_field", "source", "kind", "note", "list_separator",
    "reference_set",
]
PART = re.compile(r"(\d+)(?:-(\d+)|\.(word|score))?")
WORD_SCORE = re.compile(r"(?P<word>.*)\((?P<score>.*)\)")
# The tab, and each character `str.splitlines` ends a line at.
LINE_BREAKS = re.compile("[\t\n\x0b\x0c\r\x1c\x1d\x1e\x85" + chr(0x2028) + chr(0x2029) + "]")


class TupleError(Exception):
    """A refusal, reported as one ERROR line."""


class Rule(NamedTuple):
    """One row of the field map: the parts it reads and the field they are written under."""

    first: int
    last: int
    half: str
    field: str
    separator: str
    selected: bool


def load_field_map(path, fields=None):
    """Read the field map as {fastVEP column: [Rule, ...]}.

    `fields` is the set of Ensembl VEP field names a run keeps, or None for all of
    them. A malformed row, a field mapped twice and a kept name the map lacks are
    refused.
    """
    try:
        lines = Path(path).read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise TupleError(f"cannot read the field map {path}: {error.strerror}") from error
    except UnicodeDecodeError as error:
        raise TupleError(f"the field map {path} is not UTF-8 text") from error
    if not lines or lines[0].split("\t") != MAP_HEADER:
        raise TupleError(f"{path}:1: the header is not {' '.join(MAP_HEADER)}")
    rules = {}
    named = set()
    for number, line in enumerate(lines[1:], start=2):
        where = f"{path}:{number}"
        cells = line.split("\t")
        if len(cells) != len(MAP_HEADER):
            raise TupleError(f"{where}: {len(cells)} columns, not {len(MAP_HEADER)}")
        column, part, field, _source, kind, _note, separator, reference_set = cells
        spec = PART.fullmatch(part)
        if spec is None:
            raise TupleError(f"{where}: part {part!r} is not N, N-M, N.word or N.score")
        first = int(spec.group(1))
        last = int(spec.group(2) or first)
        half = spec.group(3) or ""
        if first < 1 or (spec.group(2) and last <= first):
            raise TupleError(f"{where}: part {part!r} must count from 1, and a range must rise")
        if kind not in ("text", "number"):
            raise TupleError(f"{where}: kind {kind!r} is not text or number")
        if not (column and field and separator and reference_set):
            raise TupleError(
                f"{where}: an empty fastvep_column, vep_field, list_separator or reference_set"
            )
        if field in named:
            raise TupleError(f"{where}: {field} is mapped twice")
        named.add(field)
        separator = "" if separator == "-" else separator
        if separator and (half or last > first):
            raise TupleError(f"{where}: a list separator needs a part written N, not {part!r}")
        rules.setdefault(column, []).append(
            Rule(first, last, half, field, separator, fields is None or field in fields)
        )
    unknown = sorted(set(fields or ()) - named)
    if unknown:
        raise TupleError(f"--fields names no row of {path}: {', '.join(unknown)}")
    return rules


def decode(text):
    """Percent-decode one value, leaving a tab and every line-breaking character encoded.

    The escapes are decoded as bytes among the value's own bytes, so the two read as
    one UTF-8 text and a byte that is not UTF-8 stays that byte in the tuple file.
    """
    if "%" in text:
        text = unquote_to_bytes(text.encode("utf-8", "surrogateescape")).decode(
            "utf-8", "surrogateescape"
        )
    return LINE_BREAKS.sub(lambda found: quote(found.group()), text)


def line_cells(line):
    """Split one line of the input into its cells.

    The line feed that ends the line, and the carriage return before it, belong to no
    cell. Any other carriage return is returned percent-encoded, as `decode` writes
    one, so a Location, Allele or Feature cell holding one keeps its tuple on one line.
    """
    text = line.removesuffix("\n").removesuffix("\r")
    if "\r" in text:
        text = text.replace("\r", quote("\r"))
    return text.split("\t")


def kept_fields(value):
    """Return the set of field names `--fields` lists, or None for a run without the option.

    An empty value, and an empty name beside a comma, names no field and is refused.
    """
    if value is None:
        return None
    names = value.split(",")
    if "" in names:
        raise TupleError(
            "--fields names no field" + (f" beside a comma: {value}" if value else "")
        )
    return set(names)


def parts_read(rules):
    """Map each part index a column's rules read to the halves read, or None for all of it."""
    read = {}
    for rule in rules:
        for index in range(rule.first, rule.last + 1):
            if not rule.half:
                read[index] = None
            elif read.get(index, ()) is not None:
                read[index] = {*read.get(index, ()), rule.half}
    return read


def decode_cell(cell, rules, read):
    """Decode one populated supplementary cell.

    `read` is the column's `parts_read`. Returns the (field, value) pairs of the kept
    fields, the populated parts and `word(score)` halves no rule reads (`3`, `1.word`),
    the number of populated values `--fields` leaves out, the number of parts a
    `word(score)` rule could not split, and whether the cell holds several entries.
    """
    pairs = []
    dropped = []
    deselected = 0
    unshaped = set()
    entries = cell.split(",")
    for number, entry in enumerate(entries):
        parts = entry.split("|")
        for rule in rules:
            chunk = parts[rule.first : rule.last + 1]
            if not any(chunk):
                continue
            if rule.half:
                shaped = WORD_SCORE.fullmatch(chunk[0])
                if shaped is None:
                    unshaped.add((number, rule.first))
                    continue
                half = shaped.group(rule.half)
                if not half:
                    continue
                value = decode(half)
            elif rule.last > rule.first:
                value = "|".join(decode(part) for part in chunk)
            elif rule.separator:
                value = rule.separator.join(decode(item) for item in chunk[0].split("&"))
            else:
                value = decode(chunk[0])
            if not rule.selected:
                deselected += 1
            elif (rule.field, value) not in pairs:
                pairs.append((rule.field, value))
        for index in range(1, len(parts)):
            if not parts[index]:
                continue
            if index not in read:
                dropped.append(str(index))
            elif read[index] is not None and (shaped := WORD_SCORE.fullmatch(parts[index])):
                dropped.extend(
                    f"{index}.{half}"
                    for half in ("word", "score")
                    if half not in read[index] and shaped.group(half)
                )
    return pairs, dropped, deselected, len(unshaped), len(entries) > 1


def place_columns(header, rules, fields_named):
    """Find the header's mapped columns and its source columns without a map row.

    A mapped column is (index, name, its rules, their `parts_read`). `fields_named`
    says the run named its fields. A missing header, a header naming none of the
    map's columns, a mapped column named more than once and the absent column of a
    named field are refused. Two columns of one name would write two cells of a row
    under one key.
    """
    if header is None:
        raise TupleError("the input has no column header line")
    columns = [
        (index, name, rules[name], parts_read(rules[name]))
        for index, name in enumerate(header)
        if name in rules
    ]
    if not columns:
        raise TupleError(
            "the field map names none of the header's columns: " + ", ".join(header)
        )
    places = {}
    for index, name, _, _ in columns:
        places.setdefault(name, []).append(str(index + 1))
    repeated = [
        f"{name} (columns {', '.join(found[:-1])} and {found[-1]})"
        for name, found in places.items()
        if len(found) > 1
    ]
    if repeated:
        raise TupleError(
            "the header names a mapped column more than once: " + ", ".join(repeated)
        )
    absent = sorted(
        f"{rule.field} ({column})"
        for column, column_rules in rules.items()
        if column not in header
        for rule in column_rules
        if fields_named and rule.selected
    )
    if absent:
        raise TupleError(
            "the header lacks the column of a field --fields names: " + ", ".join(absent)
        )
    unmapped = [name for name in header if name.startswith("FV_") and name not in rules]
    return columns, unmapped


class TupleWriter:
    """Writes each row's supplementary columns as plugin tuples and counts what it leaves out."""

    def __init__(self, handle, columns, unmapped):
        self.handle = handle
        self.columns = columns
        self.unmapped = unmapped
        self.memo = {}
        self.width = max(5, 1 + max(index for index, _, _, _ in columns))
        self.rows = self.short = self.deselected = self.unshaped = self.several = 0
        self.dropped = Counter()
        self.by_field = {
            rule.field: 0 for _, _, rules, _ in columns for rule in rules if rule.selected
        }

    def write_row(self, line):
        self.rows += 1
        cells = line_cells(line)
        if len(cells) < self.width:
            self.short += 1
        if len(cells) < 5:
            return
        key = f"{cells[1]}\t{cells[2]}\t{cells[4]}\t"
        for index, name, rules, read in self.columns:
            if index >= len(cells) or cells[index] in ("-", ""):
                continue
            cell = cells[index]
            # A source's cell repeats on every row of its allele, so the last decoding is kept.
            last_cell, decoding = self.memo.get(index, (None, None))
            if cell != last_cell:
                decoding = decode_cell(cell, rules, read)
                self.memo[index] = (cell, decoding)
            pairs, dropped, deselected, unshaped, several = decoding
            for field, value in pairs:
                self.handle.write(f"{key}{field}\t{value}\n")
                self.by_field[field] += 1
            self.dropped.update((name, part) for part in dropped)
            self.deselected += deselected
            self.unshaped += unshaped
            self.several += several

    def summary(self):
        dropped = ", ".join(
            f"{name}:{part}={count}"
            for (name, part), count in sorted(
                self.dropped.items(),
                key=lambda item: (item[0][0], int(item[0][1].partition(".")[0]), item[0][1]),
            )
        )
        if dropped:
            dropped = f" ({dropped})"
        by_field = ", ".join(f"{field}={count}" for field, count in self.by_field.items())
        return (
            f"fastvep_plugin_tuples: wrote {sum(self.by_field.values())} tuples "
            f"from {self.rows} rows ({self.short} too short); "
            f"by field: {by_field or 'none'}; "
            f"mapped columns: {', '.join(name for _, name, _, _ in self.columns)}; "
            f"dropped parts with no Ensembl VEP field: {sum(self.dropped.values())}{dropped}; "
            f"values outside --fields: {self.deselected}; "
            f"parts not written word(score): {self.unshaped}; "
            f"cells with several entries: {self.several}; "
            f"unmapped columns: {', '.join(self.unmapped) or 'none'}"
        )


def open_file(path, mode, name):
    """Open the input or the tuple file, refusing a path that cannot be opened.

    Both are UTF-8 text whatever the locale; `surrogateescape` carries a byte that is
    not UTF-8 from the input to the tuple file unchanged. The input is read with a
    line ending at its line feed alone: universal newlines end one at
    a carriage return inside a cell.
    """
    try:
        return open(
            path, mode, encoding="utf-8", errors="surrogateescape",
            newline="\n" if mode == "r" else None,
        )
    except OSError as error:
        verb = "write" if mode == "w" else "read"
        raise TupleError(f"cannot {verb} the {name} {path}: {error.strerror}") from error


def write_tuples(input_path, tuple_path, rules, fields_named):
    """Write the input's supplementary columns to the tuple file; return the summary line.

    The tuple file is opened once the header's columns are placed.
    """
    with open_file(input_path, "r", "input") as rows:
        header = None
        first_row = None
        for line in rows:
            if not line.startswith("#"):
                first_row = line
                break
            if not line.startswith("##"):
                header = line_cells(line.lstrip("#"))
        columns, unmapped = place_columns(header, rules, fields_named)
        try:
            is_input = Path(tuple_path).samefile(input_path)
        except OSError:
            is_input = False
        if is_input:
            raise TupleError(f"the tuple file {tuple_path} is the input")
        with open_file(tuple_path, "w", "tuple file") as handle:
            writer = TupleWriter(handle, columns, unmapped)
            if first_row is not None:
                writer.write_row(first_row)
            for line in rows:
                if not line.startswith("#"):
                    writer.write_row(line)
    return writer.summary()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Write the supplementary columns of fastVEP's tab output as plugin tuples."
    )
    parser.add_argument("input", help="fastVEP tab output, as written or as normalized")
    parser.add_argument("tuples", help="tuple file to write")
    parser.add_argument(
        "--field-map", default=FIELD_MAP, help="field map (default: the one beside this script)"
    )
    parser.add_argument("--fields", help="comma-separated Ensembl VEP field names to keep")
    args = parser.parse_args(argv)
    try:
        kept = kept_fields(args.fields)
        rules = load_field_map(args.field_map, kept)
        summary = write_tuples(args.input, args.tuples, rules, kept is not None)
    except TupleError as error:
        print(f"ERROR: [fastvep_plugin_tuples] {error}", file=sys.stderr)
        return 1
    print(summary, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
