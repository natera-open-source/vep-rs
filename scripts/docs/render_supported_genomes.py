#!/usr/bin/env python3
"""Write the supported-genomes page, `docs/supported-genomes.md`, from the committed genome table.

Inputs, all under the repository root:
- `tests/genomes/genomes.tsv` and `tests/genomes/bacteria.tsv.gz`: one row a ready-made Ensembl VEP
  cache, in the columns `docs/genome-checks.md` gives.
- `tests/genomes/sets.tsv`: each set's release and how many caches it holds.
- `tests/genomes/outcomes.tsv`: the meaning of every outcome but `exact`.
- `tests/golden/116/*/manifest.json`, with `json_cache/info.json` where the corpus has one: which
  genomes have a committed corpus. A corpus belongs to the species its manifest names, else the one
  its cache names, else `homo_sapiens`, the default of `--species`; its transcript set is `refseq` or
  `merged` where the manifest's flags name `--refseq` or `--merged` or the cache states that source,
  else `ensembl`.
- `docs/cache-species.md`: its level-two headings. A corpus cell links the section whose heading,
  in lower case with underscores for blanks, is the genome's name; a species without such a section
  links the heading `SPECIES_TABLE_HEADING`. A human pangenome assembly, whose species directory is
  `homo_sapiens_` and the assembly's accession, is human and no other species: its corpus links the
  heading `PANGENOME_HEADING` where the page has it. A corpus of `homo_sapiens` links nothing.

The page. A set whose name ends in `_collection` is summarised by collection, one line a set; every
other set gets one row a cache, in the table's own order. Every count on the page is counted from
the rows, and a count of four digits or more carries thousands separators. The line of four counts
gives `exact` of compared for the species caches, the bacterial and archaeal genomes, the fungal and
protist collection genomes and the human pangenome assemblies; a row whose outcome opens
`excluded:` was never compared and counts under neither figure. The sentences on what was compared
set the excluded caches apart only where a row's outcome opens `excluded:`. The sentence on
`--use_given_ref` is written only where a row of `homo_sapiens` under a transcript set of
`EDITED_FROM_ALIGNMENT` was compared.

Refused, each as one `ERROR: [render_supported_genomes] ...` line: a table file that is missing,
does not decode or whose header is not its columns, a row of another width, a set or an outcome its
side table lists twice, a set `sets.tsv` does not list or lists with another count, a row whose
cache version is not its set's, an outcome that is neither `exact` nor a row of `outcomes.tsv`, a
species set without a title in `SPECIES_SETS`, a corpus manifest or cache `info.json` that is no
JSON object, and a page that cannot be written. The same tree gives the same bytes. Stdlib only.

usage: render_supported_genomes.py [--repo DIR] [--check | --counts]

--check writes nothing and exits 1, naming the first line that differs, when the committed page is
not, byte for byte, what the tables give. --counts writes nothing and prints the page's counts,
tab-separated.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import statistics
import sys
import zlib
from collections import Counter
from pathlib import Path

COMPONENT = "render_supported_genomes"
PAGE = "docs/supported-genomes.md"
SPECIES_PAGE = "docs/cache-species.md"
SPECIES_TABLE_HEADING = "Species with a corpus"
PANGENOME_HEADING = "A human pangenome assembly"
TABLE_DIR = "tests/genomes"
CORPORA = "tests/golden/116"
COLUMNS = ["genome", "assembly", "transcripts", "set", "cache_version", "seed", "records", "outcome", "sha256"]
EXACT = "exact"
HUMAN = "homo_sapiens"
# The transcript sets of human's caches whose transcripts are edited from an alignment file. On a
# cache that names such a file `scripts/genomes/check_genome.sh` gives both engines
# `--use_given_ref`: without the option each refuses a run that names no FASTA file.
EDITED_FROM_ALIGNMENT = ("refseq", "merged")
# set -> (title on the page, the name the set's release goes by)
SPECIES_SETS = {
    "vertebrates": ("Ensembl", "Ensembl"),
    "metazoa": ("Ensembl Metazoa", "Ensembl Genomes"),
    "plants": ("Ensembl Plants", "Ensembl Genomes"),
    "fungi": ("Ensembl Fungi", "Ensembl Genomes"),
    "protists": ("Ensembl Protists", "Ensembl Genomes"),
    "pangenome": ("Human pangenome assemblies", "the rapid-release directory"),
}
# prefix of a collection set -> (title, what one of its genomes is called, the file its rows stand in)
COLLECTIONS = {
    "bacteria_": ("Bacterial and archaeal collections", "bacterial and archaeal genomes", "bacteria.tsv.gz"),
    "fungi_": ("Fungal collections", "fungal collection genomes", "genomes.tsv"),
    "protists_": ("Protist collections", "protist collection genomes", "genomes.tsv"),
}
# What a table file that does not decode raises: gzip data cut short or damaged, bytes that are not
# UTF-8, a field over the csv module's limit.
UNREADABLE = (EOFError, UnicodeDecodeError, csv.Error, gzip.BadGzipFile, zlib.error)


class Refusal(Exception):
    """A table, a record or a page this script does not render."""


def number(value: int) -> str:
    return f"{value:,}"


def read_rows(path: Path) -> list[dict[str, str]]:
    """The rows of one table file, plain or gzip, held to the nine columns."""
    if not path.is_file():
        raise Refusal(f"{path}: no such file")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            header = next(reader, None)
            if header != COLUMNS:
                raise Refusal(f"{path}: the header is not the columns {', '.join(COLUMNS)}")
            rows = []
            for line, cells in enumerate(reader, start=2):
                if len(cells) != len(COLUMNS):
                    raise Refusal(f"{path}:{line}: {len(cells)} cells where the header has {len(COLUMNS)}")
                rows.append(dict(zip(COLUMNS, cells)))
    except UNREADABLE as error:
        raise Refusal(f"{path}: {error}") from None
    return rows


def read_pairs(path: Path, columns: list[str]) -> list[dict[str, str]]:
    """The rows of a side table (`sets.tsv`, `outcomes.tsv`) under its own header, each first cell once."""
    if not path.is_file():
        raise Refusal(f"{path}: no such file")
    try:
        with open(path, encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
            if next(reader, None) != columns:
                raise Refusal(f"{path}: the header is not the columns {', '.join(columns)}")
            rows, named = [], set()
            for line, cells in enumerate(reader, start=2):
                if len(cells) != len(columns):
                    raise Refusal(f"{path}:{line}: {len(cells)} cells where the header has {len(columns)}")
                if cells[0] in named:
                    raise Refusal(f"{path}:{line}: a second row of {columns[0]} {cells[0]}")
                named.add(cells[0])
                rows.append(dict(zip(columns, cells)))
    except UNREADABLE as error:
        raise Refusal(f"{path}: {error}") from None
    return rows


def read_object(path: Path) -> dict:
    """The JSON object one file holds."""
    try:
        found = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as error:
        raise Refusal(f"{path}: {error}") from None
    if not isinstance(found, dict):
        raise Refusal(f"{path}: no JSON object")
    return found


def collection_prefix(name: str) -> str | None:
    if not name.endswith("_collection"):
        return None
    for prefix in COLLECTIONS:
        if name.startswith(prefix):
            return prefix
    raise Refusal(f"set {name}: a collection of no kind this script names")


def is_human(genome: str) -> bool:
    """Whether a species directory is human's own or that of a human pangenome assembly."""
    return genome == HUMAN or genome.startswith(HUMAN + "_")


def read_corpora(repo: Path) -> dict[tuple[str, str, str], list[str]]:
    """(genome, assembly, transcripts) -> the corpus directories that annotate it."""
    found: dict[tuple[str, str, str], list[str]] = {}
    root = repo / CORPORA
    for manifest_path in sorted(root.glob("*/manifest.json")):
        manifest = read_object(manifest_path)
        info_path = manifest_path.parent / "json_cache" / "info.json"
        info = read_object(info_path) if info_path.is_file() else {}
        species = manifest.get("species") or info.get("species") or HUMAN
        assembly = manifest.get("assembly") or info.get("assembly")
        if not isinstance(assembly, str) or not assembly:
            raise Refusal(f"{manifest_path}: names no assembly")
        if not isinstance(species, str):
            raise Refusal(f"{manifest_path}: names no species")
        flags = manifest.get("flags") or []
        if not isinstance(flags, list):
            raise Refusal(f"{manifest_path}: `flags` is no list")
        if "--refseq" in flags or info.get("source_type") == "refseq":
            transcripts = "refseq"
        elif "--merged" in flags or info.get("source_type") == "merged":
            transcripts = "merged"
        else:
            transcripts = "ensembl"
        found.setdefault((species, assembly, transcripts), []).append(manifest_path.parent.name)
    return found


def corpus_records(repo: Path, names: list[str]) -> list[int]:
    counts = []
    for name in names:
        records = read_object(repo / CORPORA / name / "manifest.json").get("records")
        if not isinstance(records, list):
            raise Refusal(f"{CORPORA}/{name}/manifest.json: `records` is no list")
        counts.append(len(records))
    return counts


def heading_anchor(heading: str) -> str:
    """The id mdBook gives a heading: lower case, blanks to hyphens, other punctuation dropped."""
    kept = [c for c in heading.strip().lower() if c.isalnum() or c in " -_"]
    return "".join(kept).replace(" ", "-")


def species_sections(repo: Path) -> dict[str, str]:
    """genome name -> anchor, for each level-two heading of the species page that is a Latin name."""
    path = repo / SPECIES_PAGE
    if not path.is_file():
        raise Refusal(f"{path}: no such file")
    anchors = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            heading = line[3:].strip()
            anchors[heading.lower().replace(" ", "_")] = heading_anchor(heading)
    return anchors


def outcome_cells(rows: list[dict[str, str]]) -> tuple[str, str]:
    """The compared outcomes by value, `exact` first, and the excluded ones, each with its count."""
    tally = Counter(row["outcome"] for row in rows)
    compared = [name for name in sorted(tally) if not name.startswith("excluded:")]
    compared.sort(key=lambda name: name != EXACT)
    excluded = [name for name in sorted(tally) if name.startswith("excluded:")]

    def cell(names: list[str]) -> str:
        parts = [f"{name if name == EXACT else f'`{name}`'} {number(tally[name])}" for name in names]
        return "; ".join(parts) if parts else "0"

    return cell(compared), cell(excluded)


def compared_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if not row["outcome"].startswith("excluded:")]


def exact_of_compared(rows: list[dict[str, str]]) -> tuple[int, int]:
    compared = compared_rows(rows)
    return sum(row["outcome"] == EXACT for row in compared), len(compared)


def median_text(values: list[int]) -> str:
    median = statistics.median(values)
    return number(int(median)) if median == int(median) else f"{median:,.1f}"


class Tables:
    """The table files of one tree, grouped the way the page prints them."""

    def __init__(self, repo: Path) -> None:
        table_dir = repo / TABLE_DIR
        self.rows = read_rows(table_dir / "genomes.tsv") + read_rows(table_dir / "bacteria.tsv.gz")
        self.sets = {row["set"]: row for row in read_pairs(table_dir / "sets.tsv", ["set", "release", "cache_version", "caches"])}
        self.meanings = {row["outcome"]: row["meaning"] for row in read_pairs(table_dir / "outcomes.tsv", ["outcome", "meaning"])}
        self.by_set: dict[str, list[dict[str, str]]] = {}
        for row in self.rows:
            self.by_set.setdefault(row["set"], []).append(row)
        for name, rows in self.by_set.items():
            if name not in self.sets:
                raise Refusal(f"set {name}: {TABLE_DIR}/sets.tsv does not list it")
            if int(self.sets[name]["caches"]) != len(rows):
                raise Refusal(f"set {name}: {len(rows)} rows where {TABLE_DIR}/sets.tsv counts {self.sets[name]['caches']}")
            if collection_prefix(name) is None and name not in SPECIES_SETS:
                raise Refusal(f"set {name}: no title in SPECIES_SETS")
        for name, listed in self.sets.items():
            if name not in self.by_set and listed["caches"] != "0":
                raise Refusal(f"set {name}: no row where {TABLE_DIR}/sets.tsv counts {listed['caches']}")
        for row in self.rows:
            version = self.sets[row["set"]]["cache_version"]
            if row["cache_version"] != version:
                raise Refusal(
                    f"{row['genome']} {row['assembly']}: cache version {row['cache_version']} where "
                    f"{TABLE_DIR}/sets.tsv gives set {row['set']} {version}"
                )
            if row["outcome"] != EXACT and row["outcome"] not in self.meanings:
                raise Refusal(f"outcome {row['outcome']}: {TABLE_DIR}/outcomes.tsv does not declare it")
        self.species_sets = [name for name in SPECIES_SETS if name in self.by_set]
        self.collections = {
            prefix: sorted(
                (name for name in self.by_set if collection_prefix(name) == prefix),
                key=lambda name: [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", name)],
            )
            for prefix in COLLECTIONS
        }

    def species_rows(self, pangenome: bool) -> list[dict[str, str]]:
        names = [name for name in self.species_sets if (name == "pangenome") == pangenome]
        return [row for name in names for row in self.by_set[name]]

    def collection_rows(self, prefix: str) -> list[dict[str, str]]:
        return [row for name in self.collections[prefix] for row in self.by_set[name]]

    def set_title(self, name: str) -> str:
        title, _ = SPECIES_SETS[name]
        return title if name == "pangenome" else f"{title} {self.sets[name]['release']}"

    def groups(self) -> list[tuple[str, list[dict[str, str]]]]:
        """Each line of the summary tables: its label and its rows."""
        lines = [(self.set_title(name), self.by_set[name]) for name in self.species_sets if name != "pangenome"]
        for prefix, (title, _, _) in COLLECTIONS.items():
            if self.collections[prefix]:
                lines.append((f"{title} ({number(len(self.collections[prefix]))})", self.collection_rows(prefix)))
        if "pangenome" in self.by_set:
            lines.append((self.set_title("pangenome"), self.by_set["pangenome"]))
        return lines

    def four_counts(self) -> list[tuple[str, int, int]]:
        """(label, exact, compared) for the four sets the release figure's genomes panel names."""
        both = self.collection_rows("fungi_") + self.collection_rows("protists_")
        return [
            ("species caches", *exact_of_compared(self.species_rows(pangenome=False))),
            ("bacterial and archaeal genomes", *exact_of_compared(self.collection_rows("bacteria_"))),
            ("fungal and protist collection genomes", *exact_of_compared(both)),
            ("human pangenome assemblies", *exact_of_compared(self.species_rows(pangenome=True))),
        ]


def render(repo: Path) -> str:
    tables = Tables(repo)
    corpora = read_corpora(repo)
    sections = species_sections(repo)
    table_anchor = heading_anchor(SPECIES_TABLE_HEADING)
    species_page = Path(SPECIES_PAGE).name
    pangenome_anchor = heading_anchor(PANGENOME_HEADING)
    other_species = sorted({genome for genome, _, _ in corpora if not is_human(genome)})
    other_corpora = [name for (genome, _, _), names in corpora.items() if not is_human(genome) for name in names]
    sizes = corpus_records(repo, other_corpora)
    some_excluded = len(compared_rows(tables.rows)) < len(tables.rows)
    given_reference = any(
        row["genome"] == HUMAN and row["transcripts"] in EDITED_FROM_ALIGNMENT for row in compared_rows(tables.rows)
    )

    out = ["# Supported genomes", ""]
    out += [
        (
            "Every genome with a ready-made Ensembl VEP cache is listed here with the outcome of its check. "
            f"Both engines annotated a small generated input for each{' compared' if some_excluded else ''}."
        ),
        "",
        (
            "Human is measured on whole genomes ([Benchmarks](benchmarks/current-release.md)). "
            f"The {number(len(other_species))} other species with a committed corpus are also compared under several "
            f"flag sets, on {number(min(sizes))} to {number(max(sizes))} records a corpus "
            f"([Other species]({species_page})). "
            f"Every cache on this page{' that was not excluded' if some_excluded else ''} was compared once "
            "on the generated records its row counts. "
            "Both engines ran under default flags. "
            + (
                "On the human RefSeq and merged caches, whose transcripts are edited from an alignment file, "
                "both were also given `--use_given_ref` ([Genome checks](genome-checks.md)). "
                if given_reference
                else ""
            )
            + "No genome other than human is timed. "
            "A genome built from a gene-model file is not measured."
        ),
        "",
        (
            "The rows are those of [`tests/genomes/genomes.tsv`](../tests/genomes/genomes.tsv) and "
            "[`tests/genomes/bacteria.tsv.gz`](../tests/genomes/bacteria.tsv.gz). "
            "[Genome checks](genome-checks.md) describes both files, the outcome values and the job that writes a row."
        ),
        "",
    ]
    counts = "; ".join(f"{label} {number(exact)} of {number(compared)}" for label, exact, compared in tables.four_counts())
    out += [f"Exact of compared: {counts}.", ""]

    out += ["## Outcomes", "", "| Outcome | Caches | Meaning |", "| --- | --- | --- |"]
    tally = Counter(row["outcome"] for row in tables.rows)
    out.append(
        f"| `{EXACT}` | {number(tally[EXACT])} | On the genome's records vep-rs matches Ensembl VEP exactly "
        "except for the documented divergences. |"
    )
    for name in sorted(name for name in tally if name != EXACT):
        out.append(f"| `{name}` | {number(tally[name])} | {tables.meanings[name]} |")
    out.append("")

    out += [
        "## Sets",
        "",
        (
            "A genome is a species under one assembly. It can ship as up to three caches: Ensembl, RefSeq and merged "
            "transcripts." + (" A cache that was excluded was never compared." if some_excluded else "")
        ),
        "",
        "| Set | Genomes (caches) | Compared, by outcome | Excluded |",
        "| --- | --- | --- | --- |",
    ]
    for label, rows in tables.groups():
        genomes = len({(row["set"], row["genome"], row["assembly"]) for row in rows})
        compared, excluded = outcome_cells(rows)
        out.append(f"| {label} | {number(genomes)} ({number(len(rows))}) | {compared} | {excluded} |")
    out.append("")

    out += [
        "## Records compared",
        "",
        "How many generated records a compared cache was annotated on, by set.",
        "",
        "| Set | Least | Median | Greatest |",
        "| --- | --- | --- | --- |",
    ]
    few = []
    for label, rows in tables.groups():
        records = [int(row["records"]) for row in compared_rows(rows)]
        if not records:
            continue
        out.append(f"| {label} | {number(min(records))} | {median_text(records)} | {number(max(records))} |")
        under_ten = sum(value < 10 for value in records)
        if under_ten:
            few.append(f"{label}: {number(under_ten)} of {number(len(records))} compared caches had fewer than ten records.")
    out.append("")
    if few:
        out += [*few, ""]

    for name in tables.species_sets:
        rows = tables.by_set[name]
        _, release_name = SPECIES_SETS[name]
        version = tables.sets[name]["cache_version"]
        out += [f"## {tables.set_title(name)}", ""]
        if name == "pangenome":
            out.append(
                f"The {number(len(rows))} assemblies of {release_name}, release `{tables.sets[name]['release']}`. "
                f"Each cache takes `--cache_version {version}`."
            )
        else:
            out.append(
                f"The {number(len(rows))} caches of {release_name} release {tables.sets[name]['release']}. "
                f"Each takes `--cache_version {version}`."
            )
        out += ["", "| Genome | Assembly and transcripts | Outcome | Corpus |", "| --- | --- | --- | --- |"]
        for row in rows:
            key = (row["genome"], row["assembly"], row["transcripts"])
            if key not in corpora:
                corpus = "-"
            elif row["genome"] == HUMAN:
                corpus = "yes"
            elif is_human(row["genome"]):
                linked = pangenome_anchor in sections.values()
                corpus = f"[yes]({species_page}#{pangenome_anchor})" if linked else "yes"
            else:
                corpus = f"[yes]({species_page}#{sections.get(row['genome'], table_anchor)})"
            out.append(f"| `{row['genome']}` | {row['assembly']}, `{row['transcripts']}` | `{row['outcome']}` | {corpus} |")
        out.append("")

    for prefix, (title, noun, file_name) in COLLECTIONS.items():
        names = tables.collections[prefix]
        if not names:
            continue
        rows = tables.collection_rows(prefix)
        out += [
            f"## {title}",
            "",
            (
                f"The {number(len(rows))} {noun} of {number(len(names))} collections are summarised one line a "
                f"collection. Each genome's own row is in [`{TABLE_DIR}/{file_name}`](../{TABLE_DIR}/{file_name})."
            ),
            "",
            "| Collection | Genomes | Compared, by outcome | Excluded |",
            "| --- | --- | --- | --- |",
        ]
        for name in names:
            compared, excluded = outcome_cells(tables.by_set[name])
            out.append(f"| `{name}` | {number(len(tables.by_set[name]))} | {compared} | {excluded} |")
        out.append("")

    out += [
        "## Where next",
        "",
        f"- [Other species]({species_page}): a worked example for each set of caches, and the species with a corpus.",
        "- [Genome checks](genome-checks.md): the table files, the outcome values and how one genome is checked again.",
        "- [Cache conversion](cache-conversion.md): converting a ready-made Ensembl VEP cache.",
        "- [Intended divergences](intended-divergences.md): where the two engines differ, and why.",
    ]
    return "\n".join(out) + "\n"


def counts_report(repo: Path) -> str:
    tables = Tables(repo)
    corpora = read_corpora(repo)
    lines = [f"{label}\t{exact}\t{compared}" for label, exact, compared in tables.four_counts()]
    lines.append(f"species with a corpus, human apart\t{len({genome for genome, _, _ in corpora if not is_human(genome)})}")
    return "\n".join(lines) + "\n"


def first_difference(committed: str, rendered: str) -> str:
    old, new = committed.splitlines(keepends=True), rendered.splitlines(keepends=True)
    for index, (left, right) in enumerate(zip(old, new), start=1):
        if left != right:
            page_line, table_line = left.splitlines()[0], right.splitlines()[0]
            if page_line == table_line:
                return f"line {index}: the page's line end is not the one the tables give"
            return f"line {index}: the page has `{page_line}`, the tables give `{table_line}`"
    return f"line {min(len(old), len(new)) + 1}: the page has {len(old)} lines, the tables give {len(new)}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent.parent)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--counts", action="store_true")
    args = parser.parse_args(argv)
    page = args.repo / PAGE
    try:
        if args.counts:
            sys.stdout.write(counts_report(args.repo))
            return 0
        rendered = render(args.repo)
        if not args.check:
            page.write_bytes(rendered.encode("utf-8"))
            return 0
        committed = page.read_bytes() if page.is_file() else b""
    except (Refusal, OSError, ValueError) as error:
        print(f"ERROR: [{COMPONENT}] {error}", file=sys.stderr)
        return 1
    if committed == rendered.encode("utf-8"):
        return 0
    difference = first_difference(committed.decode("utf-8", errors="backslashreplace"), rendered)
    print(f"ERROR: [{COMPONENT}] {PAGE} is not what the tables give: {difference}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
