#!/usr/bin/env python3
"""Write a golden corpus's synthetic conservation database for LoFTEE from its own cache.

The `master` branch of LoFTEE opens the value of `conservation_file` as a SQLite database
(`dbi:SQLite:dbname=<value>`, LoF.pm:126) and reads two tables of it, one row a query:

    phylocsf_data  `SELECT * FROM phylocsf_data WHERE transcript = ? AND exon_number = ?`
                   (LoF.pm:570), bound to the transcript's stable id and the exon number
                   of a stop-gained or frameshift row that has one (LoF.pm:276); the
                   columns read are `corresponding_orf_score` and `max_score`
                   (LoF.pm:293-296). No row writes `PHYLOCSF_TOO_SHORT` in `LoF_info`; a
                   row writes `ANN_ORF` and `MAX_ORF` there and, under a negative first
                   score, the flag `PHYLOCSF_UNLIKELY_ORF` when the second score is
                   positive and `PHYLOCSF_WEAK` otherwise.
    gerp_exons     `SELECT * FROM gerp_exons where transcript_id = ? AND exon_num = ?`
                   (gerp_dist.pl:129), bound to the transcript's stable id and the place,
                   counted from 1, of an exon the walk to the stop codon crosses whole
                   (gerp_dist.pl:59-60); the column read is `gerp` (gerp_dist.pl:133),
                   the exon's total, and an exon without a row adds nothing to the
                   walk's total (gerp_dist.pl:64).

The `grch38` branch opens the value the same way (LoF.pm:126 there) and reads one table:

    phylocsf_summary  `SELECT * FROM phylocsf_summary WHERE transcript = ? AND exon = ?`
                   (LoF.pm:556 there), bound to the same two values (LoF.pm:553-554,
                   557), read for the same two columns and written as the same items
                   and flags (LoF.pm:274-284). That branch takes every interval of
                   the walk from its bigWig and has the `gerp_exons` read commented
                   out (gerp_dist.pl:58-59 there).

This tool writes a file of either layout with made-up values and reads no database of the
plugin. Under `--branch master`, the default:

    CREATE TABLE gerp_exons (transcript_id TEXT, exon_num INTEGER, gerp REAL);
    CREATE TABLE phylocsf_data (transcript TEXT, exon_number INTEGER,
                                corresponding_orf_score REAL, max_score REAL);

Under `--branch grch38` the file holds those two tables, row for row as a `master` file
of the same cache and seed holds them, and the table that branch queries:

    CREATE TABLE phylocsf_summary (transcript TEXT, exon INTEGER,
                                   corresponding_orf_score REAL, max_score REAL);

The `grch38` branch reads neither of the first two, and the rows of `phylocsf_summary` are
drawn apart from those of `phylocsf_data`, so a run on that branch's layout that read
another table than `phylocsf_summary` would write other items.

The plugin gives both values of a query without a type (LoF.pm:571, gerp_dist.pl:130;
LoF.pm:557 on `grch38`), and the SQLite driver binds such a value as text. The key
columns therefore carry a declared type: SQLite compares an INTEGER column with the bound
text as numbers, and an exon number in a column declared without a type equals no text.

The rows cover every protein-coding transcript of the corpus's `json_cache` (the plugin
annotates no other biotype, LoF.pm:170; LoF.pm:156 on `grch38`) and each of its exons,
numbered from 1 in the transcript's own order, the order the cache lists them in. A value
is a function of the seed, the table, the transcript's stable id and the exon number
alone, so a row does not move when the cache gains or loses a transcript:

    phylocsf_data, phylocsf_summary
                   three exons in eighteen get no row, and the others fall evenly into a
                   non-negative annotated score, a negative annotated score under a
                   positive highest score, and a negative annotated score under a highest
                   score that is negative or exactly 0. A score has at most three
                   decimals, the highest score is never below the annotated one, and
                   exact zeros stand among the annotated and among the highest scores.
    gerp_exons     one exon in eight gets no row, and the others a total between -300
                   and 900 with two decimals.

Ensembl VEP prints a score as Perl prints a number, with up to fifteen significant
digits, so a value of at most three decimals reads as the same text in both engines.

The file is written under `<corpus>/data/`, or under the directory `--dir` names (a
corpus that passes `--dir_plugins {plugins}` keeps its fixtures under `plugin_data/`), and
its name carries `synthetic`. The bytes of a SQLite file depend on the library that
writes it, so a seed fixes the rows and not the bytes: `read_rows` and `read_summary_rows`
read the tables of a file back, and the tests hold a committed file to the rows its
recorded seed gives. The branch, the seed, the row counts, the file's sha256 and the
SQLite library's version go under `plugin_fixtures.<name>` in the corpus's
provenance.json, every other key of which is kept. An existing file is overwritten only
under --force.

Usage:
    loftee_conservation_database.py --corpus tests/golden/116/<corpus> \\
        --name phylocsf_gerp.synthetic.sql --seed 1 --dir plugin_data [--force]
    loftee_conservation_database.py --corpus tests/golden/116/<corpus> --branch grch38 \\
        --name loftee.synthetic.sql --seed 1 --dir plugin_data [--force]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import sqlite3
import sys
import tempfile
from pathlib import Path

GENERATOR = "scripts/fixtures/loftee_conservation_database.py"
SYNTHETIC = "synthetic"
# The one biotype the plugin annotates (LoF.pm:170).
CODING_BIOTYPE = "protein_coding"
PHYLOCSF_TABLE = "phylocsf_data"
GERP_TABLE = "gerp_exons"
SCHEMA = (
    f"CREATE TABLE {GERP_TABLE} (transcript_id TEXT, exon_num INTEGER, gerp REAL)",
    f"CREATE TABLE {PHYLOCSF_TABLE} (transcript TEXT, exon_number INTEGER, corresponding_orf_score REAL, max_score REAL)",
)
# The branches of the plugin a file is laid out for. `grch38` queries one table more
# (LoF.pm:556 there), which its layout adds to the two of `master`.
MASTER = "master"
GRCH38 = "grch38"
BRANCHES = (MASTER, GRCH38)
SUMMARY_TABLE = "phylocsf_summary"
SUMMARY_SCHEMA = f"CREATE TABLE {SUMMARY_TABLE} (transcript TEXT, exon INTEGER, corresponding_orf_score REAL, max_score REAL)"
# What the plugin writes for a stop-gained or frameshift row by its exon's lookup
# (LoF.pm:290-298): the `LoF_info` item of an exon without a row, the flags of a negative
# annotated score, and no flag under a non-negative one.
NO_ROW = "PHYLOCSF_TOO_SHORT"
NO_FLAG = "no flag"
UNLIKELY_ORF = "PHYLOCSF_UNLIKELY_ORF"
WEAK = "PHYLOCSF_WEAK"
OUTCOMES = (NO_ROW, NO_FLAG, UNLIKELY_ORF, WEAK)
# The draw of an exon's outcome, in eighteenths: three without a row, five for each other.
OUTCOME_DRAWS = 18
OUTCOME_BOUNDS = ((3, NO_ROW), (8, NO_FLAG), (13, UNLIKELY_ORF), (OUTCOME_DRAWS, WEAK))
# The thousandths a row's scores are multiples of: three, two, one or no decimals.
SCORE_STEPS = (1, 10, 100, 1000)
# The largest sizes, in thousandths, of a non-negative annotated score, of a negative one,
# of a positive highest score, and of what a highest score adds to a non-negative
# annotated one.
MAX_ANNOTATED = 400_000
MAX_NEGATIVE = 300_000
MAX_HIGHEST = 400_000
MAX_ABOVE_ANNOTATED = 200_000
# The shares of an exact zero: among the non-negative annotated scores, and among the
# highest scores that are not positive.
ANNOTATED_ZERO_ONE_IN = 8
HIGHEST_ZERO_ONE_IN = 3
GERP_MISSING_ONE_IN = 8
GERP_HUNDREDTHS = (-30_000, 90_000)

PhylocsfRow = tuple[str, int, float, float]
GerpRow = tuple[str, int, float]


def fail(message: str) -> SystemExit:
    return SystemExit(f"ERROR: [loftee_conservation_database] {message}")


def row_rng(seed: int, table: str, transcript: str, exon_number: int) -> random.Random:
    """One generator per row, seeded from the row's identity, so a value depends on
    nothing but the seed and the row."""
    return random.Random(f"{seed}:{table}:{transcript}:{exon_number}")


def phylocsf_scores(seed: int, transcript: str, exon_number: int, table: str = PHYLOCSF_TABLE) -> tuple[float, float] | None:
    """The annotated and the highest score of one exon's row in `table`, `phylocsf_data` or
    `phylocsf_summary`, or None for an exon without a row there: the outcome is drawn
    first, then one step both scores are multiples of, then the scores inside the
    outcome's bounds."""
    rng = row_rng(seed, table, transcript, exon_number)
    draw = rng.randrange(OUTCOME_DRAWS)
    outcome = next(name for bound, name in OUTCOME_BOUNDS if draw < bound)
    if outcome == NO_ROW:
        return None
    step = rng.choice(SCORE_STEPS)
    if outcome == NO_FLAG:
        annotated = 0 if rng.randrange(ANNOTATED_ZERO_ONE_IN) == 0 else step * rng.randint(1, MAX_ANNOTATED // step)
        highest = annotated + step * rng.randint(0, MAX_ABOVE_ANNOTATED // step)
    else:
        annotated = -step * rng.randint(1, MAX_NEGATIVE // step)
        if outcome == UNLIKELY_ORF:
            highest = step * rng.randint(1, MAX_HIGHEST // step)
        elif rng.randrange(HIGHEST_ZERO_ONE_IN) == 0:
            highest = 0
        else:
            highest = -step * rng.randint(1, -annotated // step)
    return annotated / 1000, highest / 1000


def phylocsf_outcome(scores: tuple[float, float] | None) -> str:
    """The entry of `OUTCOMES` the plugin's tests give an exon's lookup (LoF.pm:290-298):
    a row or none, the annotated score below 0 or not, the highest score above 0 or not."""
    if scores is None:
        return NO_ROW
    annotated, highest = scores
    if annotated >= 0:
        return NO_FLAG
    return UNLIKELY_ORF if highest > 0 else WEAK


def gerp_exon_total(seed: int, transcript: str, exon_number: int) -> float | None:
    """The total of one exon's `gerp_exons` row, two decimals within `GERP_HUNDREDTHS`, or
    None for an exon without a row."""
    rng = row_rng(seed, GERP_TABLE, transcript, exon_number)
    if rng.randrange(GERP_MISSING_ONE_IN) == 0:
        return None
    return rng.randint(*GERP_HUNDREDTHS) / 100


def read_coding_exon_counts(cache: Path) -> dict[str, int]:
    """The exon count of every protein-coding transcript of a JSON cache that has an exon,
    by stable id, in the ids' order. A transcript the cache lists in several shards counts
    once; one listed with two exon counts is refused."""
    counts: dict[str, int] = {}
    for shard in sorted(path for path in (cache / "transcripts").rglob("*") if path.is_file()):
        opener = gzip.open if shard.suffix == ".gz" else open
        with opener(shard, "rt", encoding="utf-8") as handle:
            transcripts = json.load(handle)
        for transcript in transcripts:
            if not transcript or transcript.get("biotype") != CODING_BIOTYPE:
                continue
            exons = len(transcript.get("exons") or [])
            if not exons:
                continue
            stable_id = transcript["stable_id"]
            if counts.setdefault(stable_id, exons) != exons:
                raise fail(f"{cache}: {stable_id} is listed with {counts[stable_id]} exons and with {exons}")
    return dict(sorted(counts.items()))


def database_rows(exon_counts: dict[str, int], seed: int) -> tuple[list[PhylocsfRow], list[GerpRow], dict[str, int]]:
    """The rows of `phylocsf_data` and of `gerp_exons` for the transcripts of `exon_counts`,
    each table's in transcript and exon order, with the counts recorded for the file:
    `transcripts`, `exons` and the rows of each table."""
    phylocsf: list[PhylocsfRow] = []
    gerp: list[GerpRow] = []
    for transcript in sorted(exon_counts):
        for exon_number in range(1, exon_counts[transcript] + 1):
            scores = phylocsf_scores(seed, transcript, exon_number)
            if scores is not None:
                phylocsf.append((transcript, exon_number, *scores))
            total = gerp_exon_total(seed, transcript, exon_number)
            if total is not None:
                gerp.append((transcript, exon_number, total))
    counts = {"transcripts": len(exon_counts), "exons": sum(exon_counts.values()), PHYLOCSF_TABLE: len(phylocsf), GERP_TABLE: len(gerp)}
    return phylocsf, gerp, counts


def summary_rows(exon_counts: dict[str, int], seed: int) -> list[PhylocsfRow]:
    """The rows of `phylocsf_summary` for the transcripts of `exon_counts`, in transcript
    and exon order: one per exon that draws a row in that table."""
    return [
        (transcript, exon_number, *scores)
        for transcript in sorted(exon_counts)
        for exon_number in range(1, exon_counts[transcript] + 1)
        if (scores := phylocsf_scores(seed, transcript, exon_number, table=SUMMARY_TABLE)) is not None
    ]


def write_database(out: Path, phylocsf: list[PhylocsfRow], gerp: list[GerpRow], summary: list[PhylocsfRow] | None = None) -> None:
    """Writes the two tables of `SCHEMA` with their rows to `out` and, when `summary` is
    given, the table of `SUMMARY_SCHEMA` with those rows, from a file built in a scratch
    directory beside the target."""
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=out.parent) as scratch_name:
        built = Path(scratch_name) / out.name
        connection = sqlite3.connect(built)
        try:
            for statement in SCHEMA:
                connection.execute(statement)
            connection.executemany(f"INSERT INTO {GERP_TABLE} VALUES (?, ?, ?)", gerp)
            connection.executemany(f"INSERT INTO {PHYLOCSF_TABLE} VALUES (?, ?, ?, ?)", phylocsf)
            if summary is not None:
                connection.execute(SUMMARY_SCHEMA)
                connection.executemany(f"INSERT INTO {SUMMARY_TABLE} VALUES (?, ?, ?, ?)", summary)
            connection.commit()
        finally:
            connection.close()
        built.replace(out)


def read_rows(path: Path) -> tuple[list[PhylocsfRow], list[GerpRow]]:
    """The rows of `phylocsf_data` and of `gerp_exons` in a written file, each table's in
    transcript and exon order. The file is opened read-only."""
    connection = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
    try:
        phylocsf = connection.execute(
            f"SELECT transcript, exon_number, corresponding_orf_score, max_score FROM {PHYLOCSF_TABLE} ORDER BY transcript, exon_number"
        ).fetchall()
        gerp = connection.execute(f"SELECT transcript_id, exon_num, gerp FROM {GERP_TABLE} ORDER BY transcript_id, exon_num").fetchall()
    finally:
        connection.close()
    return phylocsf, gerp


def read_summary_rows(path: Path) -> list[PhylocsfRow]:
    """The rows of `phylocsf_summary` in a file of the `grch38` layout, in transcript and
    exon order. The file is opened read-only; a file without the table is SQLite's error."""
    connection = sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True)
    try:
        return connection.execute(
            f"SELECT transcript, exon, corresponding_orf_score, max_score FROM {SUMMARY_TABLE} ORDER BY transcript, exon"
        ).fetchall()
    finally:
        connection.close()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def update_provenance(corpus: Path, out: Path, branch: str, seed: int, counts: dict[str, int]) -> None:
    """Rewrites `plugin_fixtures.<file name>` in the corpus's provenance.json and keeps every
    other key, other fixtures of the corpus included."""
    path = corpus / "provenance.json"
    record = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    if not isinstance(record, dict):
        raise fail(f"{path}: the record is not a JSON object")
    fixtures = record.setdefault("plugin_fixtures", {})
    if not isinstance(fixtures, dict):
        raise fail(f"{path}: plugin_fixtures is not a JSON object")
    fixtures[out.name] = {
        "generator": GENERATOR,
        "branch": branch,
        "seed": seed,
        "rows": counts,
        "files": {str(out.relative_to(corpus)): sha256(out)},
        "tools": {"sqlite": sqlite3.sqlite_version},
    }
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    tmp.replace(path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", required=True, help="the corpus directory holding json_cache")
    ap.add_argument("--name", required=True, help=f"the file name under <corpus>/<dir>/, which carries `{SYNTHETIC}`")
    ap.add_argument("--dir", "--data-dir", dest="dir", default="data", help="the directory inside the corpus the file is written to (default: data)")
    ap.add_argument(
        "--branch",
        default=MASTER,
        choices=BRANCHES,
        help=f"the branch of the plugin the file is laid out for (default: {MASTER}); {GRCH38} adds the table `{SUMMARY_TABLE}`",
    )
    ap.add_argument("--seed", required=True, type=int, help="the seed every value derives from")
    ap.add_argument("--force", action="store_true", help="overwrite an existing file")
    args = ap.parse_args(argv)

    corpus = Path(args.corpus)
    if "/" in args.name or not args.name:
        raise fail(f"--name is a file name, not a path: {args.name!r}")
    if SYNTHETIC not in args.name:
        raise fail(f"--name carries `{SYNTHETIC}`, since every value of the file is made up: {args.name!r}")
    if "/" in args.dir or args.dir in ("", ".."):
        raise fail(f"--dir is a directory name inside the corpus, not a path: {args.dir!r}")
    cache = corpus / "json_cache"
    if not (cache / "transcripts").is_dir():
        raise fail(f"corpus has no json_cache/transcripts, whose transcripts and exons the rows cover: {corpus}")
    out = corpus / args.dir / args.name
    if out.exists() and not args.force:
        raise fail(f"{out} exists; pass --force to overwrite it")
    exon_counts = read_coding_exon_counts(cache)
    if not exon_counts:
        raise fail(f"{cache}: no protein-coding transcript with an exon to write a row for")
    phylocsf, gerp, counts = database_rows(exon_counts, args.seed)
    summary = summary_rows(exon_counts, args.seed) if args.branch == GRCH38 else None
    summary_text = ""
    if summary is not None:
        counts[SUMMARY_TABLE] = len(summary)
        summary_text = f"{SUMMARY_TABLE} {len(summary)} rows, "
    write_database(out, phylocsf, gerp, summary)
    update_provenance(corpus, out, args.branch, args.seed, counts)
    print(
        f"wrote {out.relative_to(corpus)}: {summary_text}{PHYLOCSF_TABLE} {len(phylocsf)} rows, {GERP_TABLE} {len(gerp)} rows, "
        f"over {counts['exons']} exons of {counts['transcripts']} transcripts, seed {args.seed}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
