#!/usr/bin/env python3
"""Write a MaxEntScan matrix directory in the layout the MaxEntScan plugin reads, with
fabricated scores, and score k-mers against it as the Ensembl VEP plugin does.

The MaxEntScan plugin of Ensembl VEP (`MaxEntScan.pm`) reads a directory of eleven files
(MaxEntScan.pm:104-105, 824-826; opened at :737, :753, :831), one number (or one 7-mer) per
line with no header (:738-743, :754-759, :832-837), and this tool writes each with the line
count the plugin's base-4 hashing implies (:803-820, :851-859):

    me2x5                         16384 lines, the 5' score of each 7-mer, in the order of
    splicemodels/splice5sequences 16384 lines, the 7-mers (a permutation of all of them)
    splicemodels/me2x3acc1..5     16384 lines each, indexed by the hash of a 7-mer
    splicemodels/me2x3acc6, 8     64 lines each, indexed by the hash of a 3-mer
    splicemodels/me2x3acc7, 9     256 lines each, indexed by the hash of a 4-mer

The package's matrices are not redistributed here. This tool writes that layout from
a seed, every value log-uniform in [0.01, 100] with five decimals, so a corpus whose
reference Ensembl VEP wrote over the generated directory is reproduced by the same seed
and both engines read one copy. `splice5sequences` is shuffled by the seed, so a reader
that indexes `me2x5` by the 7-mer's hash instead of by its line in that file scores wrongly.

The plugin keeps each line with its first whitespace character removed (MaxEntScan.pm:739-741,
755-757, 833-835) and multiplies the kept text (:689, :860), so Perl's numeric conversion
decides what a line may hold: blanks before and after the number, and a decimal or an
exponent spelling of it. A second generator, seeded from the same seed, picks each line's
spelling: a value stands as its five-place decimal or with an exponent of either sign, a
line of `me2x5` stands between spaces or tabs, and the directory holds files the plugin
never opens, some named as a table in a directory the plugin does not read it from. The
values of a seed are the same under every spelling; `--plain-lines` writes the five-place
decimals alone and nothing beside the eleven tables.

The scoring functions are the plugin's (`score5` MaxEntScan.pm:665-696 and 764-797,
`score3` :698-726 and 803-892), operation for operation: the consensus factor at
positions 3 and 4 (5') or 18 and 19 (3') over the fixed background and consensus
tables, the product of the table entries in the plugin's order, `log(x) / log(2)`, and
`%.3f`. A score closer than 1e-9 to a `%.3f` rounding boundary could print differently
under a one-unit difference in the last binary place between two engines, so
`--check-kmers` refuses a fixture whose k-mers, or whose paired differences, sit there.

Usage:
    maxentscan_fixture_matrices.py --seed 1 --out <dir> [--force] [--plain-lines]
        writes the directory and prints a JSON summary (seed, line count and sha256 per
        table, size and sha256 per other file); an existing file is overwritten only
        under --force
    maxentscan_fixture_matrices.py --dir <dir> --score <kmer>...
        prints each 9-mer's score5 or 23-mer's score3 as `<kmer>\\t<%.3f>`
    maxentscan_fixture_matrices.py --dir <dir> --check-kmers <file>
        one k-mer per line, or two comma-separated k-mers whose difference is also
        checked; exits 1 naming every k-mer or pair at a rounding boundary
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import sys
from pathlib import Path

BASES = "ACGT"

# `score5_scoreconsensus` and `score3_scoreconsensus` (MaxEntScan.pm:771-792, 871-892).
BACKGROUND = {"A": 0.27, "C": 0.23, "G": 0.23, "T": 0.27}
DONOR_CONSENSUS_1 = {"A": 0.004, "C": 0.0032, "G": 0.9896, "T": 0.0032}
DONOR_CONSENSUS_2 = {"A": 0.0034, "C": 0.0039, "G": 0.0042, "T": 0.9884}
ACCEPTOR_CONSENSUS_1 = {"A": 0.9903, "C": 0.0032, "G": 0.0034, "T": 0.0030}
ACCEPTOR_CONSENSUS_2 = {"A": 0.0027, "C": 0.0037, "G": 0.9905, "T": 0.0030}

# The acceptor tables by file name, the k-mer length each is indexed by and the slice of
# the 21 non-consensus bases it scores (`score3_maxentscore`, MaxEntScan.pm:851-859).
ACCEPTOR_TABLES = (
    ("me2x3acc1", 7, slice(0, 7)),
    ("me2x3acc2", 7, slice(7, 14)),
    ("me2x3acc3", 7, slice(14, 21)),
    ("me2x3acc4", 7, slice(4, 11)),
    ("me2x3acc5", 7, slice(11, 18)),
    ("me2x3acc6", 3, slice(4, 7)),
    ("me2x3acc7", 4, slice(7, 11)),
    ("me2x3acc8", 3, slice(11, 14)),
    ("me2x3acc9", 4, slice(14, 18)),
)

LOW, HIGH = 0.01, 100.0

FIRST_TABLE = "me2x5"
SEQUENCES = "splicemodels/splice5sequences"

# The blanks written before and after a value of `me2x5`: `s/\s//` removes one, and Perl's
# numeric conversion skips the ones left.
BLANKS = (("  ", ""), ("    ", ""), ("", "  "), (" ", " "), ("\t", ""), ("\t", " "), (" ", ""))

# Files the plugin never opens: two named as a table in the directory the plugin does not
# read that table from, one named as a tenth acceptor table, and a note in each directory.
OTHER_FILES = (
    "notes.txt",
    "splice5sequences",
    "splicemodels/me2x5",
    "splicemodels/me2x3acc10",
    "splicemodels/notes.txt",
)
BOUNDARY_MARGIN = 1e-9
_LOG2 = math.log(2)


def hashseq(seq: str) -> int:
    """The base-4 hash of `seq` (`score3_hashseq`, MaxEntScan.pm:803-820): A 0, C 1, G 2, T 3,
    the first base most significant; `hashseq("CAGAAGT")` is 4619."""
    total = 0
    for base in seq.upper():
        total = total * 4 + BASES.index(base)
    return total


def all_kmers(k: int) -> list[str]:
    """Every k-mer over ACGT, in base-4 order (the order a hash of 0, 1, 2, ... names)."""
    kmers = [""]
    for _ in range(k):
        kmers = [prefix + base for prefix in kmers for base in BASES]
    return kmers


def file_lines(seed: int) -> dict[str, list[str]]:
    """The eleven files as lists of lines, by their path under the directory, each number
    as its five-place decimal."""
    rng = random.Random(seed)

    def values(count: int) -> list[str]:
        return [f"{10 ** rng.uniform(-2, 2):.5f}" for _ in range(count)]

    sequences = all_kmers(7)
    rng.shuffle(sequences)
    files = {
        "me2x5": values(len(sequences)),
        "splicemodels/splice5sequences": sequences,
    }
    for name, k, _ in ACCEPTOR_TABLES:
        files[f"splicemodels/{name}"] = values(4**k)
    return files


def spelled_lines(seed: int) -> dict[str, list[str]]:
    """The eleven files as they are written: each number of `file_lines` as its five-place
    decimal or with an exponent, the same value either way, and each line of `me2x5`
    between the blanks of one entry of `BLANKS`. The spellings come from a generator of
    their own, so the values of a seed do not depend on them."""
    rng = random.Random(f"{seed}/line spellings")
    files = file_lines(seed)
    for name, lines in files.items():
        if name == SEQUENCES:
            continue
        spelled = []
        for value in lines:
            text = value if rng.random() < 0.25 else f"{float(value):.6e}"
            if name == FIRST_TABLE:
                before, after = rng.choice(BLANKS)
                text = before + text + after
            spelled.append(text)
        files[name] = spelled
    return files


def other_files(seed: int) -> dict[str, bytes]:
    """The files of `OTHER_FILES`, by their path under the directory: fabricated numbers of
    no table's line count, and a note."""
    rng = random.Random(f"{seed}/other files")
    files = {}
    for name in OTHER_FILES:
        if name.endswith("notes.txt"):
            files[name] = b"A fabricated file. The plugin opens the eleven tables alone.\n"
        else:
            lines = [f"{10 ** rng.uniform(-2, 2):.5f}" for _ in range(rng.randint(3, 9))]
            files[name] = ("\n".join(lines) + "\n").encode("ascii")
    return files


def write_directory(out: Path, seed: int, force: bool, plain_lines: bool = False) -> dict:
    """Writes the directory under `out` and returns the summary: the seed, per table its line
    count and the sha256 of its bytes and, per other file, its size and sha256. Under
    `plain_lines` every number is its five-place decimal and no other file is written. An
    existing file stops the write unless `force`."""
    files = file_lines(seed) if plain_lines else spelled_lines(seed)
    others = {} if plain_lines else other_files(seed)
    if not force:
        present = [name for name in (*files, *others) if (out / name).exists()]
        if present:
            raise FileExistsError(
                f"{out}: {', '.join(present)} exist; pass --force to overwrite"
            )
    summary: dict[str, dict] = {}
    for name, lines in files.items():
        path = out / name
        path.parent.mkdir(parents=True, exist_ok=True)
        data = ("\n".join(lines) + "\n").encode("ascii")
        path.write_bytes(data)
        summary[name] = {"lines": len(lines), "sha256": hashlib.sha256(data).hexdigest()}
    other_summary: dict[str, dict] = {}
    for name, data in others.items():
        (out / name).write_bytes(data)
        other_summary[name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    return {"seed": seed, "files": summary, "other_files": other_summary}


class Tables:
    """The matrices as the plugin holds them after `new`: `me2x5` by the line number the
    7-mer has in `splice5sequences` (`score5_makesequencematrix` and
    `score5_makescorematrix`, MaxEntScan.pm:732-762; a repeated 7-mer keeps its last line,
    as `$matrix{$_} = $n` does) and the nine acceptor tables by line number
    (`score3_makemaxentscores`, :822-843)."""

    def __init__(self, directory: Path):
        donor_scores = _numbers(directory / "me2x5")
        self.sequence_index = {
            seq: n for n, seq in enumerate(_lines(directory / "splicemodels/splice5sequences"))
        }
        self.me2x5 = donor_scores
        self.acceptor = [_numbers(directory / "splicemodels" / name) for name, _, _ in ACCEPTOR_TABLES]

    def score5(self, seq: str) -> float:
        """`score5` (MaxEntScan.pm:665-696) of a 9-mer: the consensus factor of positions 3
        and 4 times the `me2x5` entry of the other seven bases, as log2."""
        seq = seq.upper()
        if len(seq) != 9:
            raise ValueError(f"score5 takes a 9-mer, got {seq!r}")
        consensus = (
            DONOR_CONSENSUS_1[seq[3]]
            * DONOR_CONSENSUS_2[seq[4]]
            / (BACKGROUND[seq[3]] * BACKGROUND[seq[4]])
        )
        rest = seq[0:3] + seq[5:9]
        table_score = self.me2x5[self.sequence_index[rest]]
        return math.log(consensus * table_score) / _LOG2

    def score3(self, seq: str) -> float:
        """`score3` (MaxEntScan.pm:698-726) of a 23-mer: the consensus factor of positions
        18 and 19 times the maximum-entropy ratio of the other 21 bases, as log2."""
        seq = seq.upper()
        if len(seq) != 23:
            raise ValueError(f"score3 takes a 23-mer, got {seq!r}")
        consensus = (
            ACCEPTOR_CONSENSUS_1[seq[18]]
            * ACCEPTOR_CONSENSUS_2[seq[19]]
            / (BACKGROUND[seq[18]] * BACKGROUND[seq[19]])
        )
        rest = seq[0:18] + seq[20:23]
        sc = [
            table[hashseq(rest[part])]
            for table, (_, _, part) in zip(self.acceptor, ACCEPTOR_TABLES)
        ]
        ratio = sc[0] * sc[1] * sc[2] * sc[3] * sc[4] / (sc[5] * sc[6] * sc[7] * sc[8])
        return math.log(consensus * ratio) / _LOG2

    def score(self, kmer: str) -> float:
        """`score5` of a 9-mer, `score3` of a 23-mer."""
        if len(kmer) == 9:
            return self.score5(kmer)
        if len(kmer) == 23:
            return self.score3(kmer)
        raise ValueError(f"a k-mer is scored at length 9 or 23, got {kmer!r}")


def _lines(path: Path) -> list[str]:
    return path.read_text(encoding="ascii").splitlines()


def _numbers(path: Path) -> list[float]:
    return [float(line) for line in _lines(path)]


def format_score(score: float) -> str:
    """The plugin's rendering of a score (`sprintf('%.3f', ...)`, MaxEntScan.pm:245)."""
    return f"{score:.3f}"


def boundary_distance(score: float) -> float:
    """How far `score` sits from the nearest value where `%.3f` changes its last digit."""
    scaled = score * 1000.0
    return abs(scaled - math.floor(scaled) - 0.5) / 1000.0


def at_boundary(score: float) -> bool:
    return boundary_distance(score) < BOUNDARY_MARGIN


def check_kmers(tables: Tables, specs: list[str]) -> list[str]:
    """The k-mers, and `a,b` pairs whose difference `a - b` the plugin also prints, whose
    score sits at a rounding boundary."""
    offenders = []
    for spec in specs:
        parts = spec.split(",")
        scores = [tables.score(part) for part in parts]
        if len(parts) == 2:
            scores.append(scores[0] - scores[1])
        if any(at_boundary(s) for s in scores):
            offenders.append(spec)
    return offenders


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--seed", type=int, help="the seed of the directory --out writes")
    parser.add_argument("--out", type=Path, help="the directory to write")
    parser.add_argument("--force", action="store_true", help="overwrite existing files")
    parser.add_argument(
        "--plain-lines",
        action="store_true",
        help="write each number as its five-place decimal and no file beside the eleven tables",
    )
    parser.add_argument("--dir", type=Path, help="a written directory to score against")
    parser.add_argument("--score", nargs="+", metavar="KMER", help="k-mers to score")
    parser.add_argument("--check-kmers", type=Path, help="k-mers and pairs to check")
    args = parser.parse_args(argv)

    if args.out is not None:
        if args.seed is None:
            parser.error("--out needs --seed")
        try:
            summary = write_directory(args.out, args.seed, args.force, args.plain_lines)
        except FileExistsError as e:
            print(f"ERROR: [maxentscan_fixture_matrices] {e}", file=sys.stderr)
            return 1
        json.dump(summary, sys.stdout, indent=1)
        print()
        return 0

    if args.dir is None:
        parser.error("give --out with --seed, or --dir with --score or --check-kmers")
    tables = Tables(args.dir)
    if args.score:
        for kmer in args.score:
            print(f"{kmer}\t{format_score(tables.score(kmer))}")
        return 0
    if args.check_kmers is not None:
        specs = [line.strip() for line in _lines(args.check_kmers) if line.strip()]
        offenders = check_kmers(tables, specs)
        for spec in offenders:
            print(f"ERROR: [maxentscan_fixture_matrices] rounding boundary: {spec}", file=sys.stderr)
        return 1 if offenders else 0
    parser.error("--dir needs --score or --check-kmers")
    return 2


if __name__ == "__main__":
    sys.exit(main())
