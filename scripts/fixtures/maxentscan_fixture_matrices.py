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

`--whole-numbers` writes every number as a whole number below 2**31, drawn from the seed by
a generator of its own, and nothing beside the eleven tables. Perl multiplies two whole
numbers as integers while the product fits 64 bits (`perl_multiply`), so the product of
MaxEntScan.pm:860 over such a table is exact where the same product in doubles has
rounded once per factor. The two agree to the last place or the one before it, which
`%.3f` shows only where a score sits at a rounding boundary: each 23-mer given to
`--boundary-kmers` gets nine entries, one of them of sixteen digits, found by a seeded
search so that its score prints one third decimal under Perl's arithmetic and the next
one under doubles.

The scoring functions are the plugin's (`score5` MaxEntScan.pm:665-696 and 764-797,
`score3` :698-726 and 803-892), operation for operation: the consensus factor at
positions 3 and 4 (5') or 18 and 19 (3') over the fixed background and consensus
tables, the product of the table entries in the plugin's order, multiplied and divided
as Perl multiplies and divides the text of a line, `log(x) / log(2)`, and `%.3f`. A
score closer than 1e-9 to a `%.3f` rounding boundary could print differently under a
one-unit difference in the last binary place between two engines, so `--check-kmers`
refuses a fixture whose k-mers, or whose paired differences, sit there.

Usage:
    maxentscan_fixture_matrices.py --seed 1 --out <dir> [--force] [--plain-lines]
        writes the directory and prints a JSON summary (seed, line count and sha256 per
        table, size and sha256 per other file); an existing file is overwritten only
        under --force
    maxentscan_fixture_matrices.py --seed 1 --out <dir> --whole-numbers
            [--boundary-kmers <23-mer>...]
        writes the directory with whole numbers, the named 23-mers at a rounding
        boundary
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
import re
import sys
from functools import reduce
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

# The whole numbers `--whole-numbers` draws lie in [1, 2**WHOLE_NUMBER_BITS): the product
# of two can pass 2**53, the last power of two below which every whole number is a double.
WHOLE_NUMBER_BITS = 31
_DOUBLE_INTEGERS = 2**53
_IV_MIN, _UV_MAX = -(2**63), 2**64 - 1

# A boundary k-mer's score lies within this many thousandths of zero, where a double of
# the score is fine enough to tell two neighbouring doubles of the product apart, and
# each of its two scores stands at least this far from the boundary between them, more
# than the last place of a logarithm moves either.
BOUNDARY_SPAN = 20
BOUNDARY_CLEARANCE = 3e-17


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


class PointNumber(float):
    """The number of a line written with a decimal point and no exponent, which Perl
    multiplies as a double whatever its value."""


def perl_number(text: str) -> int | float:
    """The number Perl's arithmetic makes of the text of a table line, for the spellings
    this tool writes: a whole number that fits 64 bits is an integer, and so is an
    exponent spelling whose double is a whole number that fits them; a spelling with a
    decimal point and no exponent is a `PointNumber`; any other is a double."""
    text = text.strip(" \t\n\r\f\v")
    if re.fullmatch(r"[+-]?[0-9]+", text):
        value = int(text)
        return value if _IV_MIN <= value <= _UV_MAX else float(text)
    value = float(text)
    if "e" not in text.lower():
        return PointNumber(value)
    if math.isfinite(value) and value == int(value) and _IV_MIN <= int(value) <= _UV_MAX:
        return int(value)
    return value


def _whole(number: float) -> int | None:
    """`number` as the integer Perl's multiply and divide read it as, or `None` where
    they read a double: an integer, or a product or quotient in doubles that is a whole
    number below 2**53."""
    if isinstance(number, int):
        return number
    if isinstance(number, PointNumber) or not math.isfinite(number):
        return None
    if abs(number) < _DOUBLE_INTEGERS and number == int(number):
        return int(number)
    return None


def perl_multiply(left: float, right: float) -> int | float:
    """`left * right` as Perl computes it: the exact product of two integers where it fits
    64 bits, signed or unsigned, else the product of the two doubles."""
    whole_right = _whole(right)
    whole_left = _whole(left) if whole_right is not None else None
    if whole_left is not None and whole_right is not None:
        product = whole_left * whole_right
        if _IV_MIN <= product <= _UV_MAX:
            return product
    return float(left) * float(right)


def perl_divide(left: float, right: float) -> int | float:
    """`left / right` as Perl computes it: the exact quotient of two integers where the
    dividend is past 2**53 and the divisor divides it, else the quotient of the two
    doubles; a zero divisor raises `ZeroDivisionError`, where Perl dies."""
    whole_right = _whole(right)
    whole_left = _whole(left) if whole_right is not None else None
    if whole_left is not None and whole_right is not None:
        if whole_right == 0:
            raise ZeroDivisionError("Illegal division by zero")
        if (
            abs(whole_left) >= abs(whole_right)
            and abs(whole_left) > _DOUBLE_INTEGERS
            and whole_left % whole_right == 0
        ):
            quotient = whole_left // whole_right
            return quotient if quotient >= _IV_MIN else float(quotient)
    return float(left) / float(right)


def consensus_factor(seq: str) -> float:
    """`score3_scoreconsensus` (MaxEntScan.pm:871-892) of a 23-mer: the factor of
    positions 18 and 19."""
    return (
        ACCEPTOR_CONSENSUS_1[seq[18]]
        * ACCEPTOR_CONSENSUS_2[seq[19]]
        / (BACKGROUND[seq[18]] * BACKGROUND[seq[19]])
    )


def acceptor_lines(seq: str) -> list[tuple[str, int]]:
    """The nine lines a 23-mer's `score3` reads, each as its file and its line index
    (`score3_maxentscore`, MaxEntScan.pm:851-859)."""
    rest = seq[0:18] + seq[20:23]
    return [(f"splicemodels/{name}", hashseq(rest[part])) for name, _, part in ACCEPTOR_TABLES]


def ratio_as_perl(sc: list[int | float]) -> int | float:
    """The quotient of MaxEntScan.pm:860 over the nine entries, as Perl computes it."""
    return perl_divide(reduce(perl_multiply, sc[:5]), reduce(perl_multiply, sc[5:]))


def ratio_in_doubles(sc: list[int | float]) -> float:
    """The quotient of MaxEntScan.pm:860 over the nine entries, every entry a double."""
    d = [float(entry) for entry in sc]
    return d[0] * d[1] * d[2] * d[3] * d[4] / (d[5] * d[6] * d[7] * d[8])


def _log2_score(consensus: float, ratio: float) -> float:
    """`log2($a * $c)` of MaxEntScan.pm:719 and 794-797."""
    return math.log(float(perl_multiply(consensus, ratio))) / _LOG2


def whole_number_lines(seed: int) -> dict[str, list[str]]:
    """The eleven files with every number a whole number in [1, 2**WHOLE_NUMBER_BITS),
    log-uniform, from a generator of its own; the 7-mers stand in the order `file_lines`
    gives the seed."""
    rng = random.Random(f"{seed}/whole numbers")
    files = file_lines(seed)
    for name, lines in files.items():
        if name != SEQUENCES:
            files[name] = [str(int(2 ** rng.uniform(0, WHOLE_NUMBER_BITS))) for _ in lines]
    return files


def boundary_entries(seed: int, kmer: str, fixed: dict[int, int]) -> list[int]:
    """Nine whole numbers for the lines `kmer` reads, the entries of `fixed` kept at their
    places, whose score prints one third decimal with the quotient of MaxEntScan.pm:860
    computed as Perl computes it and the next one with it computed in doubles. The first
    entry has sixteen digits, so that the score moves by less than its last place from one
    candidate to the next; the four after it are small, so that Perl's product stays an
    integer where the doubles have rounded; the last four are drawn at the size that
    brings the score near zero."""
    if 0 in fixed:
        raise ValueError(f"{kmer}: another boundary k-mer already wrote its first entry")
    consensus = consensus_factor(kmer)
    rng = random.Random(f"{seed}/boundary/{kmer}")
    drawn = [place for place in range(5, 9) if place not in fixed]
    for _ in range(2000):
        sc = [0] + [rng.randint(2, 9) for _ in range(4)] + [1] * 4
        for place, value in fixed.items():
            sc[place] = value
        others = float(reduce(perl_multiply, sc[1:5]))
        # The divisor that puts a first entry of about 2**52.5 at a score of zero, shared
        # out over the entries of the divisor still to draw.
        bits = math.log2(consensus * others) + 52.5 - math.log2(float(reduce(perl_multiply, sc[5:])))
        for place in drawn:
            sc[place] = max(1, int(2 ** (bits / len(drawn) + rng.uniform(-0.5, 0.5))))
        divisor = float(reduce(perl_multiply, sc[5:]))
        for step in range(-BOUNDARY_SPAN, BOUNDARY_SPAN):
            boundary = (step + 0.5) / 1000.0
            nearest = round(2.0**boundary * divisor / (consensus * others))
            if nearest - 2 < _DOUBLE_INTEGERS // 2 or nearest + 2 >= _DOUBLE_INTEGERS:
                continue
            for first in range(nearest - 2, nearest + 3):
                sc[0] = first
                as_perl = _log2_score(consensus, ratio_as_perl(sc))
                in_doubles = _log2_score(consensus, ratio_in_doubles(sc))
                if (
                    format_score(as_perl) != format_score(in_doubles)
                    and abs(as_perl - boundary) > BOUNDARY_CLEARANCE
                    and abs(in_doubles - boundary) > BOUNDARY_CLEARANCE
                ):
                    return list(sc)
    raise ValueError(f"{kmer}: no entries found at a rounding boundary")


def with_boundary_kmers(files: dict[str, list[str]], seed: int, kmers: list[str]) -> None:
    """Rewrites, in `files`, the nine lines each 23-mer of `kmers` reads with its
    `boundary_entries`; a line an earlier k-mer wrote keeps its number."""
    written: dict[tuple[str, int], int] = {}
    for kmer in kmers:
        if len(kmer) != 23 or set(kmer) - set(BASES):
            raise ValueError(f"a boundary k-mer is a 23-mer over ACGT, got {kmer!r}")
        lines = acceptor_lines(kmer)
        fixed = {place: written[line] for place, line in enumerate(lines) if line in written}
        for line, value in zip(lines, boundary_entries(seed, kmer, fixed)):
            written[line] = value
            files[line[0]][line[1]] = str(value)


def write_directory(
    out: Path,
    seed: int,
    force: bool,
    plain_lines: bool = False,
    whole_numbers: bool = False,
    boundary_kmers: list[str] | None = None,
) -> dict:
    """Writes the directory under `out` and returns the summary: the seed, per table its line
    count and the sha256 of its bytes and, per other file, its size and sha256. Under
    `plain_lines` every number is its five-place decimal and no other file is written;
    under `whole_numbers` every number is a whole number, the 23-mers of `boundary_kmers`
    read the entries `boundary_entries` finds, and no other file is written. An existing
    file stops the write unless `force`."""
    if whole_numbers:
        files = whole_number_lines(seed)
        with_boundary_kmers(files, seed, boundary_kmers or [])
    else:
        files = file_lines(seed) if plain_lines else spelled_lines(seed)
    others = {} if plain_lines or whole_numbers else other_files(seed)
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
        return _log2_score(consensus, table_score)

    def score3(self, seq: str) -> float:
        """`score3` (MaxEntScan.pm:698-726) of a 23-mer: the consensus factor of positions
        18 and 19 times the maximum-entropy ratio of the other 21 bases, as log2."""
        seq = seq.upper()
        if len(seq) != 23:
            raise ValueError(f"score3 takes a 23-mer, got {seq!r}")
        return _log2_score(consensus_factor(seq), ratio_as_perl(self.entries(seq)))

    def score3_in_doubles(self, seq: str) -> float:
        """`score3` with the quotient of MaxEntScan.pm:860 computed in doubles."""
        seq = seq.upper()
        return _log2_score(consensus_factor(seq), ratio_in_doubles(self.entries(seq)))

    def entries(self, seq: str) -> list[int | float]:
        """The nine entries a 23-mer's `score3` reads, in the plugin's order."""
        return [
            table[line] for table, (_, line) in zip(self.acceptor, acceptor_lines(seq))
        ]

    def score(self, kmer: str) -> float:
        """`score5` of a 9-mer, `score3` of a 23-mer."""
        if len(kmer) == 9:
            return self.score5(kmer)
        if len(kmer) == 23:
            return self.score3(kmer)
        raise ValueError(f"a k-mer is scored at length 9 or 23, got {kmer!r}")


def _lines(path: Path) -> list[str]:
    return path.read_text(encoding="ascii").splitlines()


def _numbers(path: Path) -> list[int | float]:
    return [perl_number(line) for line in _lines(path)]


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
    parser.add_argument(
        "--whole-numbers",
        action="store_true",
        help="write each number as a whole number and no file beside the eleven tables",
    )
    parser.add_argument(
        "--boundary-kmers",
        nargs="+",
        metavar="KMER",
        help="23-mers whose entries --whole-numbers writes at a rounding boundary",
    )
    parser.add_argument("--dir", type=Path, help="a written directory to score against")
    parser.add_argument("--score", nargs="+", metavar="KMER", help="k-mers to score")
    parser.add_argument("--check-kmers", type=Path, help="k-mers and pairs to check")
    args = parser.parse_args(argv)

    if args.out is not None:
        if args.seed is None:
            parser.error("--out needs --seed")
        if args.whole_numbers and args.plain_lines:
            parser.error("--whole-numbers and --plain-lines name two spellings")
        if args.boundary_kmers and not args.whole_numbers:
            parser.error("--boundary-kmers needs --whole-numbers")
        try:
            summary = write_directory(
                args.out,
                args.seed,
                args.force,
                args.plain_lines,
                args.whole_numbers,
                args.boundary_kmers,
            )
        except (FileExistsError, ValueError) as e:
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
