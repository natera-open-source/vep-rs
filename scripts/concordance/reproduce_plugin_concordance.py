#!/usr/bin/env python3
"""Reproduce plugin concordance rows: run both engines on a run of plugin_gate_runs.toml and score it.

usage: reproduce_plugin_concordance.py (--run NAME [--set NAME] | --row PLUGIN --assembly ASSEMBLY | --list)
           --work-dir DIR [--vep BINARY] [--input FILE] [--reference FILE | --run-ensembl-vep]
           [--steps check,ensembl-vep,utrannotator-sweep,vep-rs,score] [--sweep-input FILE] [--overwrite-sweep]
           [--pairing-tool BINARY] [--dry-run] [named inputs]

plugin_gate_runs.toml, beside this script, lists every run a release scores its supported plugins
on and the rows of scripts/parity/plugin_gate_rows.py each run scores. `--run` names a whole run by
its `family` or a corpus by its directory name (every set of it, or the one `--set` names); `--row`
names a derived row by its plugin label and assembly and selects the run or the corpus sets it
stands on; `--list` prints every run with its rows. A whole run takes its input from `--input`: the
genome-wide ClinVar file of the run's dataset, or one chromosome's file of a population dataset for
a whole-genome cell under the same flags, which `--dataset` then names in the printed rows.

Steps, in order, for each selected run:

- `check`: every plugin data file the run names stands under `--plugin-data` with the index files
  the run lists, and every other input the run names exists. Names alone are read.
- `ensembl-vep`, under `--run-ensembl-vep`: Ensembl VEP writes the reference. A corpus set runs
  through generate_reference_output.sh on a copy of the corpus in the work directory. A whole run
  runs in a container of this script's own, the image named by the digest that
  generate_reference_output.sh pins, with the run's `ensembl_vep` words; `--perl-env NAME=VALUE`
  sets a variable in it, `--ensembl-vep-fork` replaces the recorded fork count, and a `--tool-dir`
  holding a program the run lists under `ensembl_vep_tools` goes first on the container's PATH.
  Without `--run-ensembl-vep` the reference is the file `--reference` names, or for a corpus set
  the corpus's committed expected output.
- `utrannotator-sweep`, only where `--steps` names it: Ensembl VEP writes, into the directory
  `--utrannotator-sweep` names, every file `sweep-seed<seed>.txt` the run's `score` words name.
  Each is written by one container of the pinned image with no network, `PERL_HASH_SEED` set to
  the seed of the file's name and `PERL_PERTURB_KEYS=0`: the run's `ensembl_vep` words in one
  process, with UTRAnnotator alone of the plugins, the header lines kept and no statistics file.
  The runs read `--sweep-input`, by default the run's `--input`: the whole input, or the records
  behind the UTRAnnotator residual rows of a first score, the rows the classes read the sweep on.
  A file that stands is written over only under `--overwrite-sweep`. A run whose `score` words
  name no such file gets no command and is named on standard error.
- `vep-rs`: `--vep` (default `vep` on PATH; a relative path is read from the working directory)
  runs with the run's `vep_rs` words, a corpus set inside its corpus directory.
- `score`: score_plugin_tuples.py scores vep-rs's output against the reference under the class
  table beside it, through `--pairing-tool` where the built tools/vep-plugin-pairs is given. Without
  the tool the script's own pairing is used only where the reference and vep-rs's output each hold
  at most 20,000,000 data lines; a larger output is refused, since only the tool pairs it in
  reasonable time.

`--steps` runs a subset where the caller holds the outputs of the steps left out: vep-rs's output
at `--engine-output` (default `<work-dir>/<run>/vep-rs.default.txt`) for `score` alone, the
reference for everything after `ensembl-vep`. Without `--steps` every step runs but
`utrannotator-sweep`, whose files a caller who holds them names with `--utrannotator-sweep`.
A whole run is asked for `--input` only where a named step reads it: `ensembl-vep` under
`--run-ensembl-vep`, `vep-rs`, `score`, and `utrannotator-sweep` without `--sweep-input`. It is
asked for `--reference` or `--run-ensembl-vep` only where `score` is named.

Output. For every run scored, standard output carries one JSON line per plugin and projection in
the form of the plugin record's entries (`plugin`, `assembly`, `dataset` and `family` or `corpus`
and `expected_set`, `engine`, `projection`, `scored_on`, the raw counts `perl`, `rust` and
`intersection`, the adjusted counts, `masked`, `residual_rows`, `text_only_differences`), with
`raw_f1` and `adjusted_f1` at nine decimal places. The record's own entries also carry the digest
of a release's run manifest, which these lines leave out. A run whose plugin populates no field on
either side, as a corpus's `base` set that loads no plugin, adds no line and is named on standard
error. Standard error carries one line per derived row, summed over the row's sets. The scorer's report is written to `--report`, or to
`<work-dir>/<run>/report.json`, and its residual files to `<work-dir>/<run>/residuals`, which a
score empties first. `--rows-out FILE` writes the JSON lines to a file as well.

`--dry-run` prints every command it would run, in order, one shell line each, and runs nothing;
it checks no file and needs every named input the selected steps read.

Named inputs. A word of plugin_gate_runs.toml in braces is filled from: `--input`,
`--vep-rs-cache`, `--fasta`, `--plugin-data`, `--loftee`, `--utrannotator-sweep`, `--generated`,
`--fork`, and for Ensembl VEP `--ensembl-vep-cache`, `--ensembl-vep-dir` and `--plugin-modules`.
A corpus is read from `--corpus-root` (default tests/golden/116 of this tree).

Exit 0 when every selected step ran; exit 1 with one `ERROR: [reproduce_plugin_concordance] ...`
line on the first refusal: a run or row the data file does not hold, a named input not given, a
missing file, a step that cannot start or ends non-zero, an output too large for the script's own
pairing, a file of the sweep that stands, a score with a declared class that was not applied.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tomllib
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

COMPONENT = "reproduce_plugin_concordance"
HERE = Path(__file__).resolve().parent
TREE = HERE.parent.parent
RUNS_FILE = HERE / "plugin_gate_runs.toml"
GENERATOR = HERE / "generate_reference_output.sh"
SCORER = HERE / "score_plugin_tuples.py"
STEPS = ("check", "ensembl-vep", "utrannotator-sweep", "vep-rs", "score")
SWEEP_STEP = "utrannotator-sweep"
DEFAULT_STEPS = tuple(step for step in STEPS if step != SWEEP_STEP)
SWEEP_PLUGIN = "UTRAnnotator"
SWEEP_OUTPUT_RE = re.compile(r"\{utrannotator_sweep\}/(sweep-seed(.+)\.txt)")
ENGINE = "vep-rs"
BASE_SET = "base"
PLACEHOLDER_RE = re.compile(r"\{([a-z_]+)\}")
PLUGIN_DATA_PATH_RE = re.compile(r"\{plugin_data\}/[^,]+")
IMAGE_RE = re.compile(r'^PERL_IMAGE="([^"$]+)"$', re.MULTILINE)
DIGEST_RE = re.compile(r'^PERL_IMAGE_DIGEST="\$\{VEP_PERL_IMAGE_DIGEST:-(sha256:[0-9a-f]{64})\}"$', re.MULTILINE)
PYTHON_PAIRING_LINE_LIMIT = 20_000_000
NINE_PLACES = Decimal("0.000000001")
RAW_COUNTS = ("intersection", "perl", "rust")
ADJUSTED_COUNTS = ("adj_intersection", "adj_perl", "adj_rust")
# Where a whole run's inputs stand inside the container Ensembl VEP runs in.
CONTAINER = {
    "ensembl_vep_cache": "/work/cache",
    "plugin_data": "/work/plugin_data",
    "loftee": "/work/loftee",
    "plugin_modules": "/work/plugins",
    "ensembl_vep_dir": "/work/dir",
}
CONTAINER_OUT = "/work/perl_out"
CONTAINER_TOOLS = "/work/tools"
# Where a seeded run's inputs stand inside its container: the records' file, and the directories
# of the cache, of the FASTA and of UTRAnnotator's data file. It writes into SWEEP_OUT.
SWEEP_RECORDS = "/work/in/input.vcf"
SWEEP_CACHE = "/work/cache"
SWEEP_FASTA = "/work/ref"
SWEEP_DATA = "/work/plugin_data"
SWEEP_OUT = "/work/out"
# The named inputs a caller gives, by the option that carries each.
INPUT_OPTIONS = {
    "input": "--input",
    "fork": "--fork",
    "vep_rs_cache": "--vep-rs-cache",
    "fasta": "--fasta",
    "plugin_data": "--plugin-data",
    "loftee": "--loftee",
    "utrannotator_sweep": "--utrannotator-sweep",
    "generated": "--generated",
    "ensembl_vep_cache": "--ensembl-vep-cache",
    "ensembl_vep_dir": "--ensembl-vep-dir",
    "plugin_modules": "--plugin-modules",
}
# The named inputs this script fills itself.
DERIVED_INPUTS = ("output", "reference", "corpus", "home", "ensembl_vep_output")


class Refusal(Exception):
    """A refusal, printed as the one error line."""


def load_runs(path: Path = RUNS_FILE) -> tuple[list[dict], dict]:
    """(the runs of the data file, its `inputs` table)."""
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise Refusal(f"{path.name} cannot be read: {error}") from error
    runs, inputs = data.get("run"), data.get("inputs")
    if not isinstance(runs, list) or not runs or not isinstance(inputs, dict):
        raise Refusal(f"{path.name} holds no `run` list and `inputs` table")
    return runs, inputs


def steps_of(args: argparse.Namespace) -> tuple[str, ...]:
    """The steps `--steps` names, in the caller's order."""
    return tuple(step.strip() for step in args.steps.split(",") if step.strip())


def run_key(run: dict) -> str:
    """A run's name: its family, or its corpus and set."""
    return run["family"] if "family" in run else f"{run['corpus']}/{run['set']}"


def plugin_of(label: str) -> str:
    """The plugin a row's label names; a label `X (...)` is a row of the plugin `X`."""
    return label.split(" (")[0]


def select(runs: list[dict], args: argparse.Namespace) -> list[dict]:
    """The runs `--run` or `--row` names, in the data file's order."""
    if args.run:
        found = [run for run in runs if args.run in (run.get("family"), run.get("corpus"))]
        what = f"the run {args.run}"
    else:
        if not args.assembly:
            raise Refusal("--row takes --assembly beside it")
        found = [run for run in runs if any(row[0] == args.row and row[1] == args.assembly for row in run["rows"])]
        what = f"the row {args.row} on {args.assembly}"
    if args.set:
        found = [run for run in found if run.get("set") == args.set]
        what = f"{what}, set {args.set}"
    if not found:
        raise Refusal(f"{RUNS_FILE.name} holds no entry for {what}; --list prints every run")
    return found


def pinned_image() -> str:
    """The image Ensembl VEP runs from, `<name>@<digest>`, read from the one line of
    generate_reference_output.sh that pins each."""
    text = GENERATOR.read_text(encoding="utf-8")
    names, digests = IMAGE_RE.findall(text), DIGEST_RE.findall(text)
    if len(names) != 1 or len(digests) != 1:
        raise Refusal(f"{GENERATOR.name} does not pin one image and one digest")
    return f"{names[0].rsplit(':', 1)[0]}@{digests[0]}"


def fill(word: str, values: dict[str, str], key: str) -> str:
    """`word` with every named input replaced; refused when the caller gave none for a name."""

    def value(match: re.Match) -> str:
        name = match.group(1)
        if values.get(name) is None:
            option = INPUT_OPTIONS.get(name)
            hint = f"pass {option}" if option else "this script fills no such input"
            raise Refusal(f"the run {key} names {{{name}}}: {hint}")
        return values[name]

    return PLACEHOLDER_RE.sub(value, word)


def program(word: str) -> str:
    """The program a word names: a path as an absolute one, a bare name as PATH resolves it."""
    if os.sep in word:
        return str(Path(word).resolve())
    return shutil.which(word) or word


def f1_nine(matched: int, reference: int, engine: int, exact: bool) -> str | None:
    """2I/(A+B) at nine decimal places, rounded half up, or None where neither side holds a tuple.
    Under `exact` a figure that rounds to 1 from unequal counts reads as the last figure below it."""
    if reference + engine == 0:
        return None
    value = (Decimal(2 * matched) / Decimal(reference + engine)).quantize(NINE_PLACES, rounding=ROUND_HALF_UP)
    if exact and value == 1 and not matched == reference == engine:
        value = 1 - NINE_PLACES
    return format(value, "f")


def record_entries(run: dict, report: dict, dataset: str | None = None) -> list[dict]:
    """The report's scores of the run's rows as entries of the plugin record, with their F1; a whole
    run's entries name `dataset` where one is given, the run's own dataset otherwise. A score with a
    declared class under `not_applied` holds no adjusted count and is refused."""
    entries = []
    for score in report.get("scores", []):
        if score.get("engine") != ENGINE:
            continue
        rows = [row for row in run["rows"] if plugin_of(row[0]) == score["source"]]
        for label, assembly, scored_on, target in rows:
            if score.get("not_applied"):
                where = f"{dataset or target} ({run['family']})" if "family" in run else f"{target}/{run['set']}"
                classes = ", ".join(f"{klass} ({count})" for klass, count in score["not_applied"].items())
                raise Refusal(f"the score of {label} on {where} carries no adjusted count: declared and not applied: {classes}")
            entry = {"plugin": label, "assembly": assembly}
            if "family" in run:
                entry.update(dataset=dataset or target, family=run["family"])
            else:
                entry.update(corpus=target, expected_set=run["set"])
            entry.update(engine=ENGINE, projection=score["projection"], scored_on=scored_on)
            for name in ("perl", "rust", "intersection", "adj_perl", "adj_rust", "adj_intersection"):
                entry[name] = score[name]
            entry.update(
                masked=score["masked"],
                residual_rows=score["residual_rows"],
                text_only_differences=score["text_only_differences"],
                raw_f1=f1_nine(*(score[name] for name in RAW_COUNTS), exact=False),
                adjusted_f1=f1_nine(*(score[name] for name in ADJUSTED_COUNTS), exact=True),
            )
            entries.append(entry)
    return entries


def row_summaries(entries: list[dict]) -> list[str]:
    """One line per derived row: its per-feature counts summed over the row's sets, and its F1."""
    sums: dict[tuple, dict[str, int]] = {}
    for entry in entries:
        if entry["projection"] != "feature":
            continue
        key = (entry["plugin"], entry["assembly"], entry.get("dataset") or entry.get("corpus"))
        counts = sums.setdefault(key, dict.fromkeys((*RAW_COUNTS, *ADJUSTED_COUNTS, "sets"), 0))
        for name in (*RAW_COUNTS, *ADJUSTED_COUNTS):
            counts[name] += entry[name]
        counts["sets"] += 1
    lines = []
    for (label, assembly, target), counts in sums.items():
        raw = f1_nine(*(counts[name] for name in RAW_COUNTS), exact=False)
        adjusted = f1_nine(*(counts[name] for name in ADJUSTED_COUNTS), exact=True)
        lines.append(
            f"{label} {assembly} on {target}: raw F1 {raw}, adjusted F1 {adjusted}; per-feature tuples: "
            f"Ensembl VEP {counts['perl']}, vep-rs {counts['rust']}, shared {counts['intersection']}; "
            f"adjusted {counts['adj_perl']}, {counts['adj_rust']}, {counts['adj_intersection']}; "
            f"{counts['sets']} scored set(s)"
        )
    return lines


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def data_lines_exceed(path: Path, limit: int) -> bool:
    """Whether a default-format output holds more than `limit` data lines; the read stops there."""
    opener = gzip.open if path.suffix == ".gz" else open
    count = 0
    with opener(path, "rb") as handle:
        for line in handle:
            if not line.startswith(b"#"):
                count += 1
                if count > limit:
                    return True
    return False


class Plan:
    """The commands of one run, in order, and the files its steps read and write."""

    def __init__(self, run: dict, args: argparse.Namespace, inputs: dict):
        self.run, self.args = run, args
        self.key = run_key(run)
        self.whole = "family" in run
        self.work = Path(args.work_dir).resolve() / self.key
        unknown = sorted(
            {name for field in ("vep_rs", "ensembl_vep", "score", "vep_rs_env", "data_indexes") for word in run.get(field, [])
             for name in PLACEHOLDER_RE.findall(word)} - set(inputs)
        )
        if unknown:
            raise Refusal(f"the run {self.key} names {{{unknown[0]}}}, which {RUNS_FILE.name} lists under no input")
        # The steps that read a whole run's input; the sweep reads its own records where given.
        self.steps = steps_of(args)
        annotates = bool({"vep-rs", "score"} & set(self.steps)) or ("ensembl-vep" in self.steps and args.run_ensembl_vep)
        if self.whole and annotates and not args.input and not args.list:
            raise Refusal(f"the run {self.key} annotates a genome-wide file: pass --input")
        if self.whole and SWEEP_STEP in self.steps and not (args.sweep_input or args.input):
            raise Refusal(f"the {SWEEP_STEP} step of {self.key} reads a file of records: pass --sweep-input or --input")
        if not self.whole and args.input:
            raise Refusal(f"the run {self.key} annotates its corpus's own file; --input is for a whole run")
        if not self.whole and args.dataset:
            raise Refusal(f"the run {self.key} is scored on its corpus; --dataset names a whole run's input")
        self.vep = program(args.vep)
        self.dataset = args.dataset or run.get("dataset")
        self.corpus = None if self.whole else Path(args.corpus_root).resolve() / run["corpus"]
        self.output = Path(args.engine_output).resolve() if args.engine_output else self.work / "vep-rs.default.txt"
        self.report = Path(args.report).resolve() if args.report else self.work / "report.json"
        self.ensembl_out = self.work / "ensembl-vep"
        self.corpus_copy = self.work / "corpus"
        self.values: dict[str, str | None] = {name: getattr(args, name) for name in INPUT_OPTIONS}
        for name in ("input", "vep_rs_cache", "fasta", "plugin_data", "loftee", "utrannotator_sweep", "generated",
                     "ensembl_vep_cache", "plugin_modules"):
            if self.values[name] is not None:
                self.values[name] = str(Path(self.values[name]).resolve())
        self.values["fork"] = str(args.fork)
        self.values["ensembl_vep_dir"] = str(Path(args.ensembl_vep_dir).resolve() if args.ensembl_vep_dir else self.work / "ensembl-vep-dir")
        self.values["output"] = str(self.output)
        self.values["home"] = str(self.work / "home")
        self.values["corpus"] = None if self.corpus is None else str(self.corpus)
        self.reference = self.reference_path()
        self.values["reference"] = None if self.reference is None else str(self.reference)

    def reference_path(self) -> Path | None:
        """Ensembl VEP's output for the run; None for a whole run whose named steps read none."""
        if self.args.reference:
            return Path(self.args.reference).resolve()
        if self.args.run_ensembl_vep:
            if self.whole:
                return self.ensembl_out / "output.txt"
            return self.corpus_copy / Path(self.run["reference"].replace("{corpus}/", ""))
        if self.whole:
            if "score" not in self.steps:
                return None
            raise Refusal(f"the run {self.key} needs Ensembl VEP's output: pass --reference or --run-ensembl-vep")
        return Path(fill(self.run["reference"], self.values, self.key))

    def filled(self, field: str) -> list[str]:
        return [fill(word, self.values, self.key) for word in self.run.get(field, [])]

    def vep_rs_command(self) -> dict:
        env = dict(word.split("=", 1) for word in self.filled("vep_rs_env"))
        return {
            "argv": [self.vep, *self.filled("vep_rs")],
            "env": env,
            "cwd": self.corpus,
            "stdout": self.work / "vep-rs.stdout",
            "stderr": self.work / "vep-rs.stderr",
        }

    def score_command(self) -> dict:
        vep_rs = self.vep_rs_command()["argv"]
        argv = [
            sys.executable, str(SCORER), "--reference", str(self.reference), "--engine", f"{ENGINE}={self.output}",
            f"--run-command={shlex.join(vep_rs)}", *self.filled("score"),
        ]
        if Path(self.vep).is_file() and not self.args.dry_run:
            argv += ["--engine-binary", f"{ENGINE}={sha256_of(Path(self.vep))}"]
        if self.args.engine_commit:
            argv += ["--engine-commit", f"{ENGINE}={self.args.engine_commit}"]
        if self.args.class_table:
            argv += ["--class-table", str(Path(self.args.class_table).resolve())]
        if self.args.classes is not None:
            argv += ["--classes", self.args.classes]
        argv += ["--work-dir", str(self.work / "scratch"), "--report", str(self.report), "--residuals", str(self.work / "residuals")]
        if self.args.sort_memory:
            argv += ["--sort-memory", self.args.sort_memory]
        if self.args.sort_threads:
            argv += ["--sort-threads", str(self.args.sort_threads)]
        if self.args.pairing_tool:
            argv += ["--pairing-tool", str(Path(self.args.pairing_tool).resolve())]
            if self.args.pairing_threads:
                argv += ["--pairing-threads", str(self.args.pairing_threads)]
        return {"argv": argv, "env": {"PYTHONDONTWRITEBYTECODE": "1"}, "cwd": None,
                "stdout": self.work / "score.stdout", "stderr": self.work / "score.stderr",
                "clears": self.work / "residuals"}

    def tool_dirs(self) -> list[Path]:
        """The `--tool-dir` directories that hold a program the run lists as lacking in the image."""
        wanted = self.run.get("ensembl_vep_tools", [])
        given = [Path(directory).resolve() for directory in self.args.tool_dir]
        if self.args.dry_run:
            return given if wanted else []
        holders = []
        for tool in wanted:
            holder = next((directory for directory in given if (directory / tool).is_file()), None)
            if holder is None:
                raise Refusal(f"the run {self.key} calls {tool}, which the image lacks: pass --tool-dir with a directory holding it")
            if holder not in holders:
                holders.append(holder)
        return holders

    def ensembl_vep_commands(self) -> list[dict]:
        if not self.whole:
            return self.generator_commands()
        words = list(self.run["ensembl_vep"])
        if self.args.ensembl_vep_fork:
            words[words.index("--fork") + 1] = str(self.args.ensembl_vep_fork)
        named = {name for word in words for name in PLACEHOLDER_RE.findall(word)}
        source = Path(fill("{input}", self.values, self.key))
        in_mount = "/work/input.vcf.gz" if source.name.endswith(".gz") else "/work/input.vcf"
        mounts = [f"{source}:{in_mount}:ro", f"{self.ensembl_out}:{CONTAINER_OUT}"]
        inside = {"input": in_mount, "ensembl_vep_output": f"{CONTAINER_OUT}/output.txt"}
        if "fasta" in named:
            fasta = Path(fill("{fasta}", self.values, self.key))
            ref = "/work/ref.fa.gz" if fasta.name.endswith(".gz") else "/work/ref.fa"
            mounts += [f"{fasta}:{ref}:ro", f"{fasta}.fai:{ref}.fai:ro"]
            if fasta.name.endswith(".gz"):
                mounts.append(f"{fasta}.gzi:{ref}.gzi:ro")
            inside["fasta"] = ref
        for name, target in CONTAINER.items():
            if name in named:
                mode = "" if name == "ensembl_vep_dir" else ":ro"
                mounts.append(f"{fill('{' + name + '}', self.values, self.key)}:{target}{mode}")
                inside[name] = target
        tools = self.tool_dirs()
        mounts += [f"{directory}:{CONTAINER_TOOLS}/{index}:ro" for index, directory in enumerate(tools)]
        argv = ["docker", "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}"]
        for mount in mounts:
            argv += ["-v", mount]
        for word in (*self.run.get("ensembl_vep_env", []), *self.args.perl_env):
            argv += ["-e", word]
        argv.append(pinned_image())
        if tools:
            # The image lacks the program: its directory goes first on PATH.
            path = ":".join(f"{CONTAINER_TOOLS}/{index}" for index in range(len(tools)))
            argv += ["sh", "-c", 'PATH="$1:$PATH"; shift; exec "$@"', "sh", path]
        argv += ["vep", *(fill(word, inside, self.key) for word in words)]
        return [{"argv": argv, "env": {}, "cwd": None, "stdout": self.work / "ensembl-vep.stdout",
                 "stderr": self.work / "ensembl-vep.stderr",
                 "makes": [self.ensembl_out, Path(self.values["ensembl_vep_dir"])]}]

    def sweep_commands(self) -> list[dict]:
        """One container for every seeded output the run's `score` words name, in their order:
        the run's `ensembl_vep` words under the seed of the output's name, in one process, with
        UTRAnnotator alone of the plugins, the header lines kept and no statistics file."""
        outputs = [match for match in map(SWEEP_OUTPUT_RE.fullmatch, self.run.get("score", [])) if match]
        if not outputs:
            return []
        records = Path(self.args.sweep_input).resolve() if self.args.sweep_input else Path(fill("{input}", self.values, self.key))
        cache = Path(fill("{ensembl_vep_cache}", self.values, self.key))
        sweep = Path(fill("{utrannotator_sweep}", self.values, self.key))
        fasta = Path(fill("{fasta}", self.values, self.key))
        inside = {
            "input": SWEEP_RECORDS + (".gz" if records.name.endswith(".gz") else ""),
            "ensembl_vep_cache": SWEEP_CACHE,
            "fasta": f"{SWEEP_FASTA}/{fasta.name}",
        }
        mounts = [f"{cache}:{SWEEP_CACHE}:ro", f"{records}:{inside['input']}:ro", f"{sweep}:{SWEEP_OUT}",
                  f"{fasta.parent}:{SWEEP_FASTA}:ro"]
        reads = [cache / "homo_sapiens" / f"116_{self.run['assembly']}", records, fasta, Path(f"{fasta}.fai")]
        words: list[str] = []
        recorded = iter(self.run["ensembl_vep"])
        for word in recorded:
            if word == "--plugin":
                plugin = next(recorded)
                if plugin.split(",")[0] != SWEEP_PLUGIN:
                    continue
                data = PLUGIN_DATA_PATH_RE.search(plugin)
                if data:
                    file = Path(fill(data.group(0), self.values, self.key))
                    mounts.append(f"{file.parent}:{SWEEP_DATA}:ro")
                    reads.append(file)
                    plugin = plugin.replace(data.group(0), f"{SWEEP_DATA}/{file.name}")
                words += [word, plugin]
            elif word == "--fork":
                words += [word, "1"]
                next(recorded)
            elif word == "--stats_text":
                words.append("--no_stats")
            elif word != "--no_headers":
                words.append(word)
        if "--plugin" not in words:
            raise Refusal(f"the run {self.key} names seeded outputs and loads no {SWEEP_PLUGIN}")
        commands = []
        for output in outputs:
            name, seed = output.groups()
            inside["ensembl_vep_output"] = f"{SWEEP_OUT}/{name}"
            argv = ["docker", "run", "--rm", "--name", f"seeded-{seed}-{self.run['assembly']}", "--network", "none",
                    "--user", f"{os.getuid()}:{os.getgid()}",
                    "-e", "HOME=/tmp", "-e", f"PERL_HASH_SEED={seed}", "-e", "PERL_PERTURB_KEYS=0"]
            for mount in mounts:
                argv += ["-v", mount]
            # The image holds its plugin modules under /plugins.
            argv += ["-w", SWEEP_OUT, pinned_image(), "vep", *(fill(word, inside, self.key) for word in words),
                     "--dir_plugins", "/plugins"]
            stem = name.removesuffix(".txt")
            commands.append({"argv": argv, "env": {}, "cwd": None, "stdout": self.work / f"{stem}.stdout",
                             "stderr": self.work / f"{stem}.stderr", "makes": [sweep], "writes": sweep / name,
                             "reads": reads})
        return commands

    def generator_commands(self) -> list[dict]:
        if any("{generated}" in word for word in self.run["ensembl_vep"]):
            raise Refusal(
                f"the run {self.key} reads a directory its manifest generates, which {GENERATOR.name} cannot name; "
                "its reference is the corpus's committed expected output"
            )
        argv = ["bash", str(GENERATOR), "--corpus", str(self.corpus_copy), "--perl-cache-dir",
                fill("{ensembl_vep_cache}", self.values, self.key)]
        if self.run["set"] != BASE_SET:
            argv += ["--flag-set", self.run["set"]]
        if self.needs_plugin_modules():
            argv += ["--dir-plugins", fill("{plugin_modules}", self.values, self.key)]
        for word in (*self.run.get("ensembl_vep_env", []), *self.args.perl_env):
            argv += ["--perl-env", word]
        argv.append("--force")
        copy = {"argv": ["cp", "-R", str(self.corpus), str(self.corpus_copy)], "env": {}, "cwd": None,
                "clears": self.corpus_copy}
        return [copy, {"argv": argv, "env": {}, "cwd": None, "stdout": self.work / "ensembl-vep.stdout",
                       "stderr": self.work / "ensembl-vep.stderr"}]

    def needs_plugin_modules(self) -> bool:
        """Whether generate_reference_output.sh asks for the plugin modules' directory: the manifest
        names plugin files or a flag word carries `{plugins}`."""
        if any("{plugin_modules}" in word for word in self.run["ensembl_vep"]):
            return True
        manifest = self.corpus / "manifest.json"
        if not manifest.is_file():
            return False
        data = json.loads(manifest.read_text(encoding="utf-8"))
        return bool(data.get("plugin_files"))

    def missing(self, steps: tuple[str, ...]) -> list[str]:
        """Every file or directory the selected steps read that does not exist, by name."""
        wanted: list[Path] = []
        if "vep-rs" in steps:
            if shutil.which(self.vep) is None and not Path(self.vep).is_file():
                wanted.append(Path(self.vep))
            for field in ("vep_rs", "data_indexes"):
                for word in self.run.get(field, []):
                    wanted += [Path(fill(path, self.values, self.key)) for path in PLUGIN_DATA_PATH_RE.findall(word)]
            for name in ("input", "vep_rs_cache", "fasta", "loftee", "generated", "corpus"):
                if any(f"{{{name}}}" in word for word in self.run["vep_rs"]):
                    wanted.append(Path(fill(f"{{{name}}}", self.values, self.key)))
            if any("{fasta}" in word for word in self.run["vep_rs"]):
                wanted.append(Path(fill("{fasta}", self.values, self.key) + ".fai"))
        sweep = self.sweep_commands() if SWEEP_STEP in steps else []
        if sweep:
            if shutil.which("docker") is None:
                wanted.append(Path("docker"))
            wanted += sweep[0]["reads"]
        if "score" in steps:
            written = {str(command["writes"]) for command in sweep}
            wanted += [Path(word) for word in self.filled("score") if word.startswith("/") and word not in written]
            if self.args.pairing_tool:
                wanted.append(Path(self.args.pairing_tool))
            if "vep-rs" not in steps:
                wanted.append(self.output)
        if "ensembl-vep" in steps and self.args.run_ensembl_vep:
            if shutil.which("docker") is None:
                wanted.append(Path("docker"))
            cache = Path(fill("{ensembl_vep_cache}", self.values, self.key))
            wanted.append(cache / "homo_sapiens" / f"116_{self.run['assembly']}")
            for word in self.run["ensembl_vep"]:
                wanted += [Path(fill(path, self.values, self.key)) for path in PLUGIN_DATA_PATH_RE.findall(word)]
            for name in ("loftee", "plugin_modules"):
                if any(f"{{{name}}}" in word for word in self.run["ensembl_vep"]):
                    wanted.append(Path(fill(f"{{{name}}}", self.values, self.key)))
        elif not ({"ensembl-vep"} & set(steps) and self.args.run_ensembl_vep) and "score" in steps:
            wanted.append(self.reference)
        seen: list[str] = []
        for path in wanted:
            exists = path.exists() if path.is_absolute() else shutil.which(str(path)) is not None
            if not exists and str(path) not in seen:
                seen.append(str(path))
        return seen


def shell_line(command: dict) -> str:
    """A command as one shell line: its directory, its environment words, its words, its output files."""
    line = shlex.join(command["argv"])
    if command["env"]:
        line = "env " + " ".join(f"{name}={shlex.quote(value)}" for name, value in command["env"].items()) + " " + line
    if command["cwd"]:
        line = f"cd {shlex.quote(str(command['cwd']))} && {line}"
    for mark, name in ((">", "stdout"), ("2>", "stderr")):
        if command.get(name):
            line = f"{line} {mark} {shlex.quote(str(command[name]))}"
    return line


def execute(command: dict, step: str, key: str) -> None:
    for directory in command.get("makes", []):
        directory.mkdir(parents=True, exist_ok=True)
    cleared = command.get("clears")
    if cleared is not None and cleared.exists():
        shutil.rmtree(cleared)
    stdout, stderr = command.get("stdout"), command.get("stderr")
    handles = [path.open("wb") if path else None for path in (stdout, stderr)]
    try:
        code = subprocess.run(
            command["argv"], cwd=command["cwd"], env={**os.environ, **command["env"]},
            stdout=handles[0] or sys.stderr, stderr=handles[1], check=False,
        ).returncode
    except OSError as error:
        raise Refusal(f"the {step} step of {key} cannot start {command['argv'][0]}: {error.strerror}") from error
    finally:
        for handle in handles:
            if handle:
                handle.close()
    if code != 0:
        where = f"; its standard error is in {stderr}" if stderr else ""
        raise Refusal(f"the {step} step of {key} ended {code}{where}")


def run_plan(plan: Plan, steps: tuple[str, ...]) -> list[dict]:
    """Run (or print) the selected steps of one run; the record entries of its score."""
    args = plan.args
    commands: list[tuple[str, dict]] = []
    if "ensembl-vep" in steps and args.run_ensembl_vep:
        commands += [("ensembl-vep", command) for command in plan.ensembl_vep_commands()]
    if SWEEP_STEP in steps:
        sweep = plan.sweep_commands()
        if not sweep:
            print(f"{plan.key}: the run's score words name no seeded output; the {SWEEP_STEP} step writes nothing",
                  file=sys.stderr)
        commands += [(SWEEP_STEP, command) for command in sweep]
    if "vep-rs" in steps:
        commands.append(("vep-rs", plan.vep_rs_command()))
    if "score" in steps:
        commands.append(("score", plan.score_command()))
    if args.dry_run:
        for _step, command in commands:
            print(shell_line(command))
        return []
    standing = [command["writes"] for step, command in commands if step == SWEEP_STEP and command["writes"].exists()]
    if standing and not args.overwrite_sweep:
        raise Refusal(f"the {SWEEP_STEP} step of {plan.key} would write over {standing[0]}: pass --overwrite-sweep")
    if "check" in steps:
        missing = plan.missing(steps)
        if missing:
            raise Refusal(f"the run {plan.key} reads {len(missing)} file(s) that are missing: {'; '.join(missing)}")
    plan.work.mkdir(parents=True, exist_ok=True)
    Path(plan.values["home"]).mkdir(exist_ok=True)
    plan.output.parent.mkdir(parents=True, exist_ok=True)
    plan.report.parent.mkdir(parents=True, exist_ok=True)
    for step, command in commands:
        if step == "score" and not args.pairing_tool:
            for path in (plan.reference, plan.output):
                if data_lines_exceed(path, PYTHON_PAIRING_LINE_LIMIT):
                    raise Refusal(
                        f"{path} holds more than {PYTHON_PAIRING_LINE_LIMIT} data lines: pass --pairing-tool with the "
                        "built tools/vep-plugin-pairs, which alone pairs an output of that size"
                    )
        execute(command, step, plan.key)
    if "score" not in steps:
        return []
    shutil.rmtree(plan.work / "scratch", ignore_errors=True)
    report = json.loads(plan.report.read_text(encoding="utf-8"))
    for plugin in report.get("not_populated", []):
        print(f"{plan.key}: no field of {plugin} is populated on either side; the run adds no entry for it", file=sys.stderr)
    return record_entries(plan.run, report, plan.dataset)


def list_runs(runs: list[dict]) -> None:
    for run in runs:
        rows = ", ".join(f"{label} {assembly}" for label, assembly, _scored_on, _target in run["rows"])
        print(f"{run_key(run)}\t{run['rows'][0][2]}\t{rows}")


def parse(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--run", metavar="NAME", help="a whole run's family, or a corpus's directory name")
    parser.add_argument("--set", metavar="NAME", help="one set of the corpus (`base` is the manifest's own run)")
    parser.add_argument("--row", metavar="PLUGIN", help="a derived row's plugin label; takes --assembly")
    parser.add_argument("--assembly", help="the assembly of --row")
    parser.add_argument("--list", action="store_true", help="print every run with the rows it scores")
    parser.add_argument("--work-dir", metavar="DIR", help="where outputs, reports and scratch are written")
    parser.add_argument("--vep", default="vep", metavar="BINARY", help="the vep-rs binary (default: vep on PATH)")
    parser.add_argument("--reference", metavar="FILE", help="Ensembl VEP's stored default-format output")
    parser.add_argument("--run-ensembl-vep", action="store_true", help="run Ensembl VEP for the reference")
    parser.add_argument("--steps", default=",".join(DEFAULT_STEPS),
                        help=f"the steps to run, of {', '.join(STEPS)} (default: every step but {SWEEP_STEP})")
    parser.add_argument("--sweep-input", metavar="FILE", help=f"the VCF the {SWEEP_STEP} step reads (default: --input)")
    parser.add_argument("--overwrite-sweep", action="store_true", help=f"let the {SWEEP_STEP} step write over files that stand")
    parser.add_argument("--dry-run", action="store_true", help="print every command, run nothing")
    parser.add_argument("--engine-output", metavar="FILE", help="where vep-rs's output is written or read")
    parser.add_argument("--engine-commit", metavar="COMMIT", help="the source commit of the binary, for the report")
    parser.add_argument("--report", metavar="FILE", help="where the scorer's report is written (one run)")
    parser.add_argument("--rows-out", metavar="FILE", help="write the printed JSON lines to a file as well")
    parser.add_argument("--corpus-root", default=str(TREE / "tests/golden/116"), metavar="DIR", help="the corpora's directory")
    parser.add_argument("--input", metavar="FILE", help="the VCF a whole run annotates")
    parser.add_argument("--dataset", metavar="NAME", help="the dataset of --input in the printed rows (default: the run's)")
    parser.add_argument("--fork", type=int, default=1, metavar="N", help="vep-rs's fork count (default 1)")
    parser.add_argument("--vep-rs-cache", metavar="DIR", help="the vep-rs JSON cache of the run's assembly")
    parser.add_argument("--fasta", metavar="FILE", help="the reference FASTA of the run's assembly")
    parser.add_argument("--plugin-data", metavar="DIR", help="the plugin data directory")
    parser.add_argument("--loftee", metavar="DIR", help="LoFTEE's checkout with its stand-in files")
    parser.add_argument("--utrannotator-sweep", metavar="DIR", help="Ensembl VEP's outputs under fixed Perl hash seeds")
    parser.add_argument("--generated", metavar="DIR", help="the directory a corpus manifest's `generated` command writes")
    parser.add_argument("--ensembl-vep-cache", metavar="DIR", help="the directory holding Ensembl VEP's cache")
    parser.add_argument("--ensembl-vep-dir", metavar="DIR", help="a writable directory for Ensembl VEP's --dir")
    parser.add_argument("--plugin-modules", metavar="DIR", help="the directory of Ensembl VEP's plugin modules")
    parser.add_argument("--ensembl-vep-fork", type=int, metavar="N", help="replace Ensembl VEP's recorded fork count")
    parser.add_argument("--perl-env", action="append", default=[], metavar="NAME=VALUE", help="a variable of Ensembl VEP's run")
    parser.add_argument("--tool-dir", action="append", default=[], metavar="DIR", help="a directory of programs the image lacks")
    parser.add_argument("--pairing-tool", metavar="BINARY", help="the built tools/vep-plugin-pairs")
    parser.add_argument("--pairing-threads", type=int, metavar="N", help="threads of the pairing tool")
    parser.add_argument("--sort-memory", metavar="SIZE", help="buffer of each sort process of the scorer")
    parser.add_argument("--sort-threads", type=int, metavar="N", help="threads of each sort process of the scorer")
    parser.add_argument("--class-table", metavar="TSV", help="a class table in place of the tree's")
    parser.add_argument("--classes", metavar="A,B", help="the declared classes to apply, in place of the table's own")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse(argv)
    try:
        runs, inputs = load_runs()
        if args.list:
            list_runs(runs)
            return 0
        if bool(args.run) == bool(args.row):
            raise Refusal("name one of --run, --row or --list")
        if not args.work_dir:
            raise Refusal("pass --work-dir, where outputs, reports and scratch are written")
        steps = steps_of(args)
        unknown = [step for step in steps if step not in STEPS]
        if unknown:
            raise Refusal(f"--steps names {unknown[0]!r}; the steps are {', '.join(STEPS)}")
        if args.reference and args.run_ensembl_vep:
            raise Refusal("--reference and --run-ensembl-vep both name the reference; pass one")
        for word in args.perl_env:
            if "=" not in word or word.startswith("=") or any(char.isspace() for char in word):
                raise Refusal(f"--perl-env takes NAME=VALUE, not {word!r}")
        selected = select(runs, args)
        if len(selected) > 1 and (args.report or args.engine_output or args.reference):
            raise Refusal("--report, --engine-output and --reference name one run's file; select one run with --set")
        plans = [Plan(run, args, inputs) for run in selected]
        entries: list[dict] = []
        for plan in plans:
            entries += run_plan(plan, steps)
    except Refusal as refusal:
        print(f"ERROR: [{COMPONENT}] {refusal}", file=sys.stderr)
        return 1
    text = "".join(json.dumps(entry) + "\n" for entry in entries)
    sys.stdout.write(text)
    if args.rows_out and not args.dry_run:
        Path(args.rows_out).write_text(text, encoding="utf-8")
    for line in row_summaries(entries):
        print(line, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
