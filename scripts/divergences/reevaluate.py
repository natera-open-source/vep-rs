#!/usr/bin/env python3
"""Re-run every intended divergence on Ensembl VEP and vep-rs and say whether each is still present.

usage: reevaluate.py --image <Ensembl VEP image tag or digest> --vep <vep-rs binary>
                     --perl-cache-dir <dir> [--perl-cache-dir <dir> ...] [--plugins-dir <dir>]
                     [--entry <name> ...] [--jobs <n>] [--work-dir <dir>] [--keep]
                     [--expect-all-present] [--json <file>]
       reevaluate.py --check [--reference 116.2] [--page docs/intended-divergences.md]
       reevaluate.py --list

Every entry of the category tables of ``docs/intended-divergences.md`` has a
reproduction, ``reproductions/<entry>.toml`` beside this script: the input
record or file, the flags, the cache or data the example needs as a path under
the committed corpora or a synthetic input the reproduction constructs, the
field or row where the engines differ with both engines' recorded values, the
Ensembl VEP version the reproduction was last verified against, the pinned
image digest, and the entry's input class with its basis. The schema is
``reproductions/README.md``.

A run executes each reproduction on both engines, reads the compared value from
each output and prints one verdict per entry: ``still-present`` when the engines
differ and each prints its recorded value, ``closed`` when the engines agree,
``changed`` when they differ in another way, and ``error`` when a run or a read
failed. The summary table counts the verdicts; with ``--expect-all-present`` the
exit status is 1 unless every entry is ``still-present``.

``--check`` runs nothing: it exits 0 when every entry of the page's category
tables has a reproduction that passes the schema, whose ``verified_against``
equals ``--reference`` and the page's ``Verified against`` cell, no
reproduction lacks a row, and every reproduction names its input class with a
basis: ``observed`` for a tier 1 or tier 2 entry, the basis its ``Tier and
count`` cell; any of ``observed``, ``unobserved`` and ``malformed`` for a tier
3 entry, an ``observed`` basis there reading ``seen in <the real data its page
names>: <the figure its page measures there>``. A ``malformed`` entry stands in
the page's malformed-input section and no other entry does.

Ensembl VEP runs in Docker from ``--image`` (``--ensembl-vep <command>`` runs a
local installation instead), its cache the ``<species>[_merged|_refseq]/
<version>_<assembly>`` directory found under the ``--perl-cache-dir`` roots in
order, the first that holds every chromosome the input names; its plugin
modules are ``--plugins-dir``. vep-rs is given the JSON cache the reproduction
names: a committed one, or a path under the stage directory, where its prepare
steps write a cache or hold the cache builder to refusing one. Stdlib only;
Python 3.11 or later for ``tomllib``.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import gzip
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import tomllib
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_divergence_order as cdo

REPO = HERE.parents[1]
REPRODUCTIONS = HERE / "reproductions"
PAGE = REPO / "docs" / "intended-divergences.md"
VERDICTS = ("still-present", "closed", "changed", "error")
FORMATS = {
    "default": "default.txt",
    "tab": "tab.txt",
    "vcf": "vcf.vcf",
    "json": "json.jsonl",
}
INPUT_EXTENSIONS = {
    "vcf": "vcf",
    "ensembl": "txt",
    "hgvs": "txt",
    "spdi": "txt",
    "region": "txt",
    "id": "txt",
    "guess": "txt",
}
COMPARE_KINDS = (
    "rows",
    "field",
    "count",
    "header",
    "lines",
    "warnings",
    "stderr",
    "stdout",
    "exit",
    "stats",
    "json",
    "vcf",
)
# The `where` keys each compare kind reads; a key a kind ignores is refused rather than silently matching every row.
TABLE_WHERE_KEYS = (
    "record",
    "location",
    "allele",
    "gene",
    "feature",
    "feature_type",
    "consequence",
    "pattern",
)
RECORD_WHERE_KEYS = ("record", "id", "chrom", "pos", "location", "pattern")
WHERE_KEYS = {
    "rows": TABLE_WHERE_KEYS,
    "field": TABLE_WHERE_KEYS,
    "count": TABLE_WHERE_KEYS,
    "json": (*RECORD_WHERE_KEYS, "input"),
    "vcf": (*RECORD_WHERE_KEYS, "feature", "allele"),
}
WHERE_COLUMNS = {
    "record": "Uploaded_variation",
    "location": "Location",
    "allele": "Allele",
    "gene": "Gene",
    "feature": "Feature",
    "feature_type": "Feature_type",
    "consequence": "Consequence",
}
VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
DEFAULT_COLUMNS = ["Location", "Allele", "Feature", "Consequence"]
# Any absolute path, so a value naming the host tree and one naming a container mount read alike;
# the slash of an allele string (`-/A`, `C//G`, `*/Q`) follows `-`, `/` or `*` and is not one.
PATH = re.compile(r"(?<![\w.*/-])/(?:[\w.+@%:,=~-]+/)*[\w.+@%:,=~-]+")
CONTAINER_WORK = "/work"
TOOLS = {"vep": "vep", "filter": "filter_vep", "recode": "variant_recoder"}


class ReproductionError(Exception):
    """A reproduction file that cannot be run as written."""


@dataclasses.dataclass
class Reproduction:
    name: str
    path: Path
    data: dict

    @property
    def input(self) -> dict:
        return self.data.get("input", {})

    @property
    def cache(self) -> dict:
        return self.data.get("cache", {})

    @property
    def run(self) -> dict:
        return self.data.get("run", {})

    @property
    def compare(self) -> dict:
        return self.data["compare"]

    @property
    def expected(self) -> dict:
        return self.data["expected"]


@dataclasses.dataclass
class EngineRun:
    """One engine's run: its command, exit status, output file and the files beside it."""

    engine: str
    command: list[str]
    exit_code: int
    output: Path
    stdout: str
    stderr: str

    @property
    def warnings(self) -> Path:
        return self.output.with_name(self.output.name + "_warnings.txt")

    @property
    def stats(self) -> Path:
        return self.output.with_name(self.output.name + "_summary.txt")


@dataclasses.dataclass
class Result:
    entry: str
    verdict: str
    ensembl_vep: list[str] | None = None
    vep_rs: list[str] | None = None
    detail: str = ""


def load_reproductions(
    names: list[str] | None = None, directory: Path = REPRODUCTIONS
) -> list[Reproduction]:
    """Every reproduction file, or the named ones; a name without a file, a directory without one, or a file that is not TOML is an error."""
    found = {p.stem: p for p in sorted(directory.glob("*.toml"))}
    if names:
        missing = [n for n in names if n not in found]
        if missing:
            raise ReproductionError(
                f"no reproduction for {', '.join(missing)} under {directory}"
            )
        found = {n: found[n] for n in names}
    elif not found:
        raise ReproductionError(f"no reproduction under {directory}")
    repros = []
    for name, path in found.items():
        with open(path, "rb") as fh:
            try:
                data = tomllib.load(fh)
            except (tomllib.TOMLDecodeError, UnicodeDecodeError) as exc:
                raise ReproductionError(f"{path}: {exc}") from exc
        repros.append(Reproduction(name, path, data))
    return repros


KEYS = {
    "": (
        "entry",
        "page",
        "verified_against",
        "image",
        "input_class",
        "input_class_basis",
        "notes",
        "prepare",
        "input",
        "cache",
        "run",
        "compare",
        "expected",
    ),
    "input": ("text", "path", "format", "header", "samples", "name"),
    "cache": (
        "species",
        "assembly",
        "source",
        "version",
        "json_cache",
        "corpus",
        "fasta",
        "gxf",
        "gxf_kind",
        "perl_cache",
        "chromosomes",
        "writable",
    ),
    "run": (
        "flags",
        "ensembl_vep_flags",
        "vep_rs_flags",
        "format",
        "format_flag",
        "stats",
        "buffer_size",
        "tool",
        "args",
        "stdin",
        "perl_hash_seed",
        "perl_mode",
        "assembly_flag",
        "tools",
    ),
    "compare": (
        "kind",
        "where",
        "field",
        "columns",
        "with_key",
        "pattern",
        "path",
        "section",
        "info",
        "line",
        "mask_paths",
        "mask",
        "sorted",
        "strip",
    ),
    "expected": ("ensembl_vep", "vep_rs", "ensembl_vep_exit", "vep_rs_exit"),
}


# The keys whose value is a string; a number there would reach a path or a flag as a Python object.
STRING_KEYS = (
    ("", "notes"),
    ("input", "text"),
    ("input", "path"),
    ("input", "format"),
    ("input", "name"),
    ("cache", "species"),
    ("cache", "assembly"),
    ("cache", "source"),
    ("cache", "version"),
    ("cache", "json_cache"),
    ("cache", "corpus"),
    ("cache", "fasta"),
    ("cache", "gxf"),
    ("cache", "gxf_kind"),
    ("cache", "perl_cache"),
    ("run", "format"),
    ("run", "tool"),
    ("run", "perl_hash_seed"),
    ("run", "perl_mode"),
    ("compare", "kind"),
    ("compare", "field"),
    ("compare", "pattern"),
    ("compare", "path"),
    ("compare", "section"),
    ("compare", "info"),
)
# The keys whose value is a list of strings; a bare string there would be read one character at a time.
LIST_KEYS = (
    ("", "prepare"),
    ("input", "header"),
    ("input", "samples"),
    ("cache", "chromosomes"),
    ("run", "flags"),
    ("run", "ensembl_vep_flags"),
    ("run", "vep_rs_flags"),
    ("run", "args"),
    ("compare", "columns"),
    ("compare", "mask"),
    ("run", "tools"),
)


# The fixed words of an `unobserved` entry's basis, around the format or option contract it names; a basis may go on after the closing words.
UNOBSERVED_OPENING = "well-formed under "
UNOBSERVED_CLOSING = "; not seen in the measured datasets"
# The basis of an `observed` entry of tier 3: the real data its page names, up to the first colon a space follows, then the figure the page measures there; the name opens on neither a blank nor a colon and ends on no blank, and the rest opens on no blank.
SIGHTING = re.compile(r"^seen in (?P<where>[^\s:](?:(?!: ).)*)(?<=\S): (?P<what>\S.*)$")
SIGHTING_FORM = (
    "seen in <the real data its page names>: <the figure its page measures there>"
)


def input_class_errors(d: dict) -> list[str]:
    """What is wrong with a reproduction's `input_class` and `input_class_basis`: one of the three classes, and one line of basis in the form the class asks for."""
    errors = []
    input_class, basis = d.get("input_class"), d.get("input_class_basis")
    for key, value in (("input_class", input_class), ("input_class_basis", basis)):
        if not isinstance(value, str) or not value.strip():
            errors.append(f"{key} is missing")
    if errors:
        return errors
    if input_class not in cdo.INPUT_CLASSES:
        return [
            f"input_class {input_class!r} is not one of {', '.join(cdo.INPUT_CLASSES)}"
        ]
    if basis.splitlines() != [basis]:
        return ["input_class_basis is not one line"]
    if input_class == "observed":
        basis_cell = cdo.CELL.match(basis)
        tier_1_or_2_cell = basis_cell and basis_cell.group("tier") in ("1", "2")
        if not tier_1_or_2_cell and not SIGHTING.match(basis):
            errors.append(
                f"input_class_basis of an observed entry is neither a tier 1 or tier 2 `Tier and count` cell nor `{SIGHTING_FORM}`"
            )
    elif input_class == "unobserved":
        contract, closing, _ = basis.removeprefix(UNOBSERVED_OPENING).partition(
            UNOBSERVED_CLOSING
        )
        if not (basis.startswith(UNOBSERVED_OPENING) and closing and contract.strip()):
            errors.append(
                f"input_class_basis of an unobserved entry does not read `{UNOBSERVED_OPENING}<format or option contract>{UNOBSERVED_CLOSING}`"
            )
    elif (
        cdo.CELL.match(basis)
        or SIGHTING.match(basis)
        or basis.startswith(UNOBSERVED_OPENING)
    ):
        errors.append(
            "input_class_basis of a malformed entry reads as the basis of another class, not as the rule the input breaks"
        )
    return errors


def validate(repro: Reproduction) -> list[str]:
    """The schema violations of one reproduction, as one line each."""
    d = repro.data
    errors = []
    tables = True
    for table, allowed in KEYS.items():
        section = d if table == "" else d.get(table, {})
        if not isinstance(section, dict):
            errors.append(f"{table} is not a table")
            tables = False
            continue
        for key in section:
            if key not in allowed:
                errors.append(f"unknown key {table + '.' if table else ''}{key}")
    if not tables:
        # A section that is not a table has no keys to read; the checks below would fail on it.
        return errors
    if d.get("entry") != repro.name:
        errors.append(f"entry is {d.get('entry')!r}, the file is {repro.name}.toml")
    for key in ("page", "verified_against", "image"):
        if not isinstance(d.get(key), str) or not d[key]:
            errors.append(f"{key} is missing")
    if isinstance(d.get("image"), str) and "@sha256:" not in d["image"]:
        errors.append("image is not pinned by digest")
    errors += input_class_errors(d)
    typed = True
    for table, key in STRING_KEYS:
        value = (d if table == "" else d.get(table, {})).get(key)
        if value is not None and not isinstance(value, str):
            errors.append(f"{table + '.' if table else ''}{key} is not a string")
            typed = False
    for table, key in LIST_KEYS:
        value = (d if table == "" else d.get(table, {})).get(key)
        if value is not None and (
            not isinstance(value, list) or not all(isinstance(v, str) for v in value)
        ):
            errors.append(
                f"{table + '.' if table else ''}{key} is not a list of strings"
            )
            typed = False
    if not typed:
        # The checks below read the typed values; a list where a string belongs would not hash.
        return errors
    inp = repro.input
    if not ("text" in inp) ^ ("path" in inp):
        errors.append("input needs exactly one of text and path")
    if inp.get("format", "vcf") not in INPUT_EXTENSIONS:
        errors.append(
            f"input.format {inp.get('format')!r} is not one of {', '.join(INPUT_EXTENSIONS)}"
        )
    cache = repro.cache
    run = repro.run
    tool = run.get("tool", "vep")
    if tool not in TOOLS:
        errors.append(f"run.tool {tool!r} is not one of {', '.join(TOOLS)}")
    if tool == "vep":
        if not cache.get("assembly"):
            errors.append("cache.assembly is missing")
        if not (cache.get("json_cache") or cache.get("gxf")):
            errors.append("cache needs json_cache or gxf")
        if cache.get("source", "ensembl") not in ("ensembl", "merged", "refseq"):
            errors.append("cache.source is not ensembl, merged or refseq")
    elif not run.get("args"):
        errors.append("a tool run needs run.args")
    compare = d.get("compare")
    if not isinstance(compare, dict):
        errors.append("compare is missing")
    else:
        kind = compare.get("kind")
        if kind not in COMPARE_KINDS:
            errors.append(
                f"compare.kind {kind!r} is not one of {', '.join(COMPARE_KINDS)}"
            )
        needs_field = kind == "field" or (kind == "vcf" and not compare.get("line"))
        if needs_field and not compare.get("field"):
            errors.append(f"compare.kind {kind} needs compare.field")
        if kind in (
            "header",
            "lines",
            "warnings",
            "stderr",
            "stdout",
        ) and not compare.get("pattern"):
            errors.append(f"compare.kind {kind} needs compare.pattern")
        if kind == "json" and not compare.get("path"):
            errors.append("compare.kind json needs compare.path")
        if kind == "stats" and not compare.get("section"):
            errors.append("compare.kind stats needs compare.section")
        where = compare.get("where", {})
        if not isinstance(where, dict):
            errors.append("compare.where is not a table")
            where = {}
        allowed_where = WHERE_KEYS.get(kind, ())
        if kind == "vcf" and compare.get("line"):
            # A whole data line is matched at the record level alone; the CSQ-entry keys would select nothing.
            allowed_where = RECORD_WHERE_KEYS
        unknown = [k for k in where if k not in allowed_where]
        if unknown:
            errors.append(
                f"compare.where has keys {', '.join(unknown)} that compare.kind {kind} does not read"
            )
        mask = compare.get("mask", [])
        patterns = [
            ("compare.pattern", compare.get("pattern")),
            ("compare.where.pattern", where.get("pattern")),
            *(("compare.mask", m) for m in (mask if isinstance(mask, list) else [])),
        ]
        for label, pattern in patterns:
            if pattern is None:
                continue
            try:
                re.compile(pattern)
            except (re.error, TypeError) as exc:
                errors.append(f"{label} is not a regular expression: {exc}")
    expected = d.get("expected")
    if isinstance(expected, dict):
        for key in ("ensembl_vep_exit", "vep_rs_exit"):
            if key in expected and not isinstance(expected[key], int):
                errors.append(f"expected.{key} is not an integer")
    if (
        not isinstance(expected, dict)
        or "ensembl_vep" not in expected
        or "vep_rs" not in expected
    ):
        errors.append("expected needs ensembl_vep and vep_rs")
    else:
        if as_list(expected["ensembl_vep"]) == as_list(expected["vep_rs"]):
            errors.append(
                "expected.ensembl_vep equals expected.vep_rs: the entry records no divergence"
            )
    if run.get("perl_mode", "offline") not in ("offline", "database"):
        errors.append("run.perl_mode is not offline or database")
    fmt = run.get("format")
    if fmt is not None and fmt not in FORMATS:
        errors.append(f"run.format {fmt!r} is not one of {', '.join(FORMATS)}")
    return errors


def as_list(value) -> list[str]:
    if isinstance(value, list):
        return [str(v) for v in value]
    return [str(value)]


# --- running the engines -------------------------------------------------------


class Paths:
    """The host paths of one entry's run and, for the container, their mounts.

    A placeholder in a flag or argument names a directory each engine reads at its
    own path: ``{corpus}`` the committed corpus, ``{stage}`` the directory the
    ``prepare`` steps write, ``{plugins}`` the plugin directory (the corpus's
    ``plugin_data`` for vep-rs, the Ensembl VEP plugin modules for Ensembl VEP),
    ``{repo}`` the repository root, ``{input}`` the input file and ``{out}`` the
    output file.
    """

    def __init__(
        self,
        repro: Reproduction,
        stage: Path,
        out_dir: Path,
        plugins_dir: Path | None,
        container: bool = True,
        extras: dict[str, Path] | None = None,
    ):
        self.repro = repro
        self.stage = stage
        self.out_dir = out_dir
        self.plugins_dir = plugins_dir
        # False when Ensembl VEP is a local command: it then reads every host path as vep-rs does.
        self.container = container
        # `--extra-dir NAME=PATH`: a directory a reproduction names as `{NAME}`, such as a plugin repository outside the plugin modules.
        self.extras = extras or {}
        self.mounts: dict[str, tuple[Path, bool]] = {}
        json_cache = repro.cache.get("json_cache")
        corpus = repro.cache.get("corpus")
        if corpus:
            self.corpus = self.resolve(corpus)
        elif json_cache:
            self.corpus = self.resolve(json_cache).parent
        else:
            self.corpus = None

    def resolve(self, value: str) -> Path:
        """A host path: `{stage}/...` under the stage, `{repo}/...` or a relative path under the repository."""
        if value.startswith("{stage}"):
            return (self.stage / value[len("{stage}") :].lstrip("/")).resolve()
        if value.startswith("{repo}"):
            return (REPO / value[len("{repo}") :].lstrip("/")).resolve()
        p = Path(value)
        return p.resolve() if p.is_absolute() else (REPO / p).resolve()

    def mount(self, host: Path, container: str, writable: bool = False) -> str:
        self.mounts[container] = (
            host,
            writable or self.mounts.get(container, (host, False))[1],
        )
        return container

    def substitute(self, word: str, engine: str, input_path: Path, output: Path) -> str:
        """A flag word with its placeholders filled for one engine; a container directory is mounted when a word names it."""
        hosts = {
            "{corpus}": (self.corpus, f"{CONTAINER_WORK}/corpus"),
            "{stage}": (self.stage, f"{CONTAINER_WORK}/stage"),
            "{plugins}": (
                self.corpus / "plugin_data" if self.corpus else None,
                f"{CONTAINER_WORK}/plugins",
            ),
            "{repo}": (REPO, f"{CONTAINER_WORK}/repo"),
            "{input}": (input_path, f"{CONTAINER_WORK}/in/{input_path.name}"),
            "{out}": (output, f"{CONTAINER_WORK}/out/{output.name}"),
        }
        for name, host in self.extras.items():
            hosts[f"{{{name}}}"] = (host, f"{CONTAINER_WORK}/extra/{name}")
        for key, (host, container) in hosts.items():
            if key not in word:
                continue
            if key == "{plugins}" and engine != "vep-rs":
                host = self.plugins_dir
                if host is None:
                    raise ReproductionError(f"{key} in {word!r} needs --plugins-dir")
            if host is None:
                raise ReproductionError(
                    f"{key} in {word!r} names no directory: the reproduction has no corpus"
                )
            if engine == "vep-rs" or not self.container:
                word = word.replace(key, str(host))
                continue
            if (
                key in ("{corpus}", "{stage}", "{repo}", "{plugins}")
                or key[1:-1] in self.extras
            ):
                # The stage is the entry's own scratch, so a plugin may write its processed copy of a data file there.
                self.mount(host, container, writable=key == "{stage}")
            word = word.replace(key, container)
        unknown = re.findall(r"\{[a-z][a-z0-9_]*\}", word)
        if unknown:
            raise ReproductionError(
                f"{unknown[0]} in {word!r} names no directory: pass --extra-dir {unknown[0][1:-1]}=<path>"
            )
        return word


def write_input(repro: Reproduction, paths: Paths) -> Path:
    """The input file: the inline text under a VCF header unless it carries its own, or the named file."""
    inp = repro.input
    fmt = inp.get("format", "vcf")
    if "path" in inp:
        return paths.resolve(inp["path"])
    text = inp["text"]
    if not text.endswith("\n"):
        text += "\n"
    if fmt == "vcf" and not text.lstrip().startswith("#"):
        header = ["##fileformat=VCFv4.2", *inp.get("header", [])]
        columns = VCF_HEADER
        samples = inp.get("samples")
        if samples:
            columns += "\tFORMAT\t" + "\t".join(samples)
        text = "\n".join([*header, columns]) + "\n" + text
    name = inp.get("name", f"record.{INPUT_EXTENSIONS[fmt]}")
    path = paths.stage / "input" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def input_chromosomes(repro: Reproduction, input_path: Path) -> list[str]:
    """The chromosomes the input names, for the cache check: the reproduction's own list, else the VCF's first column."""
    listed = repro.cache.get("chromosomes")
    if listed:
        return list(listed)
    if repro.input.get("format", "vcf") != "vcf":
        return []
    chroms = []
    try:
        with open_text(input_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                chrom = line.split("\t")[0].split(" ")[0]
                if chrom not in chroms:
                    chroms.append(chrom)
    except OSError:
        return []
    return chroms


def species_dir(repro: Reproduction) -> str:
    species = repro.cache.get("species", "homo_sapiens")
    source = repro.cache.get("source", "ensembl")
    return species if source == "ensembl" else f"{species}_{source}"


def cache_relative(repro: Reproduction) -> str:
    """`<species>[_merged|_refseq]/<version>_<assembly>`, the cache directory under a root."""
    return f"{species_dir(repro)}/{repro.cache.get('version', '116')}_{repro.cache['assembly']}"


def find_perl_cache(
    repro: Reproduction, roots: list[Path], chroms: list[str]
) -> Path | None:
    """The Ensembl VEP cache directory for the entry's species, source and assembly: the first root holding the chromosomes.

    None when the reproduction prepares its own cache (`cache.perl_cache`) and no root holds one; an error otherwise.
    """
    rel = cache_relative(repro)
    tried = []
    for root in roots:
        candidate = root / rel
        tried.append(str(candidate))
        if candidate.is_dir() and all((candidate / c).is_dir() for c in chroms):
            return candidate
    if repro.cache.get("perl_cache"):
        return None
    need = f" with chromosomes {', '.join(chroms)}" if chroms else ""
    raise ReproductionError(
        f"no Ensembl VEP cache {rel}{need} under the --perl-cache-dir roots (tried {', '.join(tried) or 'none'})"
    )


def prepared_perl_cache(repro: Reproduction, paths: Paths) -> Path:
    """The cache directory under the root `cache.perl_cache` names, once the `prepare` steps have written it."""
    override = repro.cache["perl_cache"]
    relative = cache_relative(repro)
    candidate = paths.resolve(override) / relative
    if not candidate.is_dir():
        raise ReproductionError(f"cache.perl_cache {override} holds no {relative}")
    return candidate


def plain_fasta(fasta: Path, stage: Path) -> Path:
    """Ensembl VEP's FASTA: a `.gz` decompressed once into the stage with its `.fai` beside it, a plain file as it is."""
    if fasta.suffix != ".gz":
        return fasta
    target_dir = stage / "ref"
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / fasta.name[: -len(".gz")]
    if not target.exists():
        with gzip.open(fasta, "rb") as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)
        fai = fasta.with_name(target.name + ".fai")
        if not fai.exists():
            fai = fasta.with_name(fasta.name + ".fai")
        if fai.exists():
            shutil.copy(fai, target.with_name(target.name + ".fai"))
    return target


def run_prepare(
    repro: Reproduction,
    paths: Paths,
    vep: Path,
    perl_cache: Path | None,
    image: str = "",
) -> None:
    """The reproduction's `prepare` shell steps, run in the stage with the run's paths in the environment.

    `IMAGE` is the Ensembl VEP image the run uses, empty without `--image`, for a step that
    must run a tool of the release itself.
    """
    env = dict(os.environ)
    env.update(
        {
            "IMAGE": image,
            "STAGE": str(paths.stage),
            "REPO": str(REPO),
            "VEP": str(vep),
            "VEP_CACHE_BUILDER": str(vep.with_name("vep-cache-builder")),
            "CORPUS": str(paths.corpus) if paths.corpus else "",
            "PERL_CACHE": str(perl_cache) if perl_cache else "",
            "JSON_CACHE": str(paths.resolve(repro.cache["json_cache"]))
            if repro.cache.get("json_cache")
            else "",
        }
    )
    for step in repro.data.get("prepare", []):
        proc = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", step],
            cwd=paths.stage,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise ReproductionError(
                f"prepare step failed ({proc.returncode}): {step}\n{proc.stderr.strip()}"
            )


def engine_args(
    repro: Reproduction,
    engine: str,
    paths: Paths,
    input_path: Path,
    output: Path,
    perl_cache: Path | None,
    fmt: str,
) -> list[str]:
    """The argument list after the program name for one engine's run of a `vep` reproduction."""
    run = repro.run
    cache = repro.cache

    def sub(word: str) -> str:
        return paths.substitute(word, engine, input_path, output)

    args: list[str] = []
    if run.get("stdin"):
        args += ["-i", "STDIN"]
    else:
        args += ["-i", sub("{input}")]
    args += ["-o", sub("{out}")]
    gxf = cache.get("gxf")
    if engine == "vep-rs":
        args.append("--offline")
        if gxf and not cache.get("json_cache"):
            args += [f"--{cache.get('gxf_kind', 'gtf')}", sub(gxf)]
        else:
            args += ["--json_cache", str(paths.resolve(cache["json_cache"]))]
    elif gxf:
        if gxf.startswith("{"):
            file = sub(gxf)
        elif not paths.container:
            file = str(paths.resolve(gxf))
        else:
            host = paths.resolve(gxf)
            file = paths.mount(host.parent, f"{CONTAINER_WORK}/gxf") + "/" + host.name
        args += [f"--{cache.get('gxf_kind', 'gtf')}", file]
    elif run.get("perl_mode", "offline") == "database":
        # Ensembl's public database resolves the input formats its offline parsers refuse (HGVS, SPDI, identifiers).
        args += ["--database", "--db_version", str(cache.get("version", "116"))]
    else:
        assert perl_cache is not None
        version = str(cache.get("version", "116"))
        if paths.container:
            root = f"{CONTAINER_WORK}/cache"
            paths.mount(
                perl_cache,
                f"{root}/{cache_relative(repro)}",
                writable=bool(cache.get("writable")),
            )
        else:
            root = str(perl_cache.parent.parent)
        args += [
            "--offline",
            "--cache",
            "--dir_cache",
            root,
            "--cache_version",
            version,
        ]
    args += ["--species", cache.get("species", "homo_sapiens")]
    if run.get("assembly_flag", True) and "--grch37" not in run.get("flags", []):
        args += ["--assembly", cache["assembly"]]
    in_format = repro.input.get("format", "vcf")
    if in_format != "guess":
        args += ["--format", in_format]
    args += ["--buffer_size", str(run.get("buffer_size", 5000))]
    if engine == "ensembl-vep":
        args += ["--fork", "1"]
    args += ["--force_overwrite", "--quiet"]
    args.append("--stats_text" if run.get("stats") else "--no_stats")
    fasta = cache.get("fasta")
    if fasta:
        host = paths.resolve(fasta)
        if engine == "vep-rs":
            args += ["--fasta", str(host)]
        elif not paths.container:
            args += ["--fasta", str(plain_fasta(host, paths.stage))]
        else:
            plain = plain_fasta(host, paths.stage)
            paths.mount(plain, f"{CONTAINER_WORK}/ref/{plain.name}")
            fai = plain.with_name(plain.name + ".fai")
            if fai.exists():
                paths.mount(fai, f"{CONTAINER_WORK}/ref/{fai.name}")
            args += ["--fasta", f"{CONTAINER_WORK}/ref/{plain.name}"]
    engine_flags = "vep_rs_flags" if engine == "vep-rs" else "ensembl_vep_flags"
    flags = [sub(w) for w in [*run.get("flags", []), *run.get(engine_flags, [])]]
    if engine == "ensembl-vep" and "--plugin" in flags and "--dir_plugins" not in flags:
        if paths.plugins_dir is None:
            raise ReproductionError(
                "a reproduction that loads a plugin needs --plugins-dir, the Ensembl VEP plugin modules"
            )
        flags += ["--dir_plugins", sub("{plugins}")]
    args += flags
    if fmt != "default" and run.get("format_flag", True):
        args.append(f"--{fmt}")
    return args


def tool_args(
    repro: Reproduction, engine: str, paths: Paths, input_path: Path, output: Path
) -> list[str]:
    """The argument list of a `filter` or `recode` run: the reproduction's own words, placeholders filled."""
    return [paths.substitute(w, engine, input_path, output) for w in repro.run["args"]]


def tool_directories(repro: Reproduction, options: argparse.Namespace) -> list[Path]:
    """The `--tool-dir` directories that hold the tools the reproduction's `run.tools` names, one per tool, in order."""
    dirs: list[Path] = []
    for tool in repro.run.get("tools", []):
        holder = next((d for d in options.tool_dirs if (d / tool).is_file()), None)
        if holder is None:
            raise ReproductionError(
                f"the reproduction needs the tool {tool}: pass --tool-dir <directory holding it>"
            )
        if holder not in dirs:
            dirs.append(holder)
    return dirs


def docker_command(
    image: str,
    paths: Paths,
    program: str,
    args: list[str],
    env: dict[str, str],
    stdin: bool,
    name: str | None = None,
) -> list[str]:
    cmd = [
        "docker",
        "run",
        "--rm",
        *(["--name", name] if name else []),
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "-e",
        "HOME=/tmp",
    ]
    if stdin:
        cmd.append("-i")
    for key, val in env.items():
        cmd += ["-e", f"{key}={val}"]
    for container, (host, writable) in sorted(paths.mounts.items()):
        cmd += ["-v", f"{host}:{container}" + ("" if writable else ":ro")]
    cmd += ["-w", f"{CONTAINER_WORK}/out", image, program, *args]
    return cmd


def run_engine(
    engine: str,
    repro: Reproduction,
    paths: Paths,
    options: argparse.Namespace,
    perl_cache: Path | None,
    input_path: Path,
    fmt: str,
) -> EngineRun:
    """One engine's run; its outputs land under `<out_dir>/<engine>/`."""
    out_dir = paths.out_dir / engine
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / FORMATS[fmt]
    tool = repro.run.get("tool", "vep")
    if tool == "vep":
        args = engine_args(repro, engine, paths, input_path, output, perl_cache, fmt)
    else:
        args = tool_args(repro, engine, paths, input_path, output)
    seed = str(repro.run.get("perl_hash_seed", "0"))
    perl_env = {"PERL_HASH_SEED": seed, "PERL_PERTURB_KEYS": "0"} if seed else {}
    stdin = bool(repro.run.get("stdin"))
    stdin_text = input_path.read_bytes() if stdin else None
    run_env = dict(os.environ)
    container = None
    if engine == "vep-rs":
        program = [str(options.vep)] + ([] if tool == "vep" else [tool])
        command = program + args
    elif options.ensembl_vep:
        program = shlex.split(options.ensembl_vep)
        if tool != "vep":
            program[-1] = str(Path(program[-1]).with_name(TOOLS[tool]))
        command = program + args
        run_env.update(perl_env)
        tool_dirs = tool_directories(repro, options)
        if tool_dirs:
            path_entries = [str(d) for d in tool_dirs]
            if run_env.get("PATH"):
                path_entries.append(run_env["PATH"])
            run_env["PATH"] = ":".join(path_entries)
    else:
        paths.mount(out_dir, f"{CONTAINER_WORK}/out", writable=True)
        if not stdin:
            paths.mount(input_path, f"{CONTAINER_WORK}/in/{input_path.name}")
        program = TOOLS[tool]
        tool_dirs = tool_directories(repro, options)
        if tool_dirs:
            # The image lacks the tool a plugin shells out to (LoFTEE's samtools): its directory goes first on PATH.
            mounted = [
                paths.mount(d, f"{CONTAINER_WORK}/tools/{i}")
                for i, d in enumerate(tool_dirs)
            ]
            args = [
                "-c",
                'PATH="$1:$PATH"; shift; exec "$@"',
                "sh",
                ":".join(mounted),
                program,
                *args,
            ]
            program = "sh"
        # A named container can be removed when the client times out; `docker run` alone would leave it running.
        container = f"reevaluate-{repro.name}-{os.getpid()}-{uuid.uuid4().hex[:8]}"
        command = docker_command(
            options.image, paths, program, args, perl_env, stdin, container
        )
    try:
        proc = subprocess.run(
            command,
            input=stdin_text,
            capture_output=True,
            env=run_env,
            cwd=out_dir,
            timeout=options.timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        if container:
            subprocess.run(
                ["docker", "rm", "-f", container], capture_output=True, check=False
            )
        raise
    stdout = proc.stdout.decode("utf-8", "backslashreplace")
    stderr = proc.stderr.decode("utf-8", "backslashreplace")
    (out_dir / "command.txt").write_text(
        shlex.join(command) + f"\nexit {proc.returncode}\n", encoding="utf-8"
    )
    (out_dir / "stdout.txt").write_text(stdout, encoding="utf-8")
    (out_dir / "stderr.txt").write_text(stderr, encoding="utf-8")
    return EngineRun(engine, command, proc.returncode, output, stdout, stderr)


# --- reading the compared value ------------------------------------------------


def is_gzip(path: Path) -> bool:
    """Whether the file starts with the gzip magic, whatever its name."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def open_text(path: Path):
    """A text handle on a plain or gzip file; a byte that is not UTF-8 reads as its `\\xNN` escape, so the two engines' bytes stay distinct."""
    opener = gzip.open if is_gzip(path) else open
    return opener(path, "rt", encoding="utf-8", errors="backslashreplace")


def read_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    with open_text(path) as fh:
        return [line.rstrip("\n") for line in fh]


def table_rows(lines: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    """The default or tab output's column names and rows, the `Extra` column's items folded in as columns."""
    columns: list[str] = []
    rows: list[dict[str, str]] = []
    for line in lines:
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            columns = line[1:].split("\t")
            continue
        if not line.strip() or not columns:
            continue
        cells = line.split("\t")
        row = dict(zip(columns, cells))
        extra = row.get("Extra")
        if extra is not None and extra != "-":
            for item in extra.split(";"):
                key, has, val = item.partition("=")
                row.setdefault(key, val if has else "")
        row["_line"] = line
        rows.append(row)
    return columns, rows


def row_matches(row: dict[str, str], where: dict) -> bool:
    for key, want in where.items():
        if key == "pattern":
            if not re.search(want, row["_line"]):
                return False
            continue
        column = WHERE_COLUMNS.get(key)
        if column is None:
            continue
        have = row.get(column, "")
        wanted = as_list(want)
        if key == "consequence":
            if not any(sorted(have.split(",")) == sorted(w.split(",")) for w in wanted):
                return False
        elif have not in wanted:
            return False
    return True


def cell(row: dict[str, str], name: str) -> str:
    """A column or `Extra` item of a row, `-` when the row lacks it."""
    value = row.get(name)
    return "-" if value is None else value


def section_lines(lines: list[str], section: str) -> list[str]:
    """The lines of one `[section]` of a statistics report, up to the next blank line or section."""
    out: list[str] = []
    inside = False
    for line in lines:
        if line.strip() == section:
            inside = True
            continue
        if inside:
            if not line.strip() or line.startswith("["):
                break
            out.append(line)
    return out


def json_path(value, path: str) -> list:
    """The values at a dotted path: `a.b[].c` iterates a list, `a[key=val].c` keeps the members whose `key` equals `val`."""
    parts = re.findall(r"[^.\[\]]+|\[[^\]]*\]", path)
    current = [value]
    for part in parts:
        nxt = []
        if part == "[]":
            for v in current:
                if isinstance(v, list):
                    nxt.extend(v)
        elif part.startswith("["):
            key, _, want = part[1:-1].partition("=")
            for v in current:
                if isinstance(v, list):
                    nxt.extend(
                        item
                        for item in v
                        if isinstance(item, dict) and str(item.get(key)) == want
                    )
        else:
            for v in current:
                if isinstance(v, dict) and part in v:
                    nxt.append(v[part])
        current = nxt
    return current


def json_value(v) -> str:
    return (
        v
        if isinstance(v, str)
        else json.dumps(v, sort_keys=True, separators=(",", ":"))
    )


def vcf_records(lines: list[str]) -> tuple[dict[str, list[str]], list[dict]]:
    """A VCF output's INFO formats (`CSQ` and the like, from their `##INFO` lines) and its data lines as fields."""
    formats: dict[str, list[str]] = {}
    records = []
    for line in lines:
        if line.startswith("##INFO=<ID="):
            m = re.match(r'##INFO=<ID=([^,]+),.*Format: ([^"]+)"', line)
            if m:
                formats[m.group(1)] = m.group(2).split("|")
            continue
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        info = {}
        if len(fields) > 7:
            for item in fields[7].split(";"):
                key, has, val = item.partition("=")
                info[key] = val if has else ""
        records.append(
            {
                "chrom": fields[0],
                "pos": fields[1] if len(fields) > 1 else "",
                "id": fields[2] if len(fields) > 2 else "",
                "info": info,
                "_line": line,
            }
        )
    return formats, records


def record_matches(rec: dict, where: dict) -> bool:
    """Whether a VCF or JSON record meets every `where` key the record level knows."""
    have = {
        "record": rec.get("id"),
        "id": rec.get("id"),
        "chrom": rec.get("chrom"),
        "pos": rec.get("pos"),
        "location": f"{rec.get('chrom')}:{rec.get('pos')}",
    }
    for key, want in where.items():
        if key == "pattern":
            if not re.search(want, rec["_line"]):
                return False
        elif key in have and have[key] not in as_list(want):
            return False
    return True


def extract(repro: Reproduction, run: EngineRun) -> list[str]:
    """The compared value read from one engine's run, as a list of strings."""
    c = repro.compare
    kind = c["kind"]
    where = c.get("where", {})
    if kind == "exit":
        return [str(run.exit_code)]
    if kind in ("stderr", "stdout"):
        text = run.stderr if kind == "stderr" else run.stdout
        return [line for line in text.splitlines() if re.search(c["pattern"], line)]
    if kind == "warnings":
        return [
            line for line in read_lines(run.warnings) if re.search(c["pattern"], line)
        ]
    if kind == "stats":
        return section_lines(read_lines(run.stats), c["section"])
    lines = read_lines(run.output)
    if kind == "header":
        return [
            line
            for line in lines
            if line.startswith("#") and re.search(c["pattern"], line)
        ]
    if kind == "lines":
        return [line for line in lines if re.search(c["pattern"], line)]
    if kind == "json":
        values = []
        for line in lines:
            if not line.strip():
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            rec = {
                "id": str(obj.get("id", "")),
                "chrom": str(obj.get("seq_region_name", "")),
                "pos": str(obj.get("start", "")),
                "_line": line,
            }
            if not record_matches(rec, where):
                continue
            if "input" in where and obj.get("input") not in as_list(where["input"]):
                continue
            values.extend(json_value(v) for v in json_path(obj, c["path"]))
        return values
    if kind == "vcf":
        formats, records = vcf_records(lines)
        values = []
        for rec in records:
            if not record_matches(rec, where):
                continue
            if c.get("line"):
                values.append(rec["_line"])
                continue
            key = c.get("info", "CSQ")
            field = c["field"]
            if field == key:
                values.append(rec["info"].get(key, "-"))
                continue
            names = formats.get(key)
            if names is None or key not in rec["info"]:
                values.append("-")
                continue
            for entry in rec["info"][key].split(","):
                parts = entry.split("|")
                entry_row = dict(zip(names, parts))
                if "feature" in where and entry_row.get("Feature") not in as_list(
                    where["feature"]
                ):
                    continue
                if "allele" in where and entry_row.get("Allele") not in as_list(
                    where["allele"]
                ):
                    continue
                values.append(entry_row.get(field, "-"))
        return values
    _, rows = table_rows(lines)
    matched = [r for r in rows if row_matches(r, where)]
    if kind == "count":
        return [str(len(matched))]
    if kind == "field":
        field = c["field"]
        if c.get("with_key"):
            return [f"{cell(r, 'Feature')}={cell(r, field)}" for r in matched]
        return [cell(r, field) for r in matched]
    names = c.get("columns", DEFAULT_COLUMNS)
    return [" | ".join(cell(r, n) for n in names) for r in matched]


def normalise(values: list[str], compare: dict) -> list[str]:
    out = []
    for v in values:
        if compare.get("strip", True):
            v = v.rstrip()
        if compare.get("mask_paths", True):
            v = PATH.sub("<path>", v)
        for pattern in compare.get("mask", []):
            v = re.sub(pattern, "<masked>", v)
        out.append(v)
    if compare.get("sorted", True):
        out.sort()
    return out


def verdict(
    expected: dict, compare: dict, ensembl: list[str], vep_rs: list[str]
) -> str:
    want_ens = normalise(as_list(expected["ensembl_vep"]), compare)
    want_rs = normalise(as_list(expected["vep_rs"]), compare)
    if ensembl == vep_rs:
        return "closed"
    if ensembl == want_ens and vep_rs == want_rs:
        return "still-present"
    return "changed"


# --- one entry -----------------------------------------------------------------


def evaluate(repro: Reproduction, options: argparse.Namespace, work: Path) -> Result:
    errors = validate(repro)
    if errors:
        return Result(repro.name, "error", detail="; ".join(errors))
    stage = work / repro.name / "stage"
    out_dir = work / repro.name / "out"
    for d in (stage, out_dir):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    try:
        paths = Paths(
            repro,
            stage,
            out_dir,
            options.plugins_dir,
            container=not options.ensembl_vep,
            extras=getattr(options, "extras", {}),
        )
        tool = repro.run.get("tool", "vep")
        perl_cache = None
        input_path = write_input(repro, paths)
        needs_cache = (
            tool == "vep"
            and not repro.cache.get("gxf")
            and repro.run.get("perl_mode", "offline") == "offline"
        )
        if needs_cache:
            perl_cache = find_perl_cache(
                repro, options.perl_cache_dir, input_chromosomes(repro, input_path)
            )
        run_prepare(repro, paths, options.vep, perl_cache, options.image or "")
        if needs_cache and repro.cache.get("perl_cache"):
            perl_cache = prepared_perl_cache(repro, paths)
        if "path" in repro.input:
            input_path = paths.resolve(repro.input["path"])
        if not input_path.exists():
            raise ReproductionError(f"input {input_path} does not exist")
        fmt = repro.run.get("format") or {"json": "json", "vcf": "vcf"}.get(
            repro.compare["kind"], "default"
        )
        runs = {}
        for engine in ("ensembl-vep", "vep-rs"):
            runs[engine] = run_engine(
                engine, repro, paths, options, perl_cache, input_path, fmt
            )
        for engine, run in runs.items():
            failure = run_failure(repro, engine, run)
            if failure:
                raise ReproductionError(failure)
        values = {
            engine: normalise(extract(repro, run), repro.compare)
            for engine, run in runs.items()
        }
    except (
        ReproductionError,
        OSError,
        subprocess.TimeoutExpired,
        KeyError,
        TypeError,
        ValueError,
        AttributeError,
    ) as exc:
        message = str(exc).strip()
        detail = message.splitlines()[-1] if message else repr(exc)
        return Result(repro.name, "error", detail=detail)
    v = verdict(repro.expected, repro.compare, values["ensembl-vep"], values["vep-rs"])
    detail = ""
    if v != "still-present":
        detail = f"Ensembl VEP {values['ensembl-vep']!r}, vep-rs {values['vep-rs']!r}"
        for engine, run in runs.items():
            if run.exit_code != 0:
                detail += f"; {engine} exit {run.exit_code}: {last_line(run.stderr)}"
    return Result(repro.name, v, values["ensembl-vep"], values["vep-rs"], detail)


def last_line(text: str) -> str:
    text = text.strip()
    return text.splitlines()[-1] if text else ""


# The compare kinds that read the engine's failure itself, so a non-zero exit is their subject.
FAILURE_KINDS = ("exit", "stderr", "stdout")


def run_failure(repro: Reproduction, engine: str, run: EngineRun) -> str | None:
    """Why an engine's run cannot be compared: an exit other than the one recorded for it, or no output to read.

    `expected.<engine>_exit` records an exit status the entry's divergence produces (an engine that refuses
    the input); without it, a kind that reads the output takes a non-zero exit or a missing output file as a
    run that failed, never as the recorded value.
    """
    recorded = repro.expected.get(f"{engine.replace('-', '_')}_exit")
    kind = repro.compare["kind"]
    if recorded is not None:
        if run.exit_code != int(recorded):
            return f"{engine} exited {run.exit_code}, the entry records {recorded}: {last_line(run.stderr)}"
    elif kind in FAILURE_KINDS:
        return None
    elif run.exit_code != 0:
        return f"{engine} exited {run.exit_code}: {last_line(run.stderr)}"
    reads_output = kind not in (*FAILURE_KINDS, "warnings", "stats")
    if reads_output and recorded in (None, 0) and not run.output.exists():
        return f"{engine} wrote no output file"
    return None


# --- the index page --------------------------------------------------------------


def input_classes(repros: list[Reproduction]) -> dict[str, str]:
    """The input class of each reproduction that names one of the three."""
    return {
        r.name: r.data["input_class"]
        for r in repros
        if r.data.get("input_class") in cdo.INPUT_CLASSES
    }


def check(page: Path, reference: str, directory: Path = REPRODUCTIONS) -> list[str]:
    """Every entry without a reproduction, every reproduction without an entry, every schema violation of an entry's reproduction, every version that is not the reference, every `observed` basis `input_class_errors` passes that does not match `SIGHTING` under a tier 3 cell or equal the entry's cell under any other, and every input class the page's tiers and sections contradict."""
    errors = []
    sections = cdo.parse_page(page.read_text(encoding="utf-8"))
    rows = {row.key: row for _, section in cdo.leaves(sections) for row in section.rows}
    repros = {r.name: r for r in load_reproductions(directory=directory)}
    for name, row in rows.items():
        repro = repros.get(name)
        if repro is None:
            errors.append(f"{name}: no reproduction under {directory}")
            continue
        for problem in validate(repro):
            errors.append(f"{name}: {problem}")
        own = repro.data.get("verified_against")
        if own != reference:
            errors.append(
                f"{name}: verified against {own}, the reference is {reference}"
            )
        verified = row.cells[cdo.VERIFIED_CELL]
        if verified != reference:
            errors.append(
                f"{name}: the page says verified against {verified!r}, the reference is {reference}"
            )
        basis = repro.data.get("input_class_basis")
        page_cell = row.cells[cdo.COUNT_CELL]
        observed = repro.data.get("input_class") == "observed"
        # A basis `input_class_errors` refuses, which may be no string, is not held to the page as well.
        if observed and not input_class_errors(repro.data):
            parsed_cell = cdo.CELL.match(page_cell)
            if parsed_cell and parsed_cell.group("tier") == "3":
                # At tier 3 ordinary data does not carry the shape, so the cell is no sighting and the basis names the real data the entry's page records.
                if not SIGHTING.match(basis):
                    errors.append(
                        f"{name}: input_class_basis is {basis!r}: an observed entry of tier 3 reads `{SIGHTING_FORM}`"
                    )
            elif basis != page_cell:
                errors.append(
                    f"{name}: input_class_basis is {basis!r}, the page's `Tier and count` cell is {page_cell!r}"
                )
    for name in repros:
        if name not in rows:
            errors.append(f"{name}: a reproduction with no row on the page")
    errors += cdo.check_classes(sections, input_classes(list(repros.values())))
    return errors


# --- main ------------------------------------------------------------------------


def summary_table(results: list[Result]) -> str:
    counts = {v: sum(1 for r in results if r.verdict == v) for v in VERDICTS}
    lines = ["| verdict | entries |", "| --- | --- |"]
    lines += [f"| {v} | {counts[v]} |" for v in VERDICTS]
    lines.append(f"| total | {len(results)} |")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--image", help="the Ensembl VEP image, a tag or a digest")
    ap.add_argument(
        "--ensembl-vep",
        help="a local Ensembl VEP command to run in place of the image (`vep`, `filter_vep` and `variant_recoder` beside it)",
    )
    ap.add_argument("--vep", type=Path, help="the vep-rs binary")
    ap.add_argument(
        "--perl-cache-dir",
        type=Path,
        action="append",
        default=[],
        help="a root holding <species>[_merged|_refseq]/<version>_<assembly>/ (repeatable, searched in order)",
    )
    ap.add_argument(
        "--plugins-dir",
        type=Path,
        help="the Ensembl VEP plugin modules (the VEP_plugins checkout)",
    )
    ap.add_argument(
        "--tool-dir",
        type=Path,
        action="append",
        default=[],
        dest="tool_dirs",
        help="a directory holding a tool a reproduction's run.tools names (repeatable); put first on Ensembl VEP's PATH",
    )
    ap.add_argument(
        "--extra-dir",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="a directory a reproduction names as {NAME} (repeatable), such as a plugin repository outside the plugin modules",
    )
    ap.add_argument(
        "--entry", action="append", default=[], help="run this entry alone (repeatable)"
    )
    ap.add_argument("--jobs", type=int, default=1, help="entries run at once")
    ap.add_argument(
        "--work-dir",
        type=Path,
        help="where the engines write (default: a temporary directory)",
    )
    ap.add_argument("--keep", action="store_true", help="keep the engines' outputs")
    ap.add_argument(
        "--timeout", type=int, default=1800, help="seconds allowed for one engine run"
    )
    ap.add_argument(
        "--expect-all-present",
        action="store_true",
        help="exit 1 unless every entry is still-present",
    )
    ap.add_argument("--json", type=Path, help="write the results as JSON here")
    ap.add_argument(
        "--check",
        action="store_true",
        help="check the page against the reproductions and run nothing",
    )
    ap.add_argument(
        "--reference",
        default="116.2",
        help="the Ensembl VEP release every entry is verified against (--check)",
    )
    ap.add_argument("--page", type=Path, default=PAGE, help="the index page (--check)")
    ap.add_argument(
        "--reproductions",
        type=Path,
        default=REPRODUCTIONS,
        help="the reproduction directory",
    )
    ap.add_argument(
        "--list", action="store_true", help="list the reproductions and run nothing"
    )
    options = ap.parse_args(argv)

    if options.list:
        try:
            repros = load_reproductions(directory=options.reproductions)
        except ReproductionError as exc:
            print(f"ERROR: [reevaluate] {exc}", file=sys.stderr)
            return 1
        for repro in repros:
            print(
                f"{repro.name}\t{repro.data.get('verified_against', '?')}\t{repro.data.get('page', '?')}"
            )
        return 0
    if options.check:
        try:
            errors = check(options.page, options.reference, options.reproductions)
        except (OSError, ReproductionError, ValueError) as exc:
            print(f"ERROR: [reevaluate] {exc}", file=sys.stderr)
            return 1
        for err in errors:
            print(f"ERROR: [reevaluate] {err}", file=sys.stderr)
        if errors:
            return 1
        classes = input_classes(load_reproductions(directory=options.reproductions))
        print(
            f"{len(classes)} entries have a reproduction verified against {options.reference}: {cdo.class_counts(classes.values())}"
        )
        return 0
    if not options.vep or not (options.image or options.ensembl_vep):
        ap.error("--vep and --image (or --ensembl-vep) are required to run")
    options.extras = {}
    for item in options.extra_dir:
        name, sep, path = item.partition("=")
        if (
            not sep
            or not re.fullmatch(r"[a-z][a-z0-9_]*", name)
            or not Path(path).is_dir()
        ):
            print(
                f"ERROR: [reevaluate] --extra-dir takes NAME=PATH, a lowercase name and an existing directory, not {item!r}",
                file=sys.stderr,
            )
            return 1
        options.extras[name] = Path(path).resolve()
    for d in options.tool_dirs:
        if not d.is_dir():
            print(
                f"ERROR: [reevaluate] --tool-dir {d} is not a directory",
                file=sys.stderr,
            )
            return 1
    options.tool_dirs = [d.resolve() for d in options.tool_dirs]
    if options.image and not options.ensembl_vep and shutil.which("docker") is None:
        print("ERROR: [reevaluate] docker is not on PATH", file=sys.stderr)
        return 1
    try:
        repros = load_reproductions(options.entry or None, options.reproductions)
    except ReproductionError as exc:
        print(f"ERROR: [reevaluate] {exc}", file=sys.stderr)
        return 1
    work = options.work_dir or Path(tempfile.mkdtemp(prefix="reevaluate-"))
    work.mkdir(parents=True, exist_ok=True)
    results: list[Result] = []
    try:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=max(1, options.jobs)
        ) as pool:
            futures = {pool.submit(evaluate, r, options, work): r for r in repros}
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                results.append(result)
                line = f"{result.entry}\t{result.verdict}"
                if result.detail:
                    line += f"\t{result.detail}"
                print(line, flush=True)
    finally:
        # A temporary work directory holds decompressed references of a gigabyte or more: gone unless kept, an interrupted run included.
        if not options.keep and not options.work_dir:
            shutil.rmtree(work, ignore_errors=True)
    results.sort(key=lambda r: r.entry)
    print()
    print(summary_table(results))
    if options.json:
        options.json.write_text(
            json.dumps([dataclasses.asdict(r) for r in results], indent=1) + "\n"
        )
    if options.keep or options.work_dir:
        print(f"outputs under {work}")
    if options.expect_all_present and any(
        r.verdict != "still-present" for r in results
    ):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
