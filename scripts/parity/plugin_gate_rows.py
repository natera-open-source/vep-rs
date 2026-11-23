#!/usr/bin/env python3
"""Derive the plugin concordance gate's rows and the plugins not ported from the support table.

usage: plugin_gate_rows.py [--root DIR] [--rows FILE] [--not-ported FILE] [--record-input FILE]

The support table of docs/plugins.md, read through check_parity_matrix.py's reader, lists every
Ensembl VEP plugin with one of three statuses. A release scores every supported plugin against
Ensembl VEP, and this script says on what:

- a plugin `supported, verified on real data` on the genome-wide ClinVar dataset of each assembly
  its Assemblies cell lists (`clinvar`);
- a plugin `supported on the user's licensed file, verified on synthetic fixtures` on each committed
  corpus of that plugin and assembly (`synthetic_fixture`);
- TranscriptAnnotator, which reads a table of the user's own, on its committed corpora
  (`committed_table`).

An Assemblies cell is read without its parenthesised parts, and the assemblies are the distinct
`GRCh<number>` names left. A corpus is a directory `tests/golden/<set>/<corpus>` and is read through
its manifest.json alone: it belongs to the plugin its `covers` list names, on the manifest's
`assembly`. A corpus whose `covers` names several plugins of the table, or none, is no row. A
second corpus of one plugin and assembly is a row under the fixed label CORPUS_LABELS gives it.
A row scored on ClinVar takes its label from CLINVAR_LABELS where that names its plugin and
assembly, and the plugin's name otherwise.

The plugins the table lists as `not ported` are grouped by reason. The group of each plugin is a
row of scripts/parity/not_ported_groups.tsv (`plugin`, `group`); the groups, their order and their
labels are GROUPS.

Outputs, each optional and byte-identical from one run to the next:

- `--rows`: tab-separated `plugin`, `assembly`, `scored_on`, `dataset_or_corpus`, `label`, in the
  table's order of plugins, then by assembly, a labelled corpus after the plugin's own; the first
  line is the comment `# support_table_sha256: <digest>`, the sha256 of the table's lines as the
  page writes them, header and separator rows included, each ended by one newline;
- `--not-ported`: JSON, `of` (the number of plugins the table lists) and `groups` (each `label`,
  `count` and `plugins`), every group of GROUPS in order;
- `--record-input`: the gate rows and groups as JSON, one object of `support_table_sha256`, `rows`
  (each `plugin`, `assembly`, `scored_on` and `dataset_or_corpus`, where `plugin` carries the
  row's label) and `not_ported`, as `--not-ported` writes it.

Refused, with nothing written: a table check_parity_matrix.py cannot read; a status with no kind
of data here; a supported plugin with no row; a plugin scored on ClinVar on an assembly with no
genome-wide ClinVar dataset; a plugin and assembly scored on corpora with no committed corpus; two
corpora of one plugin and assembly under one label; a manifest that cannot be read; a not-ported
plugin with no group or with two; a group row for a plugin the table does not list as not ported;
a group row that is not two cells or names no group of GROUPS.

Exit 0 with a one-line summary when clean; exit 1 with one `ERROR: [plugin_gate_rows] ...` line
per refusal.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_parity_matrix as parity

COMPONENT = "plugin_gate_rows"
GROUPS_FILE = Path("scripts/parity/not_ported_groups.tsv")
GROUPS_HEADER = ["plugin", "group"]
GOLDEN = Path("tests/golden")
ASSEMBLIES_COLUMN = 3
CLINVAR, FIXTURE, TABLE = "clinvar", "synthetic_fixture", "committed_table"
NOT_PORTED = "not ported"
KIND_OF_STATUS = {
    "supported, verified on real data": CLINVAR,
    "supported on the user's licensed file, verified on synthetic fixtures": FIXTURE,
}
# No published file exists for a table of the user's own, so no ClinVar run can score the plugin.
TABLE_PLUGINS = frozenset({"TranscriptAnnotator"})
CLINVAR_DATASETS = {"GRCh37": "clinvar_grch37", "GRCh38": "clinvar_grch38"}
CORPUS_LABELS = {"GRCh38-plugin-cadd-sv-fixture": "CADD (structural variants)"}
CLINVAR_LABELS = {("LoFTEE", "GRCh37"): "LoFTEE (synthetic conservation data)"}
PARENTHESISED_RE = re.compile(r"\([^()]*\)")
ASSEMBLY_RE = re.compile(r"\bGRCh[0-9]+\b")
ROW_COLUMNS = ("plugin", "assembly", "scored_on", "dataset_or_corpus", "label")

LABEL_RETIRED = (
    "removed from or absent from Ensembl VEP's release 116 collection, "
    "or resting on a retired resource"
)
LABEL_LIVE_SERVICE = "needing a live database, a web service or an external program"
LABEL_GRCH37_ONLY = "data for GRCh37 only"
LABEL_SPECIALISED = "specialised data, rarely used"
LABEL_WRITES_OR_FILTERS = "writing files or reports, or rewriting or filtering output"
LABEL_LICENCE = "commercial or restricted licence"
LABEL_COVERED = "covered by something vep-rs has"
GROUPS = (
    ("retired", LABEL_RETIRED),
    ("live_service", LABEL_LIVE_SERVICE),
    ("grch37_only", LABEL_GRCH37_ONLY),
    ("specialised", LABEL_SPECIALISED),
    ("writes_or_filters", LABEL_WRITES_OR_FILTERS),
    ("licence", LABEL_LICENCE),
    ("covered", LABEL_COVERED),
)


def assemblies_of(cell: str) -> list[str]:
    """The distinct assembly names of an Assemblies cell, its parenthesised parts removed, sorted."""
    return sorted(set(ASSEMBLY_RE.findall(PARENTHESISED_RE.sub("", cell))))


def read_table(root: Path) -> tuple[list[dict], str, list[str]]:
    """(the support table's rows, each with its `assemblies`; the table's sha256; the reasons the
    table cannot be read). The digest covers the header row, the separator row and every plugin row."""
    table, errors = parity.read_support_table(root)
    if errors:
        return [], "", errors
    if not table:
        return [], "", [f"{parity.PLUGINS_PAGE}: the support table has no rows"]
    lines = (root / parity.PLUGINS_PAGE).read_text(encoding="utf-8").splitlines()
    # The header and the separator are the two lines above the first row.
    first, last = table[0]["lineno"] - 2, table[-1]["lineno"]
    digest = hashlib.sha256(
        "".join(f"{line}\n" for line in lines[first - 1 : last]).encode("utf-8")
    ).hexdigest()
    for entry in table:
        cells = parity.table_cells(lines[entry["lineno"] - 1])
        entry["assemblies"] = assemblies_of(cells[ASSEMBLIES_COLUMN])
    return table, digest, []


def read_corpora(root: Path, plugins: set[str]) -> tuple[dict[tuple[str, str], list[str]], list[str]]:
    """(the corpora of each plugin and assembly, by directory name; the manifests that cannot be read).
    A corpus belongs to a plugin when its manifest's own `covers` list names that plugin and no other
    of `plugins`."""
    corpora: dict[tuple[str, str], list[str]] = {}
    errors: list[str] = []
    for manifest_path in sorted((root / GOLDEN).glob("*/*/manifest.json")):
        entry = manifest_path.parent.relative_to(root).as_posix()
        manifest, reason = parity.read_corpus_manifest(root, entry)
        if reason:
            errors.append(reason if reason.startswith(entry) else f"{entry} {reason}")
            continue
        covers = manifest.get("covers")
        covered = {name for name in covers if name in plugins} if isinstance(covers, list) else set()
        if len(covered) != 1:
            continue
        assembly = manifest.get("assembly")
        if not isinstance(assembly, str):
            errors.append(f"{entry}/manifest.json names no assembly")
            continue
        corpora.setdefault((covered.pop(), assembly), []).append(manifest_path.parent.name)
    return corpora, errors


def gate_rows(table: list[dict], corpora: dict[tuple[str, str], list[str]]) -> tuple[list[dict], list[str]]:
    """(the gate's rows in their written order, each a dict of ROW_COLUMNS; the refusals)."""
    rows: list[dict] = []
    errors: list[str] = []
    for entry in table:
        name, status = entry["name"], entry["status"]
        if status == NOT_PORTED:
            continue
        where = f"{parity.PLUGINS_PAGE}:{entry['lineno']}: {name}"
        if status not in KIND_OF_STATUS:
            errors.append(f"{where} has the status {status!r}, which names no kind of data here")
            continue
        kind = TABLE if name in TABLE_PLUGINS else KIND_OF_STATUS[status]
        refused = len(errors)
        own: list[dict] = []
        for assembly in entry["assemblies"]:
            if kind == CLINVAR:
                dataset = CLINVAR_DATASETS.get(assembly)
                if dataset is None:
                    errors.append(f"{where} is scored on ClinVar, and {assembly} has no genome-wide ClinVar dataset")
                    continue
                label = CLINVAR_LABELS.get((name, assembly), name)
                own.append(dict(zip(ROW_COLUMNS, (name, assembly, kind, dataset, label))))
                continue
            found = sorted(corpora.get((name, assembly), []))
            if not found:
                errors.append(f"{where} is scored on corpora, and {GOLDEN} holds none of it on {assembly}")
            labelled: dict[str, str] = {}
            for corpus in found:
                label = CORPUS_LABELS.get(corpus, name)
                if label in labelled:
                    errors.append(
                        f"{where} has the corpora {labelled[label]} and {corpus} on {assembly} under the one "
                        f"label {label!r}; a second corpus needs a label of its own in CORPUS_LABELS"
                    )
                    continue
                labelled[label] = corpus
                own.append(dict(zip(ROW_COLUMNS, (name, assembly, kind, corpus, label))))
        if not own and len(errors) == refused:
            errors.append(f"{where} is supported and has no row: its Assemblies cell names no assembly")
        own.sort(key=lambda row: (row["assembly"], row["label"] != name, row["label"]))
        rows.extend(own)
    return rows, errors


def read_groups(root: Path) -> tuple[dict[str, list[tuple[str, int]]], list[str]]:
    """(per plugin the (group key, line number) of each of its rows in GROUPS_FILE; the reasons the
    file cannot be read)."""
    path = root / GROUPS_FILE
    if not path.is_file():
        return {}, [f"{GROUPS_FILE} is missing"]
    keys = [key for key, _label in GROUPS]
    listed: dict[str, list[tuple[str, int]]] = {}
    errors: list[str] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].split("\t") != GROUPS_HEADER:
        return {}, [f"{GROUPS_FILE}:1: the header is not {GROUPS_HEADER}"]
    for lineno, line in enumerate(lines[1:], start=2):
        cells = line.split("\t")
        if len(cells) != len(GROUPS_HEADER) or not cells[0]:
            errors.append(f"{GROUPS_FILE}:{lineno}: not a plugin and a group")
        elif cells[1] not in keys:
            errors.append(f"{GROUPS_FILE}:{lineno}: {cells[0]} names the group {cells[1]!r}, not one of {', '.join(keys)}")
        else:
            listed.setdefault(cells[0], []).append((cells[1], lineno))
    return listed, errors


def not_ported_block(table: list[dict], listed: dict[str, list[tuple[str, int]]]) -> tuple[dict, list[str]]:
    """(the `not_ported` block: `of` and every group of GROUPS with its plugins in the table's order;
    the refusals)."""
    errors: list[str] = []
    members: dict[str, list[str]] = {key: [] for key, _label in GROUPS}
    not_ported = [entry for entry in table if entry["status"] == NOT_PORTED]
    for entry in not_ported:
        name = entry["name"]
        found = listed.get(name, [])
        if not found:
            errors.append(
                f"{parity.PLUGINS_PAGE}:{entry['lineno']}: {name} is not ported and has no group in {GROUPS_FILE}"
            )
        elif len(found) > 1:
            lines = ", ".join(str(lineno) for _key, lineno in found)
            errors.append(f"{GROUPS_FILE}: {name} has {len(found)} groups (lines {lines}); a plugin has one")
        else:
            members[found[0][0]].append(name)
    names = {entry["name"] for entry in not_ported}
    for name, found in listed.items():
        if name not in names:
            errors.append(
                f"{GROUPS_FILE}:{found[0][1]}: {name} is not a plugin the support table lists as {NOT_PORTED}"
            )
    groups = [
        {"label": label, "count": len(members[key]), "plugins": members[key]} for key, label in GROUPS
    ]
    return {"of": len(table), "groups": groups}, errors


def derive(root: Path) -> tuple[list[dict], dict, str, list[str]]:
    """(the gate's rows, the `not_ported` block, the support table's sha256, every refusal)."""
    table, digest, errors = read_table(root)
    if errors:
        return [], {}, "", errors
    corpora, errors = read_corpora(root, {entry["name"] for entry in table})
    rows, row_errors = gate_rows(table, corpora)
    listed, group_errors = read_groups(root)
    block, block_errors = not_ported_block(table, listed) if not group_errors else ({}, [])
    return rows, block, digest, errors + row_errors + group_errors + block_errors


def rows_text(rows: list[dict], digest: str) -> str:
    lines = [f"# support_table_sha256: {digest}", "\t".join(ROW_COLUMNS)]
    lines.extend("\t".join(row[column] for column in ROW_COLUMNS) for row in rows)
    return "".join(f"{line}\n" for line in lines)


def json_text(value: dict) -> str:
    return json.dumps(value, indent=2, ensure_ascii=False) + "\n"


def record_input(rows: list[dict], block: dict, digest: str) -> dict:
    """The gate rows and groups as one object; a row's `plugin` carries its label."""
    return {
        "support_table_sha256": digest,
        "rows": [
            {
                "plugin": row["label"],
                "assembly": row["assembly"],
                "scored_on": row["scored_on"],
                "dataset_or_corpus": row["dataset_or_corpus"],
            }
            for row in rows
        ],
        "not_ported": block,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--root", default=".", help="repository root (default: the working directory)")
    parser.add_argument("--rows", metavar="FILE", help="write the gate's rows, tab-separated")
    parser.add_argument("--not-ported", metavar="FILE", help="write the groups of plugins not ported, JSON")
    parser.add_argument("--record-input", metavar="FILE", help="write the gate rows and groups as JSON")
    args = parser.parse_args(argv)
    rows, block, digest, errors = derive(Path(args.root))
    for error in errors:
        print(f"ERROR: [{COMPONENT}] {error}", file=sys.stderr)
    if errors:
        return 1
    outputs = (
        (args.rows, rows_text(rows, digest)),
        (args.not_ported, json_text(block)),
        (args.record_input, json_text(record_input(rows, block, digest))),
    )
    for path, text in outputs:
        if path:
            Path(path).write_text(text, encoding="utf-8")
    kinds = ", ".join(
        f"{sum(row['scored_on'] == kind for row in rows)} {kind}" for kind in (CLINVAR, FIXTURE, TABLE)
    )
    grouped = sum(group["count"] for group in block["groups"])
    print(
        f"plugin gate rows: {len(rows)} rows over {len({row['plugin'] for row in rows})} plugins ({kinds}); "
        f"{grouped} of {block['of']} plugins not ported, in {len(block['groups'])} groups"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
