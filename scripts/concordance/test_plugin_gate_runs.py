"""plugin_gate_runs.toml against the rows plugin_gate_rows.py derives, the corpora and the data layout."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts/parity"))
import plugin_gate_rows as pgr
import reproduce_plugin_concordance as driver

RUNS_FILE = HERE / "plugin_gate_runs.toml"
CLASS_TABLE = HERE / "plugin_class_table.tsv"
ORACLE = HERE / "feature_oracle.sh"
LAYOUT = ROOT / "scripts/data/plugin_data_layout.sh"
GOLDEN = ROOT / "tests/golden/116"
WORD_FIELDS = ("vep_rs", "vep_rs_env", "ensembl_vep", "score", "data_indexes")
SCOPES = ("grch37", "grch38", "assembly_agnostic")
LAYOUT_ROLES = ("file", "ancestor", "ancestor_bgzf", "gerp", "synthetic_gerp", "synthetic_conservation", "synthetic_ancestor_bgzf")
INDEX_SUFFIXES = (".tbi", ".csi", ".fai", ".gzi")
# The module a plugin word names, where it is not the support table's name of the plugin.
MODULE_NAMES = {"LoF": "LoFTEE"}
WHOLE_FAMILIES = {
    "plugins_grch37": ("GRCh37", "clinvar_grch37"),
    "plugins_grch38": ("GRCh38", "clinvar_grch38"),
    "gwas": ("GRCh38", "clinvar_grch38"),
    "loftee_grch38": ("GRCh38", "clinvar_grch38"),
    "loftee_grch37": ("GRCh37", "clinvar_grch37"),
}
# The scorer's class of a value on a record the reference dropped. It reads the lines of the
# input as the classes feature_oracle.sh lists under INPUT_READING_MASKS do.
SKIPPED_RECORD_VALUE_CLASS = "plugin_value_on_reference_skipped_record"


def runs() -> list[dict]:
    return tomllib.loads(RUNS_FILE.read_text(encoding="utf-8"))["run"]


def inputs() -> dict:
    return tomllib.loads(RUNS_FILE.read_text(encoding="utf-8"))["inputs"]


def derived() -> list[tuple[str, str, str, str]]:
    rows, _block, _digest, errors = pgr.derive(ROOT)
    assert not errors, errors
    return [(row["label"], row["assembly"], row["scored_on"], row["dataset_or_corpus"]) for row in rows]


def words(run: dict) -> list[str]:
    return [word for field in WORD_FIELDS for word in run.get(field, [])] + [run["reference"]]


def plugin_words(run: dict, field: str) -> list[str]:
    flags = run[field]
    return [flags[i + 1] for i, word in enumerate(flags) if word == "--plugin"]


def input_reading_classes() -> set[str]:
    """The classes a comparison decides from the lines of the run's input: those on the one
    line of feature_oracle.sh that lists them, which that runner's own tests hold to the
    comparator's lists, and the scorer's class of a value on a record the reference dropped."""
    (listed,) = re.findall(r'^INPUT_READING_MASKS="([^"\n]*)"$', ORACLE.read_text(encoding="utf-8"), flags=re.MULTILINE)
    return set(listed.split()) | {SKIPPED_RECORD_VALUE_CLASS}


def declared_classes(plugin: str) -> set[str]:
    """The classes the class table declares for a plugin."""
    lines = CLASS_TABLE.read_text(encoding="utf-8").splitlines()
    rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
    return {row[0] for row in rows if plugin in row[1].split(";")}


def layout(function: str, *arguments: str) -> list[str]:
    """The lines a function of the layout script prints, or none where it refuses."""
    done = subprocess.run(
        ["bash", "-c", 'source "$1"; shift; "$@"', "bash", str(LAYOUT), function, *arguments],
        capture_output=True, text=True, check=False,
    )
    return done.stdout.split() if done.returncode == 0 else []


def test_the_file_holds_the_five_whole_runs_and_every_set_once() -> None:
    keys = [driver.run_key(run) for run in runs()]
    assert len(keys) == len(set(keys))
    whole = {run["family"]: (run["assembly"], run["dataset"]) for run in runs() if "family" in run}
    assert whole == WHOLE_FAMILIES
    for run in runs():
        assert ("family" in run) != ("corpus" in run and "set" in run), driver.run_key(run)


def test_every_derived_row_stands_on_exactly_one_whole_run_or_on_the_sets_of_its_own_corpus() -> None:
    rows = derived()
    assert rows
    for row in rows:
        holders = [run for run in runs() if list(row) in run["rows"]]
        _label, _assembly, scored_on, target = row
        if scored_on == pgr.CLINVAR:
            assert [("family" in run, run.get("dataset")) for run in holders] == [(True, target)], row
        else:
            assert holders and all(run.get("corpus") == target for run in holders), row
            assert sum(run["set"] == driver.BASE_SET for run in holders) == 1, row


def test_no_run_names_a_row_the_rows_script_does_not_derive() -> None:
    rows = {tuple(row) for row in derived()}
    for run in runs():
        assert run["rows"], driver.run_key(run)
        for row in run["rows"]:
            assert tuple(row) in rows, (driver.run_key(run), row)
            assert row[1] == run["assembly"]


def test_a_whole_run_scores_the_plugins_its_score_words_name() -> None:
    for run in runs():
        named = run["score"][run["score"].index("--plugins") + 1].split(",")
        assert sorted(named) == sorted({driver.plugin_of(row[0]) for row in run["rows"]}), driver.run_key(run)


def test_a_whole_run_that_loads_revel_hands_its_table_to_the_scorer() -> None:
    """A class of REVEL reads the table the run's plugin read and sets nothing aside without
    it, so a whole run whose plugins hold REVEL names that file after ``--revel-table`` among
    its score words, and ``--revel-no-match`` where the plugin's words set ``no_match``."""
    seen = 0
    for run in runs():
        if "family" not in run:
            continue
        revel = [word.split(",") for word in plugin_words(run, "vep_rs") if word.split(",")[0] == "REVEL"]
        if not revel:
            assert "--revel-table" not in run["score"], driver.run_key(run)
            continue
        seen += 1
        (words,) = revel
        file = next(word.removeprefix("file=") for word in words[1:] if word.startswith("file=") or "=" not in word)
        assert run["score"][run["score"].index("--revel-table") + 1] == file, driver.run_key(run)
        unmatched = any(word.startswith("no_match=") and word not in ("no_match=", "no_match=0") for word in words)
        assert ("--revel-no-match" in run["score"]) == unmatched, driver.run_key(run)
    assert seen == 2


def test_a_whole_run_that_scores_a_class_that_reads_the_input_hands_it_to_the_scorer() -> None:
    """The records of a whole run's input stand under identifiers, and a row under an
    identifier states neither the alleles of its record nor the reference bases its allele
    stands in for. A class that reads either is handed them from the lines of the input and
    sets aside no row that needs them without those lines, so a whole run that scores a
    plugin the class table declares such a class for names ``{input}`` after ``--input``
    among its score words."""
    reading = input_reading_classes()
    assert SKIPPED_RECORD_VALUE_CLASS in declared_classes("UTRAnnotator")
    held = []
    for run in runs():
        if "family" not in run:
            continue
        named = run["score"][run["score"].index("--plugins") + 1].split(",")
        read = sorted({klass for plugin in named for klass in declared_classes(plugin)} & reading)
        if not read:
            continue
        key = driver.run_key(run)
        assert run["score"].count("--input") == 1, (key, read)
        assert run["score"][run["score"].index("--input") + 1] == "{input}", key
        held.append(key)
    assert held == ["plugins_grch37", "plugins_grch38", "loftee_grch38", "loftee_grch37"]


def test_the_plugin_arguments_of_the_two_engines_are_equal() -> None:
    for run in runs():
        ours, theirs = plugin_words(run, "vep_rs"), plugin_words(run, "ensembl_vep")
        assert ours == theirs, driver.run_key(run)
        if "family" in run or run["set"] != driver.BASE_SET:
            assert ours, driver.run_key(run)


def test_no_word_names_a_place_outside_the_tree_and_every_named_input_is_listed() -> None:
    listed = set(inputs())
    assert listed == set(driver.INPUT_OPTIONS) | set(driver.DERIVED_INPUTS)
    for run in runs():
        for word in words(run):
            assert not re.search(r"(^|[=,:])[/~]", word), (driver.run_key(run), word)
            assert "<" not in word and ">" not in word, (driver.run_key(run), word)
            assert set(driver.PLACEHOLDER_RE.findall(word)) <= listed, (driver.run_key(run), word)
            assert "{" not in driver.PLACEHOLDER_RE.sub("", word), (driver.run_key(run), word)


def test_a_plugin_data_path_is_the_layout_scripts_for_its_plugin_and_assembly() -> None:
    checked = 0
    for run in runs():
        if "family" not in run:
            assert not any("{plugin_data}" in word for word in words(run)), driver.run_key(run)
            continue
        assembly = run["assembly"].lower()
        indexes = set(run.get("data_indexes", []))
        for field in ("vep_rs", "ensembl_vep"):
            for word in plugin_words(run, field):
                name = MODULE_NAMES.get(word.split(",")[0], word.split(",")[0])
                paths = driver.PLUGIN_DATA_PATH_RE.findall(word)
                held = {path for role in LAYOUT_ROLES for path in layout("plugin_data_file", name, role, assembly)}
                staged = layout("plugin_staged_files", name, assembly)
                for path in paths:
                    scope, _, relative = path.removeprefix("{plugin_data}/").partition("/")
                    assert scope in (assembly, "assembly_agnostic"), path
                    assert relative in held, (name, path)
                    wanted = {f"{path}{suffix}" for suffix in INDEX_SUFFIXES if f"{relative}{suffix}" in staged}
                    assert wanted <= indexes, (name, path)
                    assert {index for index in indexes if index.startswith(f"{path}.")} == wanted, (name, path)
                    checked += 1
    assert checked >= 26


def test_every_index_stands_beside_a_data_file_the_run_names() -> None:
    for run in runs():
        named = {path for word in run["vep_rs"] + run["ensembl_vep"] for path in driver.PLUGIN_DATA_PATH_RE.findall(word)}
        for index in run.get("data_indexes", []):
            assert index.endswith(INDEX_SUFFIXES), index
            assert index.rsplit(".", 1)[0] in named, (driver.run_key(run), index)


def test_every_corpus_set_exists_with_its_expected_folder() -> None:
    for run in runs():
        if "family" in run:
            assert run["reference"] == "{reference}"
            continue
        corpus = GOLDEN / run["corpus"]
        manifest = json.loads((corpus / "manifest.json").read_text(encoding="utf-8"))
        assert manifest["assembly"] == run["assembly"]
        if run["set"] == driver.BASE_SET:
            folder = "expected"
        else:
            (entry,) = [entry for entry in manifest["flag_sets"] if entry["name"] == run["set"]]
            folder = entry.get("expected_dir") or f"expected_{run['set']}"
        assert run["reference"] == f"{{corpus}}/{folder}/default.txt.gz", driver.run_key(run)
        assert (corpus / folder / "default.txt.gz").is_file(), driver.run_key(run)


def test_every_set_of_a_gate_corpus_has_a_run() -> None:
    for corpus in {run["corpus"] for run in runs() if "corpus" in run}:
        manifest = json.loads((GOLDEN / corpus / "manifest.json").read_text(encoding="utf-8"))
        own = {entry["name"] for entry in manifest.get("flag_sets", [])} | {driver.BASE_SET}
        assert {run["set"] for run in runs() if run.get("corpus") == corpus} == own, corpus


def test_a_corpus_run_carries_the_words_its_committed_provenance_records() -> None:
    names = {"<out>": "{ensembl_vep_output}", "<cache>": "{ensembl_vep_cache}", "{plugins}": "{plugin_modules}"}
    exact = 0
    for run in runs():
        if "family" in run:
            continue
        provenance = json.loads((GOLDEN / run["corpus"] / "provenance.json").read_text(encoding="utf-8"))
        assert run["ensembl_vep_image"] == provenance["vep_image"], driver.run_key(run)
        block = provenance if run["set"] == driver.BASE_SET else provenance["flag_sets"][run["set"]]
        command = block.get("vep_command")
        if command is None or "[," in command:
            # The set records its flags alone, or one command stands for every set: its plugin words are held.
            recorded = block.get("flags") or run["ensembl_vep"]
            assert plugin_words(run, "ensembl_vep") == plugin_words({"flags": recorded}, "flags"), driver.run_key(run)
            continue
        recorded = command.split()[1:-1]
        for old, new in names.items():
            recorded = [word.replace(old, new) for word in recorded]
        assert run["ensembl_vep"] == recorded, driver.run_key(run)
        exact += 1
    assert exact >= 60


def spliceai_table_words(run: dict) -> list[str]:
    """The score words that hand the scorer the files a run's SpliceAI plugin words name after
    ``snv=`` and ``indel=``, in that order; none for a run that does not load the plugin."""
    words = [word.split(",") for word in plugin_words(run, "vep_rs") if word.split(",")[0] == "SpliceAI"]
    if not words:
        return []
    (named,) = words
    files = {key: value for key, _, value in (word.partition("=") for word in named[1:])}
    return ["--spliceai-table", files["snv"], "--spliceai-table", files["indel"]]


def test_a_run_that_loads_spliceai_hands_its_score_files_to_the_scorer() -> None:
    """The class ``spliceai_multi_entry_cell`` reads the score files the run's plugin read and
    sets nothing aside without them, so every run whose plugins hold SpliceAI names both files
    after ``--spliceai-table`` among its score words, and no other run names the option."""
    seen = 0
    for run in runs():
        words = spliceai_table_words(run)
        if not words:
            assert "--spliceai-table" not in run["score"], driver.run_key(run)
            continue
        seen += 1
        assert run["score"][-len(words):] == words and run["score"].count("--spliceai-table") == 2, driver.run_key(run)
        assert all(word.startswith("{corpus}/plugin_data/") for word in words[1::2]), driver.run_key(run)
    assert seen == 8


def test_a_corpus_run_of_vep_rs_reads_the_corpus_and_runs_inside_it() -> None:
    for run in runs():
        if "family" in run:
            assert "vep_rs_env" not in run
            assert run["vep_rs"][:4] == ["-i", "{input}", "-o", "{output}"]
            assert run["vep_rs"][run["vep_rs"].index("--fork") + 1] == "{fork}"
            continue
        assert run["vep_rs_env"] == ["HOME={home}", "VEP_DIR_PLUGINS={corpus}"], driver.run_key(run)
        assert run["vep_rs"][:2] == ["-i", "{corpus}/variants.vcf"] and run["vep_rs"][-2:] == ["-o", "{output}"]
        assert run["score"] == ["--plugins", driver.plugin_of(run["rows"][0][0]), "--json-cache", "{corpus}/json_cache", *spliceai_table_words(run)]


def test_every_run_names_the_image_the_generator_pins() -> None:
    pinned = driver.pinned_image()
    assert re.fullmatch(r"ensemblorg/ensembl-vep@sha256:[0-9a-f]{64}", pinned)
    assert {run["ensembl_vep_image"] for run in runs()} == {pinned}
    assert [driver.run_key(run) for run in runs() if run.get("ensembl_vep_tools")] == ["loftee_grch38", "loftee_grch37"]
