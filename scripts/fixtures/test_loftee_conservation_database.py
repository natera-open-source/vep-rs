"""Tests for the synthetic LoFTEE conservation database generator on hand-built caches, and
of each committed corpus file against the rows its recorded seed gives."""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import random
import sqlite3
import sys
from collections import Counter
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import loftee_conservation_database as lcd

GOLDEN = HERE.parents[1] / "tests" / "golden" / "116"
# Each committed file: its corpus, its name, and the branch of the plugin whose layout it has.
COMMITTED = (
    ("GRCh37-plugin-loftee-fixture", "phylocsf_gerp.synthetic.sql", "master"),
    ("GRCh38-plugin-loftee-fixture", "loftee.synthetic.sql", "grch38"),
)
# The two statements of the plugin's `master` branch (LoF.pm:570, gerp_dist.pl:129).
PHYLOCSF_QUERY = "SELECT * FROM phylocsf_data WHERE transcript = ? AND exon_number = ?;"
GERP_QUERY = "SELECT * FROM gerp_exons where transcript_id = ? AND exon_num = ?; "
# The one statement of its `grch38` branch (LoF.pm:556 there).
SUMMARY_QUERY = "SELECT * FROM phylocsf_summary WHERE transcript = ? AND exon = ?;"


def exons(count: int) -> list[dict]:
    return [{"start": 100 * n, "end": 100 * n + 50} for n in range(1, count + 1)]


def transcript(stable_id: str, count: int, biotype: str = "protein_coding") -> dict:
    return {"stable_id": stable_id, "biotype": biotype, "exons": exons(count)}


CACHE_SHARDS = {
    "21/1-1000000.json": [
        transcript("ENST_B", 3),
        transcript("ENST_A", 5),
        None,
        transcript("ENST_DECAY", 4, "nonsense_mediated_decay"),
        {"stable_id": "ENST_NO_EXONS", "biotype": "protein_coding"},
    ],
    "21/1000001-2000000.json.gz": [transcript("ENST_A", 5), transcript("ENST_INTRON", 2, "retained_intron")],
    "22/1-1000000.json": [transcript("ENST_C", 1)],
}


def write_cache(corpus: Path, shards: dict[str, list]) -> Path:
    """A `json_cache` under `corpus` holding the named shards, a `.gz` one compressed."""
    cache = corpus / "json_cache"
    for name, items in shards.items():
        path = cache / "transcripts" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith(".gz"):
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(json.dumps(items))
        else:
            path.write_text(json.dumps(items), encoding="utf-8")
    return cache


def printed(value: float) -> str:
    """A number as Perl prints it: up to fifteen significant digits."""
    return f"{value:.15g}"


def test_values_derive_from_the_seed_and_the_row_alone():
    """A value is a function of the seed, the table, the transcript and the exon number;
    the outcomes of the PhyloCSF rows stand in their shares, a score has at most three
    decimals in Perl's print, the highest score is never below the annotated one, exact
    zeros stand among both and none is negative; a GERP total has two decimals inside the
    bounds; another seed moves the values."""
    keys = [(f"ENST_T{t}", n) for t in range(300) for n in range(1, 11)]
    scores = {key: lcd.phylocsf_scores(7, *key) for key in keys}
    assert scores == {key: lcd.phylocsf_scores(7, *key) for key in keys}
    outcomes = Counter(lcd.phylocsf_outcome(value) for value in scores.values())
    assert set(outcomes) == set(lcd.OUTCOMES)
    assert abs(outcomes[lcd.NO_ROW] / len(keys) - 3 / 18) < 0.03
    for name in (lcd.NO_FLAG, lcd.UNLIKELY_ORF, lcd.WEAK):
        assert abs(outcomes[name] / len(keys) - 5 / 18) < 0.03, name
    rows = [value for value in scores.values() if value is not None]
    for annotated, highest in rows:
        assert highest >= annotated
        for score in (annotated, highest):
            text = printed(score)
            assert float(text) == score and len(text.partition(".")[2]) <= 3 and "e" not in text, text
            assert math.copysign(1.0, score) == 1.0 or score < 0, "no negative zero"
    assert {len(printed(annotated).partition(".")[2]) for annotated, _ in rows} == {0, 1, 2, 3}
    assert any(annotated == 0 for annotated, _ in rows) and any(highest == 0 for _, highest in rows)
    assert all(highest > 0 for annotated, highest in rows if lcd.phylocsf_outcome((annotated, highest)) == lcd.UNLIKELY_ORF)
    assert all(annotated < 0 and highest <= 0 for annotated, highest in rows if lcd.phylocsf_outcome((annotated, highest)) == lcd.WEAK)
    assert lcd.phylocsf_outcome(None) == lcd.NO_ROW
    assert [lcd.phylocsf_outcome(pair) for pair in ((0.0, 0.0), (1.5, 2.0), (-0.001, 0.001), (-2.0, 0.0), (-2.0, -0.5))] == [
        lcd.NO_FLAG, lcd.NO_FLAG, lcd.UNLIKELY_ORF, lcd.WEAK, lcd.WEAK,
    ]

    totals = {key: lcd.gerp_exon_total(7, *key) for key in keys}
    present = [value for value in totals.values() if value is not None]
    assert abs((len(keys) - len(present)) / len(keys) - 1 / 8) < 0.03
    low, high = (bound / 100 for bound in lcd.GERP_HUNDREDTHS)
    assert (low, high) == (-300.0, 900.0) and all(low <= value <= high for value in present)
    assert all(float(printed(value)) == value and len(printed(value).partition(".")[2]) <= 2 for value in present)
    key = next(key for key, value in totals.items() if value is not None)
    rng = random.Random(f"7:gerp_exons:{key[0]}:{key[1]}")
    assert rng.randrange(8) != 0 and totals[key] == rng.randint(-30_000, 90_000) / 100

    other_scores = {key: lcd.phylocsf_scores(8, *key) for key in keys}
    assert sum(scores[key] != other_scores[key] for key in keys) > len(keys) * 0.8

    summary = {key: lcd.phylocsf_scores(7, *key, table=lcd.SUMMARY_TABLE) for key in keys}
    assert summary == {key: lcd.phylocsf_scores(7, *key, table="phylocsf_summary") for key in keys}
    assert sum(summary[key] != scores[key] for key in keys) > len(keys) * 0.8, "each table draws its own values"
    summary_outcomes = Counter(lcd.phylocsf_outcome(value) for value in summary.values())
    assert abs(summary_outcomes[lcd.NO_ROW] / len(keys) - 3 / 18) < 0.03
    for name in (lcd.NO_FLAG, lcd.UNLIKELY_ORF, lcd.WEAK):
        assert abs(summary_outcomes[name] / len(keys) - 5 / 18) < 0.03, name
    assert sum(totals[key] != lcd.gerp_exon_total(8, *key) for key in keys) > len(keys) * 0.8


def test_rows_cover_the_exons_of_the_protein_coding_transcripts(tmp_path: Path):
    """The transcripts are the cache's protein-coding ones that have an exon, each once
    whatever the shards that list it, a compressed shard read as a plain one; the rows of
    both tables stand in transcript and exon order, one per exon that draws a row; a row
    is the same over a cache that holds more transcripts; a transcript listed with two
    exon counts is refused."""
    cache = write_cache(tmp_path / "corpus", CACHE_SHARDS)
    counts = lcd.read_coding_exon_counts(cache)
    assert counts == {"ENST_A": 5, "ENST_B": 3, "ENST_C": 1} and list(counts) == ["ENST_A", "ENST_B", "ENST_C"]
    phylocsf, gerp, recorded = lcd.database_rows(counts, 7)
    keys = [("ENST_A", n) for n in range(1, 6)] + [("ENST_B", n) for n in range(1, 4)] + [("ENST_C", 1)]
    assert phylocsf == [(*key, *scores) for key in keys if (scores := lcd.phylocsf_scores(7, *key)) is not None]
    assert gerp == [(*key, total) for key in keys if (total := lcd.gerp_exon_total(7, *key)) is not None]
    assert recorded == {"transcripts": 3, "exons": 9, "phylocsf_data": len(phylocsf), "gerp_exons": len(gerp)}
    assert 0 < len(phylocsf) <= 9 and 0 < len(gerp) <= 9

    wider_phylocsf, wider_gerp, wider_recorded = lcd.database_rows({**counts, "ENST_D": 40}, 7)
    assert set(phylocsf) < set(wider_phylocsf) and set(gerp) < set(wider_gerp)
    assert wider_recorded["transcripts"] == 4 and wider_recorded["exons"] == 49
    assert len(wider_phylocsf) < 49 and len(wider_gerp) < 49, "some exons draw no row"

    clash = write_cache(tmp_path / "clash", {**CACHE_SHARDS, "22/1000001-2000000.json": [transcript("ENST_C", 2)]})
    with pytest.raises(SystemExit, match="ENST_C is listed with 1 exons and with 2"):
        lcd.read_coding_exon_counts(clash)


def test_main_writes_the_two_tables_and_records_the_file(tmp_path: Path, capsys):
    """The tool writes the two tables under the named directory with the declared column
    types; the plugin's statements, their values bound as text, return the row of an exon
    that has one and none otherwise; `read_rows` reads the generator's rows back;
    provenance records the generator, the seed, the counts, the digest and the library's
    version and keeps its other keys; an existing file is kept without --force; a name
    without `synthetic`, a path as a name or a directory, a corpus without a cache and a
    cache without a protein-coding exon are refused."""
    corpus = tmp_path / "corpus"
    corpus.mkdir()
    name = "phylocsf_gerp.synthetic.sql"
    args = ["--corpus", str(corpus), "--name", name, "--seed", "7", "--dir", "plugin_data"]
    with pytest.raises(SystemExit, match="corpus has no json_cache/transcripts"):
        lcd.main(args)
    write_cache(corpus, CACHE_SHARDS)
    kept = {"vep_image": "an image", "plugin_fixtures": {"other.synthetic.txt.gz": {"seed": 3}}}
    (corpus / "provenance.json").write_text(json.dumps(kept), encoding="utf-8")
    phylocsf, gerp, counts = lcd.database_rows(lcd.read_coding_exon_counts(corpus / "json_cache"), 7)

    assert lcd.main(args) == 0
    assert capsys.readouterr().out.strip() == (
        f"wrote plugin_data/{name}: phylocsf_data {len(phylocsf)} rows, gerp_exons {len(gerp)} rows, over 9 exons of 3 transcripts, seed 7"
    )
    out = corpus / "plugin_data" / name
    assert sorted(path.name for path in out.parent.iterdir()) == [name]
    assert lcd.read_rows(out) == (phylocsf, gerp)

    connection = sqlite3.connect(out)
    try:
        declared = {
            table: [(column[1], column[2]) for column in connection.execute(f"PRAGMA table_info({table})")]
            for table in ("gerp_exons", "phylocsf_data")
        }
        assert declared == {
            "gerp_exons": [("transcript_id", "TEXT"), ("exon_num", "INTEGER"), ("gerp", "REAL")],
            "phylocsf_data": [("transcript", "TEXT"), ("exon_number", "INTEGER"), ("corresponding_orf_score", "REAL"), ("max_score", "REAL")],
        }
        assert [row[0] for row in connection.execute("SELECT name FROM sqlite_master ORDER BY name")] == ["gerp_exons", "phylocsf_data"]
        for key in [("ENST_A", n) for n in range(1, 6)] + [("ENST_B", n) for n in range(1, 4)] + [("ENST_C", 1)]:
            bound = (key[0], str(key[1]))
            scores = lcd.phylocsf_scores(7, *key)
            assert connection.execute(PHYLOCSF_QUERY, bound).fetchone() == (None if scores is None else (*key, *scores))
            total = lcd.gerp_exon_total(7, *key)
            assert connection.execute(GERP_QUERY, bound).fetchone() == (None if total is None else (*key, total))
        for bound in (("ENST_A", "6"), ("ENST_A", "1-2"), ("ENST_DECAY", "1"), ("ENST_NO_EXONS", "1")):
            assert connection.execute(PHYLOCSF_QUERY, bound).fetchone() is None
            assert connection.execute(GERP_QUERY, bound).fetchone() is None
        stored = {type(value) for row in connection.execute("SELECT gerp FROM gerp_exons") for value in row}
        assert stored == {float}
    finally:
        connection.close()

    record = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))
    assert record["vep_image"] == "an image" and record["plugin_fixtures"]["other.synthetic.txt.gz"] == {"seed": 3}
    assert record["plugin_fixtures"][name] == {
        "generator": "scripts/fixtures/loftee_conservation_database.py",
        "branch": "master",
        "seed": 7,
        "rows": counts,
        "files": {f"plugin_data/{name}": hashlib.sha256(out.read_bytes()).hexdigest()},
        "tools": {"sqlite": sqlite3.sqlite_version},
    }

    with pytest.raises(SystemExit, match="exists; pass --force to overwrite it"):
        lcd.main(args)
    assert lcd.main([*args[:5], "8", *args[6:], "--force"]) == 0
    capsys.readouterr()
    assert lcd.read_rows(out) == lcd.database_rows(lcd.read_coding_exon_counts(corpus / "json_cache"), 8)[:2] != (phylocsf, gerp)
    assert json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]["seed"] == 8

    with pytest.raises(SystemExit, match="--name carries `synthetic`"):
        lcd.main(["--corpus", str(corpus), "--name", "phylocsf_gerp.sql", "--seed", "7"])
    with pytest.raises(SystemExit, match="--name is a file name, not a path"):
        lcd.main(["--corpus", str(corpus), "--name", "plugin_data/x.synthetic.sql", "--seed", "7"])
    with pytest.raises(SystemExit, match="--dir is a directory name inside the corpus"):
        lcd.main(["--corpus", str(corpus), "--name", name, "--seed", "7", "--dir", "a/b"])
    empty = tmp_path / "empty"
    write_cache(empty, {"21/1-1000000.json": [transcript("ENST_DECAY", 4, "nonsense_mediated_decay"), {"stable_id": "ENST_NO_EXONS", "biotype": "protein_coding"}]})
    with pytest.raises(SystemExit, match="no protein-coding transcript with an exon"):
        lcd.main(["--corpus", str(empty), "--name", name, "--seed", "7"])


def test_the_grch38_layout_adds_the_table_that_branch_queries(tmp_path: Path, capsys):
    """Under `--branch grch38` the file holds three tables: `phylocsf_summary`, which that
    branch's statement reads by `exon`, its values bound as text, and the two tables of
    `master` with the rows a `master` file holds, which that branch never reads; the rows
    of `phylocsf_summary` are drawn apart from those of `phylocsf_data`, so a run that read
    the other table would write other items; provenance records the branch and the third
    table's count; a `master` file holds no `phylocsf_summary`."""
    corpus = tmp_path / "corpus"
    write_cache(corpus, CACHE_SHARDS)
    name = "loftee.synthetic.sql"
    args = ["--corpus", str(corpus), "--name", name, "--seed", "7", "--dir", "plugin_data"]
    exon_counts = lcd.read_coding_exon_counts(corpus / "json_cache")
    phylocsf, gerp, counts = lcd.database_rows(exon_counts, 7)
    summary = lcd.summary_rows(exon_counts, 7)
    keys = [("ENST_A", n) for n in range(1, 6)] + [("ENST_B", n) for n in range(1, 4)] + [("ENST_C", 1)]
    assert summary == [(*key, *scores) for key in keys if (scores := lcd.phylocsf_scores(7, *key, table="phylocsf_summary")) is not None]
    assert 0 < len(summary) <= 9 and summary != phylocsf

    assert lcd.main([*args, "--branch", "grch38"]) == 0
    assert capsys.readouterr().out.strip() == (
        f"wrote plugin_data/{name}: phylocsf_summary {len(summary)} rows, phylocsf_data {len(phylocsf)} rows, "
        f"gerp_exons {len(gerp)} rows, over 9 exons of 3 transcripts, seed 7"
    )
    out = corpus / "plugin_data" / name
    assert lcd.read_summary_rows(out) == summary
    assert lcd.read_rows(out) == (phylocsf, gerp)

    connection = sqlite3.connect(out)
    try:
        assert [row[0] for row in connection.execute("SELECT name FROM sqlite_master ORDER BY name")] == [
            "gerp_exons", "phylocsf_data", "phylocsf_summary",
        ]
        assert [(column[1], column[2]) for column in connection.execute("PRAGMA table_info(phylocsf_summary)")] == [
            ("transcript", "TEXT"), ("exon", "INTEGER"), ("corresponding_orf_score", "REAL"), ("max_score", "REAL"),
        ]
        for key in keys:
            scores = lcd.phylocsf_scores(7, *key, table="phylocsf_summary")
            assert connection.execute(SUMMARY_QUERY, (key[0], str(key[1]))).fetchone() == (None if scores is None else (*key, *scores))
        for bound in (("ENST_A", "6"), ("ENST_A", "1-2"), ("ENST_DECAY", "1"), ("ENST_NO_EXONS", "1")):
            assert connection.execute(SUMMARY_QUERY, bound).fetchone() is None
    finally:
        connection.close()

    entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    assert entry == {
        "generator": "scripts/fixtures/loftee_conservation_database.py",
        "branch": "grch38",
        "seed": 7,
        "rows": {**counts, "phylocsf_summary": len(summary)},
        "files": {f"plugin_data/{name}": hashlib.sha256(out.read_bytes()).hexdigest()},
        "tools": {"sqlite": sqlite3.sqlite_version},
    }

    assert lcd.main([*args, "--force"]) == 0
    capsys.readouterr()
    assert lcd.read_rows(out) == (phylocsf, gerp)
    with pytest.raises(sqlite3.OperationalError, match="no such table: phylocsf_summary"):
        lcd.read_summary_rows(out)
    assert json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]["branch"] == "master"


@pytest.mark.parametrize(("corpus_name", "name", "branch"), COMMITTED)
def test_a_committed_corpus_file_holds_the_rows_of_its_recorded_seed(corpus_name: str, name: str, branch: str):
    """The database a LoFTEE corpus commits holds exactly the rows the generator writes from
    the corpus's own cache under the seed its provenance records, whatever library wrote
    its bytes; the recorded counts and digest are the file's; the rows of the table its
    branch of the plugin queries carry each of the four outcomes. An entry that records
    no branch is a `master` file."""
    corpus = GOLDEN / corpus_name
    entry = json.loads((corpus / "provenance.json").read_text(encoding="utf-8"))["plugin_fixtures"][name]
    committed = corpus / "plugin_data" / name
    assert entry["generator"] == lcd.GENERATOR
    assert entry.get("branch", "master") == branch
    exon_counts = lcd.read_coding_exon_counts(corpus / "json_cache")
    phylocsf, gerp, counts = lcd.database_rows(exon_counts, entry["seed"])
    assert lcd.read_rows(committed) == (phylocsf, gerp)
    queried = phylocsf
    if branch == "grch38":
        queried = lcd.summary_rows(exon_counts, entry["seed"])
        assert lcd.read_summary_rows(committed) == queried
        counts = {**counts, "phylocsf_summary": len(queried)}
    assert entry["rows"] == counts
    assert entry["files"] == {f"plugin_data/{name}": hashlib.sha256(committed.read_bytes()).hexdigest()}
    keys = {(row[0], row[1]) for row in queried}
    assert len(keys) == len(queried) < counts["exons"], "some exons have no row"
    assert {lcd.phylocsf_outcome((row[2], row[3])) for row in queried} == {lcd.NO_FLAG, lcd.UNLIKELY_ORF, lcd.WEAK}
