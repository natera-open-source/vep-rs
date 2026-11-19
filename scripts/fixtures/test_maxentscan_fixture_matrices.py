"""Tests for the MaxEntScan fixture generator and its scoring reference."""

from __future__ import annotations

import hashlib
import math
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import maxentscan_fixture_matrices as mfm  # noqa: E402

# Every table of these tests is a made-up stand-in that the generator or a test writes.
EXPECTED_LINES = {
    "me2x5": 16384,
    "splicemodels/splice5sequences": 16384,
    "splicemodels/me2x3acc1": 16384,
    "splicemodels/me2x3acc2": 16384,
    "splicemodels/me2x3acc3": 16384,
    "splicemodels/me2x3acc4": 16384,
    "splicemodels/me2x3acc5": 16384,
    "splicemodels/me2x3acc6": 64,
    "splicemodels/me2x3acc7": 256,
    "splicemodels/me2x3acc8": 64,
    "splicemodels/me2x3acc9": 256,
}


@pytest.fixture(scope="module")
def seed_one(tmp_path_factory) -> tuple[Path, dict]:
    out = tmp_path_factory.mktemp("maxentscan") / "fordownload_synthetic"
    summary = mfm.write_directory(out, seed=1, force=False)
    return out, summary


def test_hashseq_is_the_plugins_base_four_hash():
    # MaxEntScan.pm:805 documents this value for its own `score3_hashseq`.
    assert mfm.hashseq("CAGAAGT") == 4619
    assert mfm.hashseq("AAA") == 0
    assert mfm.hashseq("TTTT") == 255
    assert mfm.hashseq("cag") == mfm.hashseq("CAG")


def test_all_kmers_enumerate_in_hash_order():
    kmers = mfm.all_kmers(3)
    assert len(kmers) == 64
    assert kmers[0] == "AAA" and kmers[-1] == "TTT"
    assert [mfm.hashseq(k) for k in kmers] == list(range(64))


def test_files_have_the_exact_line_counts_and_value_range(seed_one):
    out, summary = seed_one
    for name, count in EXPECTED_LINES.items():
        lines = (out / name).read_text().splitlines()
        assert len(lines) == count, name
        assert summary["files"][name]["lines"] == count
        if name.endswith("splice5sequences"):
            continue
        values = [float(line) for line in lines]
        assert all(mfm.LOW <= v <= mfm.HIGH for v in values), name


def test_plain_lines_are_the_five_place_decimals_and_the_eleven_tables_alone(tmp_path):
    summary = mfm.write_directory(tmp_path, seed=1, force=False, plain_lines=True)
    assert summary["other_files"] == {}
    assert sorted(p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("*") if p.is_file()) == sorted(
        EXPECTED_LINES
    )
    for name in EXPECTED_LINES:
        lines = (tmp_path / name).read_text().splitlines()
        assert lines == mfm.file_lines(1)[name], name
        if not name.endswith("splice5sequences"):
            assert all(line == f"{float(line):.5f}" for line in lines), name


def test_every_spelling_of_a_line_is_the_value_of_its_plain_line(seed_one):
    out, _ = seed_one
    plain = mfm.file_lines(1)
    for name in EXPECTED_LINES:
        lines = (out / name).read_text().split("\n")[:-1]
        if name.endswith("splice5sequences"):
            assert lines == plain[name]
            continue
        assert [float(line) for line in lines] == [float(line) for line in plain[name]], name
        # Each spelling is written in every number table: the decimal, and an exponent of
        # either sign.
        bare = [line.strip(" \t") for line in lines]
        assert any(line == f"{float(line):.5f}" for line in bare), name
        assert any("e-" in line for line in bare) and any("e+" in line for line in bare), name
        assert all(line == f"{float(line):.5f}" or line == f"{float(line):.6e}" for line in bare), name


def test_the_first_table_stands_between_blanks_and_the_others_do_not(seed_one):
    out, _ = seed_one
    lines = (out / "me2x5").read_text().split("\n")[:-1]
    around = {(line[: len(line) - len(line.lstrip(" \t"))], line[len(line.rstrip(" \t")) :]) for line in lines}
    assert around == set(mfm.BLANKS)
    # One removal of whitespace, which the plugin makes (MaxEntScan.pm:739-741, 755-757, 833-835),
    # leaves blanks on most lines.
    assert sum(len(before) + len(after) >= 2 for before, after in mfm.BLANKS) >= 5
    for name in EXPECTED_LINES:
        if name != "me2x5":
            assert all(line == line.strip() for line in (out / name).read_text().split("\n")[:-1]), name


def test_files_the_plugin_never_opens_stand_beside_the_tables(seed_one):
    out, summary = seed_one
    assert sorted(summary["other_files"]) == sorted(mfm.OTHER_FILES)
    assert not set(mfm.OTHER_FILES) & set(EXPECTED_LINES)
    for name, entry in summary["other_files"].items():
        data = (out / name).read_bytes()
        assert entry == {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    # A file named as a table, in the directory the plugin does not read that table from
    # (MaxEntScan.pm:104-105), is not that table: it holds another line count.
    assert len((out / "splicemodels/me2x5").read_text().splitlines()) != EXPECTED_LINES["me2x5"]
    assert len((out / "splice5sequences").read_text().splitlines()) != 16384


def test_scores_do_not_depend_on_the_spelling(seed_one, tmp_path):
    out, _ = seed_one
    mfm.write_directory(tmp_path, seed=1, force=False, plain_lines=True)
    spelled, plain = mfm.Tables(out), mfm.Tables(tmp_path)
    assert spelled.me2x5 == plain.me2x5
    assert spelled.acceptor == plain.acceptor
    assert spelled.sequence_index == plain.sequence_index


def test_splice5sequences_is_a_shuffled_permutation_of_every_7mer(seed_one):
    out, _ = seed_one
    sequences = (out / "splicemodels/splice5sequences").read_text().splitlines()
    assert sorted(sequences) == mfm.all_kmers(7)
    assert sequences != mfm.all_kmers(7)


def test_summary_digests_are_the_files_bytes(seed_one):
    out, summary = seed_one
    assert summary["seed"] == 1
    for name, entry in summary["files"].items():
        assert entry["sha256"] == hashlib.sha256((out / name).read_bytes()).hexdigest()


def test_same_seed_same_bytes_other_seed_other_bytes(tmp_path):
    a = mfm.file_lines(7)
    b = mfm.file_lines(7)
    c = mfm.file_lines(8)
    assert a == b
    assert mfm.spelled_lines(7) == mfm.spelled_lines(7)
    assert mfm.spelled_lines(7)["me2x5"] != mfm.spelled_lines(8)["me2x5"]
    assert mfm.other_files(7) == mfm.other_files(7)
    assert a["me2x5"] != c["me2x5"]
    assert a["splicemodels/splice5sequences"] != c["splicemodels/splice5sequences"]


def test_existing_files_need_force(tmp_path):
    mfm.write_directory(tmp_path, seed=1, force=False)
    with pytest.raises(FileExistsError):
        mfm.write_directory(tmp_path, seed=2, force=False)
    before = (tmp_path / "me2x5").read_bytes()
    mfm.write_directory(tmp_path, seed=2, force=True)
    assert (tmp_path / "me2x5").read_bytes() != before


def write_tables(directory: Path, me2x5, sequences, acceptor):
    """Writes a directory from explicit lists: `me2x5` numbers in `sequences` order, and the
    nine acceptor tables as lists of numbers in hash order."""
    (directory / "splicemodels").mkdir(parents=True)
    (directory / "me2x5").write_text("".join(f"{v}\n" for v in me2x5))
    (directory / "splicemodels/splice5sequences").write_text("".join(f"{s}\n" for s in sequences))
    for (name, _, _), table in zip(mfm.ACCEPTOR_TABLES, acceptor):
        (directory / "splicemodels" / name).write_text("".join(f"{v}\n" for v in table))


def test_score5_reads_the_consensus_and_the_other_seven_bases_by_sequence_line(tmp_path):
    # Every 7-mer scores 1.0 except AAACCCG (the bases at 0, 1, 2, 5, 6, 7, 8 of
    # AAAGTCCCG), placed on the third line of a shuffled sequence file.
    sequences = mfm.all_kmers(7)
    sequences.remove("AAACCCG")
    sequences.insert(2, "AAACCCG")
    me2x5 = [1.0] * 16384
    me2x5[2] = 8.0
    acceptor = [[1.0] * 4**k for _, k, _ in mfm.ACCEPTOR_TABLES]
    write_tables(tmp_path, me2x5, sequences, acceptor)
    tables = mfm.Tables(tmp_path)

    consensus = 0.9896 * 0.9884 / (0.23 * 0.27)
    assert tables.score5("AAAGTCCCG") == math.log(consensus * 8.0) / math.log(2)
    # The same seven bases around a non-consensus pair read the same table entry.
    weak = 0.004 * 0.0034 / (0.27 * 0.27)
    assert tables.score5("AAAAACCCG") == math.log(weak * 8.0) / math.log(2)
    # Another 7-mer reads the 1.0 entry: the consensus factor alone.
    assert tables.score5("CCCGTAAAA") == math.log(consensus) / math.log(2)


def test_score3_multiplies_five_tables_over_four_at_the_plugins_offsets(tmp_path):
    # rest = seq[0:18] + seq[20:23]; the nine tables read rest[0:7], [7:14], [14:21],
    # [4:11], [11:18], [4:7], [7:11], [11:14], [14:18] (MaxEntScan.pm:698-726, 851-859, 867).
    # Each table carries 1.0 everywhere but the entry its slice of this k-mer hashes to.
    seq = "ACGTACGTACGTACGTACAGTTT"
    rest = seq[0:18] + seq[20:23]
    marks = [2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0]
    acceptor = []
    for (_, k, part), mark in zip(mfm.ACCEPTOR_TABLES, marks):
        table = [1.0] * 4**k
        table[mfm.hashseq(rest[part])] = mark
        acceptor.append(table)
    write_tables(tmp_path, [1.0] * 16384, mfm.all_kmers(7), acceptor)
    tables = mfm.Tables(tmp_path)

    consensus = 0.9903 * 0.9905 / (0.27 * 0.23)
    ratio = 2.0 * 3.0 * 5.0 * 7.0 * 11.0 / (13.0 * 17.0 * 19.0 * 23.0)
    assert tables.score3(seq) == math.log(consensus * ratio) / math.log(2)


def test_score_dispatches_on_length_and_refuses_other_lengths(seed_one):
    out, _ = seed_one
    tables = mfm.Tables(out)
    assert tables.score("CAGGTAAGT") == tables.score5("CAGGTAAGT")
    assert tables.score("T" * 18 + "AGGTA") == tables.score3("T" * 18 + "AGGTA")
    with pytest.raises(ValueError):
        tables.score("ACGT")


def test_format_score_is_three_decimals():
    assert mfm.format_score(8.1514) == "8.151"
    assert mfm.format_score(-13.86049) == "-13.860"
    assert mfm.format_score(0.0) == "0.000"


def test_boundary_guard_names_scores_and_pairs_at_a_rounding_boundary():
    assert mfm.boundary_distance(1.2345) == pytest.approx(0.0, abs=1e-12)
    assert mfm.at_boundary(1.2345)
    assert not mfm.at_boundary(1.2340)
    assert not mfm.at_boundary(1.2349)


def test_check_kmers_reports_the_offending_specs(seed_one, monkeypatch):
    out, _ = seed_one
    tables = mfm.Tables(out)
    scores = {"CAGGTAAGT": 1.0, "AAGGTAAGT": 1.0005, "T" * 18 + "AGGTA": 2.0005}
    monkeypatch.setattr(tables, "score", lambda k: scores[k])
    assert mfm.check_kmers(tables, ["CAGGTAAGT"]) == []
    assert mfm.check_kmers(tables, ["AAGGTAAGT"]) == ["AAGGTAAGT"]
    # A pair is checked on both scores and on their difference.
    assert mfm.check_kmers(tables, ["CAGGTAAGT,T" + "T" * 17 + "AGGTA"]) == [
        "CAGGTAAGT,T" + "T" * 17 + "AGGTA"
    ]
    scores["AAGGTAAGT"] = 1.5005
    scores["T" * 18 + "AGGTA"] = 1.2
    assert mfm.check_kmers(tables, ["AAGGTAAGT,CAGGTAAGT"]) == ["AAGGTAAGT,CAGGTAAGT"]


def test_cli_writes_scores_and_checks(tmp_path, capsys):
    out = tmp_path / "fordownload_synthetic"
    assert mfm.main(["--seed", "1", "--out", str(out)]) == 0
    captured = capsys.readouterr().out
    assert '"seed": 1' in captured
    assert (out / "me2x5").exists()
    assert mfm.main(["--seed", "1", "--out", str(out)]) == 1
    assert mfm.main(["--dir", str(out), "--score", "CAGGTAAGT"]) == 0
    line = capsys.readouterr().out.strip()
    kmer, score = line.split("\t")
    assert kmer == "CAGGTAAGT"
    assert score == mfm.format_score(mfm.Tables(out).score5("CAGGTAAGT"))
    spec = tmp_path / "kmers.txt"
    spec.write_text("CAGGTAAGT\nCAGGTAAGT,AAGGTAAGT\n\n")
    assert mfm.main(["--dir", str(out), "--check-kmers", str(spec)]) == 0
