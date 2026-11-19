"""Tests for check_divergence_order.py on a fabricated page and table, and on the committed index page's own shape."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import check_divergence_order as cdo

COLUMNS = ("entry", "kind", "category", "tier", "records_total", "records_with_shape", "per_million", "corpora_with_shape", "flag_sets_reaching", "basis", "requires", "note")


def table_row(entry, tier, total, with_shape, corpora, flag_sets, basis="population cut records"):
    per_million = f"{with_shape / total * 1_000_000:.3f}" if total else "0.000"
    return dict(zip(COLUMNS, (entry, "difference", "c", str(tier), str(total), str(with_shape), per_million, str(corpora), str(flag_sets), basis, "-", "-")))


ROWS = [
    table_row("first_rate", 1, 15_646_481, 2_155_480, 102, 379),
    table_row("second_rate", 1, 15_646_481, 114, 13, 379),
    table_row("many_sets", 2, 15_646_481, 13_491_934, 103, 15),
    table_row("fewer_sets", 2, 15_646_481, 1_543, 20, 15),
    table_row("mt_shape", 2, 487, 247, 7, 15, basis="corpora records on MT (no chromosome 21 cut carries one)"),
    table_row("alpha_synthetic", 3, 15_646_481, 0, 6, 379),
    table_row("beta_flag", 3, 379, 1, 1, 1, basis="committed flag sets (the shape is in the command line)"),
    table_row("plugin_a", 2, 15_646_481, 6_347, 62, 6),
    table_row("plugin_b", 3, 15_646_481, 0, 0, 23),
    table_row("plugin_c", 2, 1_638, 45, 1, 6, basis="rows of the committed tables plugin_data/x.tsv.gz"),
]
TABLE = {r["entry"]: r for r in ROWS}
# The entries of the malformed-input section: tier 3, their shape counted like any other entry's.
MALFORMED_ROWS = [
    table_row("delta_broken", 3, 15_646_481, 0, 0, 379),
    table_row("epsilon_broken", 3, 379, 0, 0, 4, basis="committed flag sets (a malformed or synthetic value in the command line)"),
    table_row("plugin_d", 3, 690, 0, 0, 14, basis="rows of the committed tables plugin_data/y.tsv.gz"),
    table_row("plugin_e", 3, 15_646_481, 0, 0, 23),
]
SPLIT_TABLE = {**TABLE, **{r["entry"]: r for r in MALFORMED_ROWS}}
CLASSES = {key: "unobserved" if r["tier"] == "3" else "observed" for key, r in TABLE.items()}
SPLIT_CLASSES = {**CLASSES, **{r["entry"]: "malformed" for r in MALFORMED_ROWS}}


def write_table(path: Path, rows=ROWS) -> Path:
    lines = ["# a comment line", "# cut\tdataset\tbuild\tpath\tsha256\trecords", "\t".join(COLUMNS)]
    lines += ["\t".join(r[c] for c in COLUMNS) for r in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def write_registry(directory: Path, classes=CLASSES) -> Path:
    """A registry of one record per entry, each naming its input class alone."""
    directory.mkdir()
    for key, input_class in classes.items():
        (directory / f"{key}.toml").write_text(f'entry = "{key}"\ninput_class = "{input_class}"\n')
    return directory


def row(key: str, cell: str | None = None, verified: str = "116.2") -> str:
    return f"| [`{key}`](x.md#{key}) | one sentence | {cell or cdo.render_count(SPLIT_TABLE[key])} | {verified} | [Page](x.md) |"


def line_of(text: str, needle: str) -> int:
    """The number of the first line of `text` that holds `needle`."""
    return next(n for n, line in enumerate(text.splitlines(), 1) if needle in line)


HEADER = cdo.ENTRY_TABLE_HEADER + "\n| --- | --- | --- | --- | --- |"
GLOSSARY = "## How to read the counts\n\n- **The registry**: named by identifiers such as `reference_skipped_record`.\n\n| Defect | What goes wrong | Where | Classes | Change |\n| --- | --- | --- | --- | --- |\n| A | s | w | 1, 2 | c |\n"


def page(order=("first_rate", "second_rate", "many_sets", "fewer_sets", "mt_shape", "alpha_synthetic", "beta_flag"), plugins=None, cells=None) -> str:
    cells = cells or {}
    text = "# Intended divergences\n\n## How the entries are ordered {#differences-outside-the-classes}\n\nprose\n\n## SNP and indel calling\n\n" + HEADER + "\n"
    text += "\n".join(row(k, cells.get(k)) for k in order) + "\n\n## Plugins\n\n"
    for title, keys in (plugins if plugins is not None else [("Several plugins", ["plugin_a"]), ("REVEL", ["plugin_c"]), ("SpliceAI", ["plugin_b"])]):
        text += f"### {title}\n\n" + HEADER + "\n" + "\n".join(row(k, cells.get(k)) for k in keys) + "\n\n"
    return text + GLOSSARY


def split_page(snp=("delta_broken", "epsilon_broken"), plugins=(("Several plugins", ["plugin_e"]), ("SpliceAI", ["plugin_d"])), order=("snp", "plugins"), direct=(), main: str | None = None) -> str:
    """A page with the malformed-input section after its last category: the same categories one heading level down, the plugin groups under `####`."""

    def entry_table(keys) -> str:
        return HEADER + "\n" + "\n".join(row(k) for k in keys) + "\n\n"

    categories = {
        "snp": "### SNP and indel calling\n\n" + entry_table(snp),
        "plugins": "### Plugins\n\n" + "".join(f"#### {title}\n\n" + entry_table(keys) for title, keys in plugins),
    }
    section = f"## {cdo.MALFORMED_SECTION}\n\nprose\n\n"
    if direct:
        section += entry_table(direct)
    return (main or page()).removesuffix(GLOSSARY) + section + "".join(categories[name] for name in order) + GLOSSARY


def test_render_count_states_the_tier_the_figure_and_the_reach():
    assert cdo.render_count(TABLE["first_rate"]) == "1 · 137,761 per million · no flag needed"
    assert cdo.render_count(TABLE["second_rate"]) == "1 · 7.3 per million · no flag needed"
    assert cdo.render_count(TABLE["many_sets"]) == "2 · 862,298 per million · 15 flag sets, 103 corpora"
    assert cdo.render_count(TABLE["fewer_sets"]) == "2 · 98.6 per million · 15 flag sets, 20 corpora"
    assert cdo.render_count(TABLE["mt_shape"]) == "2 · 247 of 487 MT records · 15 flag sets, 7 corpora"
    assert cdo.render_count(TABLE["alpha_synthetic"]) == "3 · none on the cuts · 379 flag sets, 6 corpora"
    assert cdo.render_count(TABLE["beta_flag"]) == "3 · 1 of 379 flag sets · 1 flag set, 1 corpus"
    assert cdo.render_count(TABLE["plugin_c"]) == "2 · 45 of 1,638 table rows · 6 flag sets, 1 corpus"
    assert cdo.render_count(table_row("rare", 1, 15_646_481, 4, 9, 379)) == "1 · 0.26 per million · no flag needed"
    assert cdo.render_count(table_row("one_set", 2, 15_646_481, 1, 1, 1)) == "2 · 0.064 per million · 1 flag set, 1 corpus"
    assert cdo.render_count(table_row("no_table", 3, 0, 0, 0, 13, basis="rows of the committed tables plugin_data/*README*")) == "3 · no table rows to count · 13 flag sets, 0 corpora"
    assert cdo.render_count(table_row("malformed", 3, 379, 1, 1, 124, basis="committed flag sets (a malformed or synthetic value in the command line)")) == "3 · 1 of 379 flag sets · 124 flag sets, 1 corpus"
    assert cdo.render_count(table_row("comma", 2, 672, 457, 6, 2, basis="epigenome names of the committed regulatory caches")) == "2 · 457 of 672 epigenome names · 2 flag sets, 6 corpora"


def test_render_count_counts_a_tool_entry_in_tool_runs():
    """The committed runs of `vep filter` and `vep recode` are not flag sets, so both the figure and the reach of a tool entry say tool runs."""
    filt = table_row("filter_shape", 3, 28, 0, 0, 1, basis="committed runs of the filter tool (outside the flag sets)")
    recode = table_row("recode_shape", 3, 4, 0, 0, 4, basis="committed runs of the recoder tool (outside the flag sets)")
    assert cdo.render_count(filt) == "3 · 0 of 28 tool runs · 1 tool run, 0 corpora"
    assert cdo.render_count(recode) == "3 · 0 of 4 tool runs · 4 tool runs, 0 corpora"
    reached = table_row("recode_hit", 2, 4, 1, 1, 4, basis="committed runs of the recoder tool (outside the flag sets)")
    assert cdo.render_count(reached) == "2 · 1 of 4 tool runs · 4 tool runs, 1 corpus"
    assert cdo.cell_order_key("recode_hit", cdo.render_count(reached)) == cdo.order_key(reached) == (2, -4, -1, "recode_hit")


def test_a_rate_counted_from_output_rows_says_so_and_orders_like_any_other_cell():
    """A shape that is a property of the annotation is counted from the output rows of the entry's stated command over the cuts: the cell says `by output rows`, the order rule reads it as any tier 2 cell, and a page stating the bare rate is reported."""
    basis = "output rows over the population cuts (records with a value of two or more items)"
    item_order = table_row("item_order", 2, 15_646_481, 2, 1, 2, basis=basis)
    cell = "2 · 0.13 per million by output rows · 2 flag sets, 1 corpus"
    assert cdo.render_count(item_order) == cell
    assert cdo.render_count(table_row("pair_order", 2, 15_646_481, 8_984, 2, 4, basis=basis)) == "2 · 574.2 per million by output rows · 4 flag sets, 2 corpora"
    assert cdo.render_count(table_row("no_row", 3, 15_646_481, 0, 0, 2, basis=basis)) == "3 · none on the cuts · 2 flag sets, 0 corpora"
    assert cdo.cell_order_key("item_order", cell) == cdo.order_key(item_order) == (2, -2, -1, "item_order")
    assert cdo.cell_order_key("k", "1 · 0.13 per million by output rows · no flag needed") == (1, -0.13)
    table = {**TABLE, "item_order": item_order}
    plugins = [("Several plugins", ["plugin_a"]), ("NearestExonJB", ["item_order"]), ("REVEL", ["plugin_c"]), ("SpliceAI", ["plugin_b"])]
    assert cdo.check(page(plugins=plugins, cells={"item_order": cell}), table) == []
    bare = "2 · 0.13 per million · 2 flag sets, 1 corpus"
    errors = cdo.check(page(plugins=plugins, cells={"item_order": bare}), table)
    assert len(errors) == 1 and f"item_order states {bare!r}, the table gives {cell!r}" in errors[0]


def test_render_count_refuses_an_unknown_basis_and_a_tier_outside_the_rule():
    with pytest.raises(ValueError, match="no unit"):
        cdo.render_count(table_row("odd", 2, 10, 1, 1, 1, basis="something else"))
    for tier in (0, 4):
        with pytest.raises(ValueError, match="is not 1, 2 or 3"):
            cdo.render_count(table_row("odd", tier, 10, 1, 1, 1))
        with pytest.raises(ValueError, match="is not 1, 2 or 3"):
            cdo.order_key(table_row("odd", tier, 10, 1, 1, 1))


def test_cell_order_key_reads_the_rule_back_from_a_rendered_cell():
    by_cell = sorted(TABLE, key=lambda k: cdo.cell_order_key(k, cdo.render_count(TABLE[k])))
    assert by_cell == sorted(TABLE, key=lambda k: cdo.order_key(TABLE[k]))
    for key, row in TABLE.items():
        if row["tier"] != "1":
            assert cdo.cell_order_key(key, cdo.render_count(row)) == cdo.order_key(row), key
    assert cdo.cell_order_key("first_rate", cdo.render_count(TABLE["first_rate"])) == (1, -137761.0)
    assert cdo.cell_order_key("k", "2 · 7,566 of 77,456 cache transcripts · 368 flag sets, 50 corpora") == (2, -368, -50, "k")
    assert cdo.cell_order_key("k", "3 · 1 of 379 flag sets · 1 flag set, 1 corpus") == (3, 0, "k")
    flags_tier2 = table_row("flag_two", 2, 379, 5, 1, 5, basis="committed flag sets (the shape is in the command line)")
    assert cdo.cell_order_key("flag_two", cdo.render_count(flags_tier2)) == cdo.order_key(flags_tier2)
    with pytest.raises(ValueError, match="not a rendered tier"):
        cdo.cell_order_key("k", "tier 2")
    with pytest.raises(ValueError, match="no rate per million"):
        cdo.cell_order_key("k", "1 · none on the cuts · no flag needed")
    with pytest.raises(ValueError, match="no flag sets and corpora"):
        cdo.cell_order_key("k", "2 · 1.0 per million · no flag needed")


def test_order_key_follows_the_rule():
    keys = sorted(TABLE, key=lambda k: cdo.order_key(TABLE[k]))
    assert keys == ["first_rate", "second_rate", "many_sets", "fewer_sets", "mt_shape", "plugin_a", "plugin_c", "alpha_synthetic", "beta_flag", "plugin_b"]


def test_a_page_in_the_tables_order_passes():
    assert cdo.check(page(), TABLE) == []


def test_split_cells_keeps_an_escaped_pipe_and_a_row_without_its_closing_pipe():
    assert cdo.split_cells("| a | b \\| c | d |") == ["a", "b | c", "d"]
    assert cdo.split_cells("| a | b | c | d") == ["a", "b", "c", "d"]
    assert cdo.split_cells("| a \\|") == ["a |"]


def test_parse_page_reads_an_indented_row_and_a_header_spaced_otherwise():
    text = page().replace(row("second_rate"), "   " + row("second_rate")).replace(HEADER, "|Entry|In one sentence|Tier and count|Verified against|Page|\n|:---|---|---|---|---|", 1)
    sections = cdo.parse_page(text)
    assert [r.key for r in sections[1].rows] == ["first_rate", "second_rate", "many_sets", "fewer_sets", "mt_shape", "alpha_synthetic", "beta_flag"]
    every_key = ["first_rate", "second_rate", "many_sets", "fewer_sets", "mt_shape", "alpha_synthetic", "beta_flag"]
    for lead in ("    ", "\t"):
        assert [r.key for r in cdo.parse_page(page().replace(row("second_rate"), lead + row("second_rate")))[1].rows] == every_key
    bare = page().replace(row("second_rate"), row("second_rate").strip("|").strip())
    assert [r.key for r in cdo.parse_page(bare)[1].rows] == every_key
    blank = page().replace(row("many_sets"), "\n" + row("many_sets"))
    assert [r.key for r in cdo.parse_page(blank)[1].rows] == ["first_rate", "second_rate"]
    under_a_heading = cdo.parse_page("## A\n\n" + HEADER + "\n" + row("first_rate") + "\n### B\nprose\n")
    assert [r.key for r in under_a_heading[0].rows] == ["first_rate"] and under_a_heading[0].subsections[0].rows == []


def test_parse_page_ignores_tables_with_another_header_and_rows_outside_a_table():
    sections = cdo.parse_page(page())
    assert [s.title for s in sections] == ["How the entries are ordered {#differences-outside-the-classes}", "SNP and indel calling", "Plugins", "How to read the counts"]
    assert [r.key for r in sections[1].rows] == ["first_rate", "second_rate", "many_sets", "fewer_sets", "mt_shape", "alpha_synthetic", "beta_flag"]
    assert [s.title for s in sections[2].subsections] == ["Several plugins", "REVEL", "SpliceAI"]
    assert sections[3].rows == []
    assert [s.title for s in cdo.parse_page("##  Plugins  \n")] == ["Plugins"]


def test_a_missing_entry_an_extra_entry_and_a_repeat_are_reported():
    text = page(order=("first_rate", "second_rate", "many_sets", "fewer_sets", "mt_shape", "alpha_synthetic"))
    assert cdo.check(text, TABLE) == ["beta_flag is in the table and not on the page"]
    extra = page().replace(row("beta_flag"), row("beta_flag") + "\n| [`unknown_key`](x.md) | s | 3 · none on the cuts · 1 flag set, 1 corpus | 116.2 | [Page](x.md) |")
    assert cdo.check(extra, TABLE) == ["line 18: unknown_key has no row in the table"]
    twice = page().replace(row("beta_flag"), row("beta_flag") + "\n" + row("first_rate"))
    errors = cdo.check(twice, TABLE)
    assert errors[0] == "line 18: first_rate is listed twice"


def test_a_row_out_of_the_tables_order_is_reported():
    swapped = page(order=("second_rate", "first_rate", "many_sets", "fewer_sets", "mt_shape", "alpha_synthetic", "beta_flag"))
    assert cdo.check(swapped, TABLE) == ["line 11: under 'SNP and indel calling' second_rate precedes first_rate, which the rule places first"]
    tier_before_rate = page(order=("first_rate", "second_rate", "many_sets", "fewer_sets", "alpha_synthetic", "mt_shape", "beta_flag"))
    assert "alpha_synthetic precedes mt_shape" in cdo.check(tier_before_rate, TABLE)[0]
    names = page(order=("first_rate", "second_rate", "many_sets", "fewer_sets", "mt_shape", "beta_flag", "alpha_synthetic"))
    assert "beta_flag precedes alpha_synthetic" in cdo.check(names, TABLE)[0]


def test_a_stale_tier_or_count_cell_is_reported():
    stale = page(cells={"many_sets": "2 · 862,298 per million · 14 flag sets, 103 corpora"})
    assert cdo.check(stale, TABLE) == ["line 13: many_sets states '2 · 862,298 per million · 14 flag sets, 103 corpora', the table gives '2 · 862,298 per million · 15 flag sets, 103 corpora'"]
    wrong_tier = page(cells={"first_rate": "2 · 137,761 per million · 379 flag sets, 102 corpora"})
    assert len(cdo.check(wrong_tier, TABLE)) == 1


def test_subsections_must_follow_several_plugins_in_alphabetical_order():
    unsorted = page(plugins=[("Several plugins", ["plugin_a"]), ("SpliceAI", ["plugin_b"]), ("REVEL", ["plugin_c"])])
    assert cdo.check(unsorted, TABLE) == ["line 19: the subsections of 'Plugins' are not in alphabetical order"]
    case_ignored = page(plugins=[("Several plugins", ["plugin_a"]), ("dbNSFP", ["plugin_c"]), ("REVEL", ["plugin_b"])])
    assert cdo.check(case_ignored, TABLE) == []
    first_not_general = page(plugins=[("REVEL", ["plugin_c"]), ("Several plugins", ["plugin_a"]), ("SpliceAI", ["plugin_b"])])
    assert cdo.check(first_not_general, TABLE) == ["line 19: 'Several plugins' is not the first subsection of 'Plugins'"]
    # A subsection headed twice would hide the order of its rows, here `plugin_b` before `plugin_c`.
    twice = page(plugins=[("Several plugins", ["plugin_a"]), ("REVEL", ["plugin_b"]), ("REVEL", ["plugin_c"])])
    assert cdo.check(twice, TABLE) == ["line 19: 'Plugins' lists a subsection twice"]
    # A subsection that lists no entry takes no part in the order.
    prose = page().replace("### REVEL\n", "### Zeta notes\n\nprose\n\n### REVEL\n")
    assert "### Zeta notes" in prose and cdo.check(prose, TABLE) == []


def test_a_section_cannot_list_entries_both_directly_and_in_subsections():
    text = page().replace("## Plugins\n\n", "## Plugins\n\n" + HEADER + "\n" + row("plugin_a") + "\n\n", 1)
    text = text.replace("### Several plugins\n\n" + HEADER + "\n" + row("plugin_a") + "\n\n", "")
    assert "lists entries both directly and in subsections" in cdo.check(text, TABLE)[0]
    # A subsection that lists no entry leaves its category listing directly.
    noted = page().replace("## Plugins\n", "### Notes\n\nprose\n\n## Plugins\n")
    assert "### Notes" in noted and cdo.check(noted, TABLE) == []


def test_a_page_with_the_malformed_input_section_passes():
    """The section holds the malformed entries under the page's own categories, one heading level down, and the check reads them like any other row."""
    text = split_page()
    assert cdo.check(text, SPLIT_TABLE, "116.2", SPLIT_CLASSES) == []
    sections = cdo.parse_page(text)
    malformed = next(top for top in sections if cdo.is_malformed_section(top))
    assert [s.title for s in malformed.subsections] == ["SNP and indel calling", "Plugins"]
    assert [r.key for r in malformed.subsections[0].rows] == ["delta_broken", "epsilon_broken"]
    assert [(s.title, [r.key for r in s.rows]) for s in malformed.subsections[1].subsections] == [("Several plugins", ["plugin_e"]), ("SpliceAI", ["plugin_d"])]
    inside = [f"{cdo.MALFORMED_SECTION}, {title}" for title in ("SNP and indel calling", "Several plugins", "SpliceAI")]
    assert [name for name, _ in cdo.leaves(sections)] == ["SNP and indel calling", "Several plugins", "REVEL", "SpliceAI", *inside]
    with_id = text.replace(f"## {cdo.MALFORMED_SECTION}", f"## {cdo.MALFORMED_SECTION} {{#malformed-input}}")
    assert cdo.check(with_id, SPLIT_TABLE, "116.2", SPLIT_CLASSES) == []
    assert cdo.is_malformed_section(cdo.parse_page(f"## {cdo.MALFORMED_SECTION} {{#malformed-input}}  \n")[0])
    section = f"## {cdo.MALFORMED_SECTION}\n\n"
    for orphan in (section + "#### B\n", "## A\n\n### B\n\n" + section + "#### D\n"):
        with pytest.raises(ValueError, match="a `####` heading outside a `###` subsection"):
            cdo.parse_page(orphan)
    assert cdo.check("### A\n", {}) == ["line 1: a subsection before any section"]
    assert [s.title for s in cdo.walk(cdo.parse_page(section + "### B\n\n#### C\n\n##### D\n")[0])] == [cdo.MALFORMED_SECTION, "B", "C"]


def test_the_order_rule_applies_inside_the_malformed_input_section():
    swapped = split_page(snp=("epsilon_broken", "delta_broken"))
    assert cdo.check(swapped, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(swapped, 'epsilon_broken')}: under '{cdo.MALFORMED_SECTION}, SNP and indel calling' epsilon_broken precedes delta_broken, which the rule places first"]
    unsorted = split_page(plugins=(("SpliceAI", ["plugin_d"]), ("REVEL", ["plugin_e"])))
    assert cdo.check(unsorted, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(unsorted, '### Plugins')}: the subsections of 'Plugins' are not in alphabetical order"]
    not_first = split_page(plugins=(("REVEL", ["plugin_d"]), ("Several plugins", ["plugin_e"])))
    assert cdo.check(not_first, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(not_first, '### Plugins')}: 'Several plugins' is not the first subsection of 'Plugins'"]
    both = split_page(plugins=(("SpliceAI", ["plugin_d"]),)).replace("### Plugins\n\n", "### Plugins\n\n" + HEADER + "\n" + row("plugin_e") + "\n\n")
    assert cdo.check(both, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(both, '### Plugins')}: 'Plugins' lists entries both directly and in subsections"]
    # An id on a subsection heading leaves its place in the alphabet alone.
    with_group_ids = split_page(plugins=(("Several plugins", ["plugin_e"]), ("REVEL", ["plugin_d"]))).replace("#### Several plugins\n", "#### Several plugins {#m}\n").replace("### Several plugins\n", "### Several plugins {#s}\n")
    assert cdo.check(with_group_ids, SPLIT_TABLE, classes=SPLIT_CLASSES) == []


def test_a_fourth_level_heading_outside_the_malformed_input_section_splits_no_table():
    """Outside the malformed-input section a `####` line is body text: it opens no section, and the order rule reads a subsection's rows across it."""
    assert [s.title for s in cdo.walk(cdo.parse_page("## A\n\n### B\n\n#### C\n")[0])] == ["A", "B"]
    assert cdo.parse_page("## A\n\n#### B\n")[0].subsections == []
    assert [s.title for s in cdo.parse_page("#### B\n\n## A\n")] == ["A"]
    assert [s.title for s in cdo.walk(cdo.parse_page(f"## {cdo.MALFORMED_SECTION}\n\n## Z\n\n### B\n\n#### C\n")[1])] == ["Z", "B"]
    hidden = page(plugins=[("Several plugins", ["plugin_a"]), ("REVEL", ["plugin_b"])]).replace(row("plugin_b") + "\n\n", row("plugin_b") + "\n\n#### Later\n\n" + HEADER + "\n" + row("plugin_c") + "\n\n")
    assert cdo.check(hidden, TABLE) == [f"line {line_of(hidden, 'plugin_b')}: under 'REVEL' plugin_b precedes plugin_c, which the rule places first"]
    # Directly under a row, with no blank line between, the line is read as a row and refused.
    glued = page().replace(row("beta_flag") + "\n", row("beta_flag") + "\n#### Later\n")
    assert cdo.check(glued, TABLE) == [f"line {line_of(glued, '#### Later')}: an entry row with 1 cells"]


def test_the_malformed_input_section_lists_its_entries_under_the_pages_categories_in_their_order_and_last():
    heading = f"## {cdo.MALFORMED_SECTION}"
    out_of_order = f"the categories of '{cdo.MALFORMED_SECTION}' do not stand in the order they have above it"
    reordered = split_page(order=("plugins", "snp"))
    assert cdo.check(reordered, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(reordered, heading)}: {out_of_order}"]
    with_ids = split_page(order=("plugins", "snp"), main=page().replace("## SNP and indel calling\n", "## SNP and indel calling {#snp}\n"))
    assert cdo.check(with_ids, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(with_ids, heading)}: {out_of_order}"]
    ids_inside = reordered.replace("### SNP and indel calling\n", "### SNP and indel calling {#snp-malformed}\n")
    assert ids_inside != reordered and cdo.check(ids_inside, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(ids_inside, heading)}: {out_of_order}"]
    direct = split_page(snp=("epsilon_broken",), direct=("delta_broken",))
    assert cdo.check(direct, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(direct, heading)}: '{cdo.MALFORMED_SECTION}' lists entries outside a category"]
    assert cdo.MALFORMED_SECTION in [name for name, _ in cdo.leaves(cdo.parse_page(direct))]
    for subsection in ("", "### Sub\n\n"):
        late = split_page(snp=("delta_broken",)).replace(GLOSSARY, "## Late\n\n" + subsection + HEADER + "\n" + row("epsilon_broken") + "\n\n" + GLOSSARY)
        assert cdo.check(late, SPLIT_TABLE, classes={**SPLIT_CLASSES, "epsilon_broken": "unobserved"}) == [f"line {line_of(late, '## Late')}: 'Late' lists entries after '{cdo.MALFORMED_SECTION}', the last section that does"]
    twice = split_page(snp=("delta_broken",)).replace("### Plugins\n\n####", "### SNP and indel calling\n\n" + HEADER + "\n" + row("epsilon_broken") + "\n\n### Plugins\n\n####")
    assert cdo.check(twice, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(twice, heading)}: '{cdo.MALFORMED_SECTION}' lists a category twice"]
    invented = split_page().replace("### SNP and indel calling\n\n" + HEADER, "### SNP and indel calling\n\n#### One\n\n" + HEADER)
    assert cdo.check(invented, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(invented, '### SNP and indel calling')}: 'SNP and indel calling' lists entries in subsections, which the category above does not have"]
    # Each `###` heading of the section repeats the heading of a `##` section above it, so one they lack, a misspelling among them, is refused.
    renamed = split_page().replace("### SNP and indel calling\n\n" + HEADER + "\n" + row("delta_broken"), "### SNP and indel caling\n\n" + HEADER + "\n" + row("delta_broken"))
    assert cdo.check(renamed, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(renamed, '### SNP and indel caling')}: 'SNP and indel caling' is no category of the tables above '{cdo.MALFORMED_SECTION}'"]
    # A heading that lists no entry is held to the rule too, and a `##` section below the section is no category of the tables above it.
    for title in ("Notes", "How to read the counts"):
        closing = split_page().replace(GLOSSARY, f"### {title}\n\n" + GLOSSARY)
        assert cdo.check(closing, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(closing, f'### {title}')}: {title!r} is no category of the tables above '{cdo.MALFORMED_SECTION}'"]
    # The order is held past the second category.
    third = split_page().replace(GLOSSARY, "### How the entries are ordered\n\n" + GLOSSARY)
    assert cdo.check(third, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(third, heading)}: {out_of_order}"]
    # The page heads the section once at most: a second heading of the title would hold prose no check reads.
    headed_twice = split_page().replace(GLOSSARY, heading + "\n\nprose\n\n" + GLOSSARY)
    second = [n for n, text in enumerate(headed_twice.splitlines(), 1) if text == heading][1]
    assert cdo.check(headed_twice, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {second}: '{cdo.MALFORMED_SECTION}' is headed twice"]


def test_a_category_of_the_malformed_input_section_lists_its_entries_as_the_category_above_does():
    """A category that groups its entries above groups them in the section, under `####` headings that repeat headings of its subsections; a `####` heading without entries is held to the names too."""
    above = f"in the tables above '{cdo.MALFORMED_SECTION}'"
    misspelt = split_page(plugins=(("Several plugins", ["plugin_e"]), ("SpliceAIs", ["plugin_d"])))
    assert cdo.check(misspelt, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(misspelt, '#### SpliceAIs')}: 'SpliceAIs' is no subsection of 'Plugins' {above}"]
    noted = split_page().replace("### Plugins\n\n####", "#### Notes\n\n### Plugins\n\n####")
    assert cdo.check(noted, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(noted, '#### Notes')}: 'Notes' is no subsection of 'SNP and indel calling' {above}"]
    flat = split_page(plugins=()).replace("### Plugins\n\n", "### Plugins\n\n" + HEADER + "\n" + row("plugin_d") + "\n" + row("plugin_e") + "\n\n")
    assert cdo.check(flat, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(flat, '### Plugins')}: 'Plugins' lists entries directly, which the category above lists in subsections"]
    # A subsection above that lists no entry leaves the category listing directly.
    noted_above = split_page(main=page().replace("## Plugins\n", "### Notes\n\nprose\n\n## Plugins\n"))
    assert "### Notes\n\nprose\n\n## Plugins" in noted_above and cdo.check(noted_above, SPLIT_TABLE, classes=SPLIT_CLASSES) == []
    # A category of the section need not repeat every group of the category above, and one that lists no entry repeats none.
    assert cdo.check(split_page(plugins=(("SpliceAI", ["plugin_d", "plugin_e"]),)), SPLIT_TABLE, classes=SPLIT_CLASSES) == []
    assert cdo.check_layout(cdo.parse_page(split_page(plugins=()))) == []
    # A group repeats a subsection of its own category above, one that lists no entry there included, and none of another category.
    emptied = split_page(plugins=(("Several plugins", ["plugin_e"]), ("Zeta notes", ["plugin_d"])), main=page().replace("### REVEL\n", "### Zeta notes\n\nprose\n\n### REVEL\n"))
    assert "#### Zeta notes" in emptied and cdo.check(emptied, SPLIT_TABLE, classes=SPLIT_CLASSES) == []
    borrowed = split_page(main=page().replace("## Plugins\n", "### Notes\n\nprose\n\n## Plugins\n")).replace(GLOSSARY, "#### Notes\n\n" + GLOSSARY)
    assert cdo.check(borrowed, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(borrowed, '#### Notes')}: 'Notes' is no subsection of 'Plugins' {above}"]
    # A group headed twice in the section is refused as it is above it.
    twice = split_page(plugins=(("SpliceAI", ["plugin_e"]), ("SpliceAI", ["plugin_d"])))
    assert cdo.check(twice, SPLIT_TABLE, classes=SPLIT_CLASSES) == [f"line {line_of(twice, '### Plugins')}: 'Plugins' lists a subsection twice"]


def test_every_malformed_entry_stands_in_the_malformed_input_section_with_tier_3():
    text = split_page()
    outside = f"outside '{cdo.MALFORMED_SECTION}'"
    stray = {**SPLIT_CLASSES, "alpha_synthetic": "malformed"}
    assert cdo.check(text, SPLIT_TABLE, classes=stray) == [f"line {line_of(text, 'alpha_synthetic')}: alpha_synthetic is malformed and stands under 'SNP and indel calling', {outside}"]
    in_a_subsection = {**SPLIT_CLASSES, "plugin_b": "malformed"}
    assert cdo.check(text, SPLIT_TABLE, classes=in_a_subsection) == [f"line {line_of(text, 'plugin_b')}: plugin_b is malformed and stands under 'SpliceAI', {outside}"]
    # A shape the measured data carries is observed whatever rule it breaks, so a tier 2 row has no place in the section.
    main = page(order=("first_rate", "second_rate", "many_sets", "fewer_sets", "alpha_synthetic", "beta_flag"))
    seen = split_page(snp=("mt_shape", "delta_broken", "epsilon_broken"), main=main)
    assert cdo.check(seen, SPLIT_TABLE, classes={**SPLIT_CLASSES, "mt_shape": "malformed"}) == [f"line {line_of(seen, 'mt_shape')}: mt_shape is tier 2 with input_class 'malformed': a tier 1 or tier 2 entry is observed"]


def test_every_tier_1_or_2_entry_is_observed_and_a_tier_3_entry_may_be():
    """Tier 1 or tier 2 makes an entry observed, never the converse: a tier 3 entry whose page records its shape in real data is observed in its own table."""
    text = split_page()
    for key, tier in (("first_rate", 1), ("many_sets", 2), ("plugin_c", 2)):
        errors = cdo.check(text, SPLIT_TABLE, classes={**SPLIT_CLASSES, key: "unobserved"})
        assert errors == [f"line {line_of(text, key)}: {key} is tier {tier} with input_class 'unobserved': a tier 1 or tier 2 entry is observed"]
    for key in ("beta_flag", "plugin_b"):
        assert cdo.CELL.match(cdo.render_count(SPLIT_TABLE[key])).group("tier") == "3"
        assert cdo.check(text, SPLIT_TABLE, classes={**SPLIT_CLASSES, key: "observed"}) == []


def test_no_observed_or_unobserved_entry_stands_in_the_malformed_input_section():
    text = split_page()
    for key, input_class in (("delta_broken", "unobserved"), ("plugin_d", "unobserved"), ("plugin_d", "observed")):
        errors = cdo.check(text, SPLIT_TABLE, classes={**SPLIT_CLASSES, key: input_class})
        assert errors == [f"line {line_of(text, key)}: {key} stands under '{cdo.MALFORMED_SECTION}' with input_class {input_class!r}"]


def test_every_entry_the_opening_of_the_malformed_input_section_names_stands_outside_it():
    """The section's opening prose points at entries of the tables above that cover malformed input too; a linked key that is no entry, or a key of the section itself, linked or bare, is refused."""

    def opening(prose: str) -> str:
        return split_page().replace(f"## {cdo.MALFORMED_SECTION}\n\nprose\n", f"## {cdo.MALFORMED_SECTION}\n\n{prose}\n")

    def errors(prose: str) -> list[str]:
        return cdo.check(opening(prose), SPLIT_TABLE, classes=SPLIT_CLASSES)

    def pointing(key: str) -> str:
        return f"One entry above, [`{key}`](x.md#{key}), covers such input too."

    text = opening(pointing("beta_flag"))
    line = line_of(text, "One entry above")
    inside = f"which the opening of '{cdo.MALFORMED_SECTION}' names as an entry of the tables above, stands in that section"
    unknown = f"which the opening of '{cdo.MALFORMED_SECTION}' names, is no entry of the page"
    # An entry of the tables above is a row of a category or of one of its subsections.
    for key in ("beta_flag", "plugin_b"):
        assert errors(pointing(key)) == []
    assert [top.named for top in cdo.parse_page(text)] == [[], [], [], [("beta_flag", line, True)], []]
    # An entry of the section is one of its category's rows or of a `####` group's.
    for key in ("delta_broken", "plugin_d"):
        assert errors(pointing(key)) == [f"line {line}: {key}, {inside}"]
    assert errors(pointing("no_such_entry")) == [f"line {line}: no_such_entry, {unknown}"]
    # Every pointer of a line is read, and one in a later paragraph of the opening.
    assert errors(pointing("beta_flag") + " So does [`plugin_d`](x.md).") == [f"line {line}: plugin_d, {inside}"]
    assert errors(pointing("beta_flag") + "\n\nSo does [`no_such_entry`](x.md).") == [f"line {line + 2}: no_such_entry, {unknown}"]
    # A bare key is read too, so an entry of the section is not named without its link; a bare word that is no entry is prose.
    assert errors("The entry `delta_broken` covers such input too.") == [f"line {line}: delta_broken, {inside}"]
    assert errors("The entries `beta_flag` and `no_such_entry`, the class `malformed` and [`reevaluate.py`](x.md) are named.") == []
    # A key in brackets with no target in parentheses behind it is bare.
    assert errors("The entry [`no_such_entry`][1] is named.") == []
    assert errors("The entry [`delta_broken`][1] is named.") == [f"line {line}: delta_broken, {inside}"]
    # A key linked in another section's prose, or under a category of this one, is no pointer.
    stray = "[`no_such_entry`](x.md) in prose\n\n"
    elsewhere = split_page().replace("prose\n\n## SNP and indel calling", stray + "## SNP and indel calling")
    below = split_page().replace("### SNP and indel calling\n\n" + HEADER, "### SNP and indel calling\n\n" + stray + HEADER)
    for unread in (elsewhere, below):
        assert stray in unread and cdo.check(unread, SPLIT_TABLE, classes=SPLIT_CLASSES) == []
        assert [section.named for top in cdo.parse_page(unread) for section in cdo.walk(top) if section.named] == []
    repeated = cdo.parse_page(f"## A\n\n### {cdo.MALFORMED_SECTION}\n\n{stray}")
    assert [section.named for section in cdo.walk(repeated[0])] == [[], []]


def test_an_entry_without_an_input_class_is_reported_when_the_classes_are_read():
    text = split_page()
    classes = {key: value for key, value in SPLIT_CLASSES.items() if key != "plugin_b"}
    assert cdo.check(text, SPLIT_TABLE, classes=classes) == [f"line {line_of(text, 'plugin_b')}: plugin_b has no input class in the registry"]
    assert cdo.check(text, SPLIT_TABLE) == []
    assert cdo.check_classes(cdo.parse_page(text), {}) == []
    # A row whose cell states no tier is held to its section alone.
    no_tier = page(cells={"many_sets": "tier 2"})
    assert cdo.check_classes(cdo.parse_page(no_tier), {**CLASSES, "many_sets": "unobserved"}) == []
    assert cdo.check_classes(cdo.parse_page(no_tier), {**CLASSES, "many_sets": "malformed"}) == [f"line {line_of(no_tier, 'many_sets')}: many_sets is malformed and stands under 'SNP and indel calling', outside '{cdo.MALFORMED_SECTION}'"]
    inside = text.replace(row("delta_broken"), row("delta_broken", "tier 3"))
    assert cdo.check_classes(cdo.parse_page(inside), {**SPLIT_CLASSES, "delta_broken": "unobserved"}) == [f"line {line_of(inside, 'delta_broken')}: delta_broken stands under '{cdo.MALFORMED_SECTION}' with input_class 'unobserved'"]
    assert cdo.check(no_tier, TABLE, classes=CLASSES) == [f"line {line_of(no_tier, 'many_sets')}: many_sets states 'tier 2', the table gives '2 · 862,298 per million · 15 flag sets, 103 corpora'"]


def test_read_classes_reads_the_input_class_of_every_record(tmp_path):
    registry = write_registry(tmp_path / "registry", SPLIT_CLASSES)
    assert cdo.read_classes(registry) == SPLIT_CLASSES
    record = registry / "first_rate.toml"
    record.write_text('entry = "first_rate"\n')
    with pytest.raises(ValueError, match=r"first_rate\.toml: input_class is None, not one of observed, unobserved, malformed"):
        cdo.read_classes(registry)
    record.write_text('entry = "first_rate"\ninput_class = "synthetic"\n')
    with pytest.raises(ValueError, match=r"first_rate\.toml: input_class is 'synthetic', not one of"):
        cdo.read_classes(registry)
    for not_toml in (b"entry = [\n", b'input_class = "observed\xff"\n'):
        record.write_bytes(not_toml)
        with pytest.raises(ValueError, match=r"first_rate\.toml: "):
            cdo.read_classes(registry)
    assert cdo.read_classes(tmp_path / "absent") == {}


def test_a_malformed_row_or_a_headless_table_is_one_error():
    headless = HEADER + "\n" + row("first_rate") + "\n"
    assert cdo.check(headless, TABLE) == ["line 1: an entry table before any heading"]
    bad_cell = page().replace("[`first_rate`](x.md#first_rate)", "[first_rate](x.md)")
    assert cdo.check(bad_cell, TABLE) == ["line 11: the first cell is not a linked key in code font: '[first_rate](x.md)'"]
    trailing = page().replace("[`first_rate`](x.md#first_rate)", "[`first_rate`](x.md#first_rate), class 3")
    assert cdo.check(trailing, TABLE) == ["line 11: the first cell is not a linked key in code font: '[`first_rate`](x.md#first_rate), class 3'"]
    short = page().replace(row("first_rate"), "| [`first_rate`](x.md) | s | t |")
    assert cdo.check(short, TABLE) == ["line 11: an entry row with 3 cells"]
    four = page().replace(row("first_rate"), "| [`first_rate`](x.md) | s | 1 · 137,761 per million · no flag needed | [Page](x.md) |")
    assert cdo.check(four, TABLE) == ["line 11: an entry row with 4 cells"]


def test_the_verified_against_cell_is_a_release_and_matches_the_reference():
    """Every row names the Ensembl VEP release its reproduction was last run against; `--reference` pins the one release."""
    assert cdo.check(page(), TABLE, "116.2") == []
    words = page().replace(row("many_sets"), row("many_sets", verified="yes"))
    assert cdo.check(words, TABLE) == ["line 13: many_sets is verified against 'yes', not a release number"]
    empty = page().replace(row("many_sets"), row("many_sets", verified=""))
    assert cdo.check(empty, TABLE) == ["line 13: many_sets is verified against '', not a release number"]
    older = page().replace(row("fewer_sets"), row("fewer_sets", verified="115.2"))
    assert cdo.check(older, TABLE) == []
    assert cdo.check(older, TABLE, "116.2") == ["line 14: fewer_sets is verified against 115.2, the reference is 116.2"]
    assert cdo.check(page().replace(row("beta_flag"), row("beta_flag", verified="117")), TABLE) == []


def test_read_table_needs_the_columns_and_unique_keys(tmp_path):
    table = write_table(tmp_path / "t.tsv")
    assert set(cdo.read_table(str(table))) == set(TABLE)
    (tmp_path / "dup.tsv").write_text("\t".join(COLUMNS) + "\n" + "\t".join(ROWS[0][c] for c in COLUMNS) + "\n" + "\t".join(ROWS[0][c] for c in COLUMNS) + "\n")
    with pytest.raises(ValueError, match="appears twice"):
        cdo.read_table(str(tmp_path / "dup.tsv"))
    (tmp_path / "narrow.tsv").write_text("entry\ttier\nx\t1\n")
    with pytest.raises(ValueError, match="lacks"):
        cdo.read_table(str(tmp_path / "narrow.tsv"))
    (tmp_path / "short.tsv").write_text("\t".join(COLUMNS) + "\nx\t1\n")
    with pytest.raises(ValueError, match="a row with 2 cells under a header of 12"):
        cdo.read_table(str(tmp_path / "short.tsv"))
    (tmp_path / "empty.tsv").write_text("# nothing\n")
    with pytest.raises(ValueError, match="no header"):
        cdo.read_table(str(tmp_path / "empty.tsv"))


def test_main_exits_0_on_a_matching_page_and_1_otherwise(tmp_path, capsys, monkeypatch):
    table = write_table(tmp_path / "t.tsv")
    registry_args = ["--reproductions", str(write_registry(tmp_path / "registry"))]
    good = tmp_path / "good.md"
    good.write_text(page())
    summary = "10 entries match the table's 10 rows and stand in its order: 7 observed, 3 unobserved, 0 malformed"
    assert cdo.main([str(good), str(table), *registry_args]) == 0
    assert summary in capsys.readouterr().out
    wider = write_registry(tmp_path / "wider", {**CLASSES, "record_without_a_row": "malformed"})
    assert cdo.main([str(good), str(table), "--reproductions", str(wider)]) == 0
    assert summary in capsys.readouterr().out
    monkeypatch.setattr(cdo, "REPRODUCTIONS", tmp_path / "registry")
    assert cdo.main([str(good), str(table)]) == 0
    assert summary in capsys.readouterr().out
    bad = tmp_path / "bad.md"
    bad.write_text(page(order=("second_rate", "first_rate", "many_sets", "fewer_sets", "mt_shape", "alpha_synthetic", "beta_flag")))
    assert cdo.main([str(bad), str(table), *registry_args]) == 1
    assert "ERROR: [divergence-order]" in capsys.readouterr().err
    older = tmp_path / "older.md"
    older.write_text(page().replace(row("beta_flag"), row("beta_flag", verified="115.2")))
    assert cdo.main([str(older), str(table), *registry_args]) == 0
    capsys.readouterr()
    assert cdo.main(["--reference", "116.2", str(older), str(table), *registry_args]) == 1
    assert "the reference is 116.2" in capsys.readouterr().err
    (tmp_path / "empty.tsv").write_text("# nothing\n")
    assert cdo.main([str(good), str(tmp_path / "empty.tsv"), *registry_args]) == 1
    odd = [dict(r, basis="something else") if r["entry"] == "mt_shape" else r for r in ROWS]
    assert cdo.main([str(good), str(write_table(tmp_path / "odd.tsv", odd)), *registry_args]) == 1
    assert "no unit for the basis" in capsys.readouterr().err
    assert cdo.main([str(tmp_path / "missing.md"), str(table), *registry_args]) == 1
    assert "ERROR: [divergence-order]" in capsys.readouterr().err
    assert cdo.main([str(good), str(tmp_path), *registry_args]) == 1
    assert "ERROR: [divergence-order]" in capsys.readouterr().err


def test_main_reads_each_entrys_input_class_from_the_registry(tmp_path, capsys):
    """The command holds the page to the registry's classes: a class the page's sections contradict, a record without one and an entry without a record each exit 1."""
    table = write_table(tmp_path / "t.tsv", [*ROWS, *MALFORMED_ROWS])
    split = tmp_path / "split.md"
    split.write_text(split_page())
    registry = write_registry(tmp_path / "registry", SPLIT_CLASSES)
    assert cdo.main([str(split), str(table), "--reproductions", str(registry)]) == 0
    assert "14 entries match the table's 14 rows and stand in its order: 7 observed, 3 unobserved, 4 malformed" in capsys.readouterr().out
    moved = write_registry(tmp_path / "moved", {**SPLIT_CLASSES, "delta_broken": "unobserved"})
    assert cdo.main([str(split), str(table), "--reproductions", str(moved)]) == 1
    assert f"delta_broken stands under '{cdo.MALFORMED_SECTION}' with input_class 'unobserved'" in capsys.readouterr().err
    (moved / "delta_broken.toml").write_text('entry = "delta_broken"\ninput_class = "synthetic"\n')
    assert cdo.main([str(split), str(table), "--reproductions", str(moved)]) == 1
    assert "delta_broken.toml: input_class is 'synthetic', not one of observed, unobserved, malformed" in capsys.readouterr().err
    (moved / "delta_broken.toml").unlink()
    assert cdo.main([str(split), str(table), "--reproductions", str(moved)]) == 1
    assert "delta_broken has no input class in the registry" in capsys.readouterr().err


def test_the_committed_index_page_parses_and_lists_every_entry_once():
    text = (HERE.parents[1] / "docs" / "intended-divergences.md").read_text(encoding="utf-8")
    sections = cdo.parse_page(text)
    listing = cdo.leaves(sections)
    keys = [r.key for _, s in listing for r in s.rows]
    assert len(keys) >= 145 and len(keys) == len(set(keys))
    assert all(cdo.RELEASE.match(r.cells[cdo.VERIFIED_CELL]) for _, s in listing for r in s.rows), "every row names the release it was verified against"
    assert not any(k.endswith(("_sentence", "_sentences")) for k in keys), "one row per mechanism"
    for name, s in listing:
        keys_in_order = [r.key for r in s.rows]
        assert keys_in_order == sorted(keys_in_order, key=lambda k: cdo.cell_order_key(k, next(r.cells[cdo.COUNT_CELL] for r in s.rows if r.key == k))), name
    plugins = next(top for top in sections if top.title == "Plugins")
    titles = [s.title for s in plugins.subsections]
    assert titles[0] == cdo.FIRST_SUBSECTION and titles[1:] == sorted(titles[1:], key=str.casefold)
    assert cdo.check_layout(sections) == []
    assert cdo.check_named(sections) == []


def spell(count: int) -> str:
    """A count below a thousand in the page's words: `twenty-two`, `one hundred and ninety-seven`."""
    units = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen")
    tens = ("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety")
    if count < 20:
        return units[count]
    if count < 100:
        return tens[count // 10 - 2] + (f"-{units[count % 10]}" if count % 10 else "")
    return f"{units[count // 100]} hundred" + (f" and {spell(count % 100)}" if count % 100 else "")


def test_spell_writes_a_count_as_the_page_does():
    assert [spell(n) for n in (0, 7, 15, 19, 20, 22, 40, 99, 100, 197, 213, 990)] == ["zero", "seven", "fifteen", "nineteen", "twenty", "twenty-two", "forty", "ninety-nine", "one hundred", "one hundred and ninety-seven", "two hundred and thirteen", "nine hundred and ninety"]


def test_the_committed_index_page_places_every_entry_by_its_input_class_and_counts_its_malformed_input_section():
    """The committed page and registry agree on every entry's class, and the page's opening counts the entries of its malformed-input section, the defects that produce a class apart from the differences outside the classes."""
    text = (HERE.parents[1] / "docs" / "intended-divergences.md").read_text(encoding="utf-8")
    sections = cdo.parse_page(text)
    classes = cdo.read_classes()
    assert {r.key for _, s in cdo.leaves(sections) for r in s.rows} == set(classes)
    assert cdo.check_classes(sections, classes) == []
    malformed = next(top for top in sections if cdo.is_malformed_section(top))
    inside = [r for s in cdo.walk(malformed) for r in s.rows]
    assert inside and {r.key for r in inside} == {key for key, input_class in classes.items() if input_class == "malformed"}
    # A row whose page cell names a class is one of the defects that produce the concordance classes.
    defects = [r for r in inside if re.search(r"\), class(?:es)? [0-9]", r.cells[4])]
    assert all(r.cells[4].startswith("[Defect ") for r in defects)
    counted = f"{spell(len(inside) - len(defects)).capitalize()} of these differences"
    if defects:
        counted += f", and {spell(len(defects))} of the defects,"
    assert f"{counted} are reached by malformed input alone" in " ".join(text.split())
