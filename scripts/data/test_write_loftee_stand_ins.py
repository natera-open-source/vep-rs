"""Tests for write_loftee_stand_ins.py on a made-up module tree: the files written and the subroutines
each defines, what a written subroutine does under perl, the opening comment, a second run, a dry
run and each refusal."""

from __future__ import annotations

import builtins
import io
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "write_loftee_stand_ins.py"
sys.path.insert(0, str(HERE))
import write_loftee_stand_ins as writer

PREFIX = "ERROR: [write_loftee_stand_ins] "
ALPHA, BETA = "maxEntScan/alpha.pl", "maxEntScan/beta.pl"
STOP_LINE = (
    "ERROR: [loftee_stand_in] weigh_item was called: MaxEntScan scoring is not available in the stand-in "
    f"for {ALPHA}, which holds no MaxEntScan content. Run LoFTEE with its splice predictions off, or "
    "with its own maxEntScan directory.\n"
)
PERL = shutil.which("perl")
needs_perl = pytest.mark.skipif(PERL is None, reason="no perl on PATH")

# A made-up module. Its constructor calls two subroutines on paths under the directory; its helper
# file requires two files there and calls two more subroutines behind the `turbo` switch.
MODULE = r"""package Gadget;

require "helper.pl";

sub new
{
    my ($class, $home) = @_;
    my $gadget = { turbo => 0 };
    my @shelves = &load_shelves("$home/maxEntScan/shelves/");
    my %labels = &load_labels($home . '/maxEntScan/labels.txt');
    $gadget->{shelves} = \@shelves;
    $gadget->{labels} = \%labels;
    return bless $gadget, $class;
}

sub describe {
    my ($gadget, $item) = @_;
    return "plain $item" unless $gadget->{turbo};
    return "rated " . rate_item($gadget, $item);
}

1;
"""
HELPER = r"""require "maxEntScan/alpha.pl";
require 'maxEntScan/beta.pl';

sub rate_item {
    my ($gadget, $item) = @_;
    return &weigh_item($item) + &grade_item($item, $gadget->{shelves});
}

1;
"""
SWITCH_ON = 'my $gadget = Gadget->new($ARGV[0]); $gadget->{turbo} = 1; '


def plant(root: Path, module: str = MODULE, helper: str = HELPER) -> Path:
    tree = root / "gadget"
    tree.mkdir()
    (tree / "Gadget.pm").write_text(module, encoding="utf-8")
    (tree / "helper.pl").write_text(helper, encoding="utf-8")
    return tree


def run(capsys, tree: Path, *flags: str) -> tuple[int, str, str]:
    code = writer.main([str(tree), *flags])
    captured = capsys.readouterr()
    return code, captured.out, captured.err


def listing(tree: Path) -> list[str]:
    return sorted(path.relative_to(tree).as_posix() for path in tree.rglob("*"))


def state(tree: Path) -> dict[str, tuple[int, int]]:
    """The size and the modification time of every path of the tree."""
    found = {path: os.lstat(tree / path) for path in listing(tree)}
    return {path: (status.st_size, status.st_mtime_ns) for path, status in found.items()}


def defined(path: Path) -> list[str]:
    return re.findall(r"^sub (\w+) ", path.read_text(encoding="utf-8"), flags=re.MULTILINE)


def perl(tree: Path, program: str) -> subprocess.CompletedProcess:
    command = [PERL, "-I", str(tree), "-MGadget", "-e", program, str(tree)]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def record_opens(monkeypatch) -> list[str]:
    """Every path handed to `open` and `os.open` from here on."""
    opened: list[str] = []

    def recording(real):
        def wrapper(file, *args, **kwargs):
            opened.append(str(file) if isinstance(file, int) else os.fsdecode(file))
            return real(file, *args, **kwargs)

        return wrapper

    monkeypatch.setattr(builtins, "open", recording(builtins.open))
    monkeypatch.setattr(io, "open", recording(io.open))
    monkeypatch.setattr(os, "open", recording(os.open))
    return opened


def test_writes_a_file_per_require_and_every_subroutine_in_the_first(tmp_path, capsys):
    tree = plant(tmp_path)
    code, out, err = run(capsys, tree)
    assert (code, err) == (0, "")
    assert out == f"LoFTEE stand-ins: wrote {ALPHA} (4 subroutines), {BETA} (0 subroutines); left none\n"
    assert listing(tree) == ["Gadget.pm", "helper.pl", "maxEntScan", ALPHA, BETA]
    assert defined(tree / ALPHA) == ["load_shelves", "load_labels", "weigh_item", "grade_item"]
    assert defined(tree / BETA) == []
    text = (tree / ALPHA).read_text(encoding="utf-8")
    for name in ("load_shelves", "load_labels"):
        assert f"sub {name} {{ return (); }}\n" in text
    for name in ("weigh_item", "grade_item"):
        assert f"sub {name} {{ $stand_in_stop->('{name}'); }}\n" in text


def test_the_script_runs_from_the_command_line(tmp_path):
    tree = plant(tmp_path)
    result = subprocess.run([sys.executable, str(SCRIPT), str(tree)], capture_output=True, text=True, check=False)
    assert (result.returncode, result.stderr) == (0, "")
    assert result.stdout.startswith(f"LoFTEE stand-ins: wrote {ALPHA} (4 subroutines), {BETA}")
    assert (tree / ALPHA).is_file() and (tree / BETA).is_file()


@needs_perl
def test_constructor_called_subroutines_return_an_empty_list_and_the_module_runs(tmp_path, capsys):
    tree = plant(tmp_path)
    assert run(capsys, tree)[0] == 0
    program = (
        'my $gadget = Gadget->new($ARGV[0]); my @direct = Gadget::load_shelves("anything"); '
        'print join("|", $gadget->describe("cup"), scalar(@{$gadget->{shelves}}), '
        'scalar(keys %{$gadget->{labels}}), scalar(@direct));'
    )
    result = perl(tree, program)
    assert (result.returncode, result.stdout, result.stderr) == (0, "plain cup|0|0|0", "")


@needs_perl
def test_switch_only_subroutine_ends_the_process_with_the_message(tmp_path, capsys):
    tree = plant(tmp_path)
    assert run(capsys, tree)[0] == 0
    # A `die` would be trapped by the eval and the program would print; the stand-in exits.
    result = perl(tree, SWITCH_ON + 'eval { $gadget->describe("cup") }; print "survived";')
    assert (result.returncode, result.stdout, result.stderr) == (1, "", STOP_LINE)


@needs_perl
def test_message_reaches_standard_error_after_the_program_reopens_it(tmp_path, capsys):
    tree = plant(tmp_path)
    assert run(capsys, tree)[0] == 0
    reopen = "close(STDERR); my $held = ''; open(STDERR, '>', \\$held); "
    result = perl(tree, SWITCH_ON + reopen + '$gadget->describe("cup"); print "survived";')
    assert (result.returncode, result.stdout, result.stderr) == (1, "", STOP_LINE)


def test_each_file_opens_with_the_three_statements(tmp_path, capsys):
    tree = plant(tmp_path)
    assert run(capsys, tree)[0] == 0
    for path, line in ((ALPHA, 1), (BETA, 2)):
        opening = (tree / path).read_text(encoding="utf-8").splitlines()[:3]
        assert opening == [
            "# This file is a stand-in written by write_loftee_stand_ins.py.",
            f"# It replaces {path}, which helper.pl line {line} requires.",
            "# It holds no MaxEntScan content: no model, no table and no line of MaxEntScan's code.",
        ]


def test_second_run_leaves_its_own_files_and_writes_a_missing_one(tmp_path, capsys, monkeypatch):
    tree = plant(tmp_path)
    assert run(capsys, tree)[0] == 0
    before = {path: ((tree / path).read_bytes(), (tree / path).stat().st_mtime_ns) for path in (ALPHA, BETA)}
    opened = record_opens(monkeypatch)
    again = run(capsys, tree)
    monkeypatch.undo()
    assert again == (0, f"LoFTEE stand-ins: wrote none; left {ALPHA}, {BETA}\n", "")
    assert str(tree / "helper.pl") in opened
    assert [path for path in opened if path.startswith(str(tree / "maxEntScan"))] == []
    assert before == {path: ((tree / path).read_bytes(), (tree / path).stat().st_mtime_ns) for path in before}
    (tree / BETA).unlink()
    assert run(capsys, tree) == (0, f"LoFTEE stand-ins: wrote {BETA} (0 subroutines); left {ALPHA}\n", "")
    assert (tree / BETA).read_bytes() == before[BETA][0]
    assert (tree / ALPHA).stat().st_mtime_ns == before[ALPHA][1]


def test_foreign_file_refuses_and_is_never_opened(tmp_path, capsys, monkeypatch):
    tree = plant(tmp_path)
    directory = tree / "maxEntScan"
    directory.mkdir()
    # One foreign file under a name the script would write, one under another name; neither readable.
    for name in ("alpha.pl", "notes.txt"):
        (directory / name).write_text("not the script's own\n", encoding="utf-8")
        (directory / name).chmod(0o000)
    before = state(tree)
    opened = record_opens(monkeypatch)
    code, out, err = run(capsys, tree)
    monkeypatch.undo()
    assert (code, out) == (1, "")
    assert err == (
        f"{PREFIX}{tree} holds {ALPHA}, which this script did not write "
        "(2 such under maxEntScan); nothing is written\n"
    )
    assert str(tree / "Gadget.pm") in opened and str(tree / "helper.pl") in opened
    assert [path for path in opened if path.startswith(str(directory))] == []
    assert before == state(tree)


def test_checkout_without_a_call_site_refuses(tmp_path, capsys):
    module = "package Gadget;\n\nsub new { return bless {}, shift; }\n\n1;\n"
    tree = plant(tmp_path, module=module, helper="sub rate_item { return 2; }\n\n1;\n")
    code, out, err = run(capsys, tree)
    assert (code, out) == (1, "")
    assert err == f"{PREFIX}{tree}: no Perl file requires a path under maxEntScan/: not a LoFTEE checkout\n"
    assert listing(tree) == ["Gadget.pm", "helper.pl"]


def test_missing_checkout_refuses(tmp_path, capsys):
    absent = tmp_path / "absent"
    assert run(capsys, absent) == (1, "", f"{PREFIX}{absent} is not a directory\n")
    assert not absent.exists()


def test_path_the_module_reads_itself_refuses(tmp_path, capsys):
    reading = '    open(my $table, "<", "$home/maxEntScan/table.txt");\n    return bless $gadget, $class;\n'
    tree = plant(tmp_path, module=MODULE.replace("    return bless $gadget, $class;\n", reading))
    code, out, err = run(capsys, tree)
    assert (code, out) == (1, "")
    assert err == (
        f"{PREFIX}{tree}: Gadget.pm line 13 reads maxEntScan/table.txt outside a call to a subroutine "
        "the stand-ins define; a stand-in cannot supply that file\n"
    )
    assert listing(tree) == ["Gadget.pm", "helper.pl"]


def test_dry_run_prints_the_files_and_subroutines_and_writes_nothing(tmp_path, capsys):
    tree = plant(tmp_path)
    code, out, err = run(capsys, tree, "--dry-run")
    assert (code, err) == (0, "")
    assert listing(tree) == ["Gadget.pm", "helper.pl"]
    assert run(capsys, tree)[0] == 0
    sizes = {path: (tree / path).stat().st_size for path in (ALPHA, BETA)}
    assert out.splitlines() == [
        f"{ALPHA}: would write {sizes[ALPHA]} bytes; replaces the file helper.pl line 1 requires",
        "  load_shelves: returns an empty list to the constructor (called at Gadget.pm line 9)",
        "  load_labels: returns an empty list to the constructor (called at Gadget.pm line 10)",
        "  weigh_item: ends the run (called at helper.pl line 6)",
        "  grade_item: ends the run (called at helper.pl line 6)",
        f"{BETA}: would write {sizes[BETA]} bytes; replaces the file helper.pl line 2 requires",
        "  no subroutine",
        "maxEntScan/shelves/: named at Gadget.pm line 9, read by no stand-in, not written",
        "maxEntScan/labels.txt: named at Gadget.pm line 10, read by no stand-in, not written",
    ]
    assert run(capsys, tree, "--dry-run")[1].splitlines()[0].startswith(f"{ALPHA}: would leave {sizes[ALPHA]} bytes")


def test_defined_imported_and_commented_names_are_not_stood_in(tmp_path, capsys):
    helper = r"""use Made::Up qw(tally);
require "maxEntScan/alpha.pl";

=pod

require "maxEntScan/in_pod.pl";

=cut

sub own_rule { return 2; }

sub rate_item {
    my ($gadget, $item) = @_;
    # &ghost_rule("maxEntScan/in_comment.pl");
    my $note = "text with &quoted_rule inside";
    return &tally($item) + &own_rule() + &weigh_item($item);
}

1;
__END__
require "maxEntScan/after_end.pl";
"""
    tree = plant(tmp_path, helper=helper)
    assert run(capsys, tree)[0] == 0
    assert listing(tree) == ["Gadget.pm", "helper.pl", "maxEntScan", ALPHA]
    assert defined(tree / ALPHA) == ["load_shelves", "load_labels", "weigh_item"]


def test_subroutine_also_called_outside_the_constructor_ends_the_run(tmp_path, capsys):
    again = '    return "plain $item" unless $gadget->{turbo};\n    my %more = &load_labels("maxEntScan/more.txt");\n'
    tree = plant(tmp_path, module=MODULE.replace('    return "plain $item" unless $gadget->{turbo};\n', again))
    assert run(capsys, tree)[0] == 0
    text = (tree / ALPHA).read_text(encoding="utf-8")
    assert "sub load_shelves { return (); }\n" in text
    assert "sub load_labels { $stand_in_stop->('load_labels'); }\n" in text
