#!/usr/bin/env python3
"""Tests for reference_release.py, the derivation of the comparators' --reference-release
from a reference set's provenance.

Run with: ``python3 -m pytest -q scripts/concordance/test_reference_release.py``
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from reference_release import (
    CACHE_VERSION_RELEASE,
    REFERENCE_RELEASES,
    main,
    provenance_file,
    release_from_constants,
    release_from_image,
    release_from_provenance,
)

HERE = Path(__file__).resolve().parent


class ReleaseFromImageTests(unittest.TestCase):
    def test_the_two_registry_releases_are_read_from_the_tag(self) -> None:
        self.assertEqual(release_from_image("ensemblorg/ensembl-vep:release_116.2"), "116.2")
        self.assertEqual(release_from_image("ensemblorg/ensembl-vep:release_115.2"), "115.2")
        self.assertEqual(REFERENCE_RELEASES, ("115.2", "116.2"))

    def test_a_release_without_a_registry_is_refused_not_guessed(self) -> None:
        """116.1 shares 116.2's skip and overlap semantics, but the registries are named by
        the release they were verified against; the caller chooses with --reference-release."""
        with self.assertRaisesRegex(ValueError, "116.1"):
            release_from_image("ensemblorg/ensembl-vep:release_116.1")

    def test_a_digest_only_or_untagged_image_names_no_release(self) -> None:
        with self.assertRaises(ValueError):
            release_from_image("ensemblorg/ensembl-vep@sha256:5c57abdc40b637cac198370b4c114777fa24de33")
        with self.assertRaises(ValueError):
            release_from_image("ensemblorg/ensembl-vep:latest")


class ReleaseFromProvenanceTests(unittest.TestCase):
    def test_every_image_key_the_reference_writers_use_is_read(self) -> None:
        for key in ("perl_docker_image", "perl_image", "perl_vep_docker_image", "vep_image_tag"):
            with self.subTest(key=key):
                self.assertEqual(release_from_provenance({key: "ensemblorg/ensembl-vep:release_116.2"}), "116.2")
        self.assertEqual(
            release_from_provenance({"perl_reference": {"image": "ensemblorg/ensembl-vep:release_116.2"}}), "116.2"
        )

    def test_the_image_tag_wins_over_the_cache_version(self) -> None:
        self.assertEqual(
            release_from_provenance({"perl_docker_image": "ensemblorg/ensembl-vep:release_115.2", "cache_version": 116}),
            "115.2",
        )

    def test_a_digest_only_image_falls_through_to_the_cache_version(self) -> None:
        self.assertEqual(
            release_from_provenance({"perl_image": "ensemblorg/ensembl-vep@sha256:abc", "cache_version": 116}), "116.2"
        )
        self.assertEqual(release_from_provenance({"cache_version": "115"}), "115.2")
        self.assertEqual(release_from_provenance({"perl_reference": {"cache_version": 116}}), "116.2")

    def test_the_reported_version_line_and_a_recorded_command_are_read(self) -> None:
        self.assertEqual(release_from_provenance({"vep_version_reported": "ensembl-vep          : 116.2"}), "116.2")
        self.assertEqual(
            release_from_provenance({"vep_command": "vep -i in.vcf --offline --cache --cache_version 115 --tab"}), "115.2"
        )

    def test_a_perl_reference_block_resolves_by_image_then_prefix_then_cache_version(self) -> None:
        block = {"image": "ensemblorg/ensembl-vep:release_116.2", "digest": "sha256:abc",
                 "cache_version": 116, "prefix": "population_r116.2_c116_20261020"}
        self.assertEqual(release_from_provenance(block), "116.2")
        self.assertEqual(release_from_provenance({"prefix": "population_r116.2_c116_20261020", "cache_version": 115}), "116.2")
        self.assertEqual(release_from_provenance({"perl_reference_prefix": "population_r115.2_c115_20260930"}), "115.2")
        self.assertEqual(release_from_provenance({"cache_version": 116}), "116.2")
        with self.assertRaises(ValueError):
            release_from_provenance({"digest": "sha256:abc"})

    def test_neither_source_is_a_loud_failure(self) -> None:
        with self.assertRaisesRegex(ValueError, "neither an image tag"):
            release_from_provenance({"fork": 4, "buffer_size": 5000})
        with self.assertRaisesRegex(ValueError, "no release"):
            release_from_provenance({"cache_version": 114})

    def test_a_golden_corpus_provenance_resolves(self) -> None:
        """A golden corpus records its reference as a digest-only image and a command line;
        the command's --cache_version names the release, and it must be the release the
        corpus directory is named for (`tests/golden/<major>/`). Skipped in a tree without
        the corpora (they are read from the public tree there)."""
        corpora = sorted((HERE.parent.parent / "tests" / "golden").glob("*/GRCh37/provenance.json"))
        if not corpora:
            self.skipTest("no golden corpus is in this tree")
        for corpus in corpora:
            with self.subTest(corpus=corpus.parent.parent.name):
                self.assertEqual(
                    release_from_provenance(json.loads(corpus.read_text())),
                    CACHE_VERSION_RELEASE[corpus.parent.parent.name],
                )


class ProvenanceDirAndConstantsTests(unittest.TestCase):
    def test_provenance_json_is_preferred_over_sidecars(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "output.provenance.json").write_text('{"cache_version": 115}')
            self.assertEqual(provenance_file(root).name, "output.provenance.json")
            (root / "provenance.json").write_text('{"cache_version": 116}')
            self.assertEqual(provenance_file(root).name, "provenance.json")

    def test_the_three_ground_truth_layouts(self) -> None:
        """A reference set is laid out one of three ways: one record at the set root over
        suite directories holding only outputs (`<set>/provenance.json`, `<set>/<suite>/output.txt`);
        a sidecar beside each output (`<set>/<suite>/output.provenance.json`); or neither, which
        names no release and must fail rather than fall through to a guess."""
        with tempfile.TemporaryDirectory() as d:
            root = Path(d) / "canonical_set"
            for suite in ("clinvar_a", "clinvar_b", "gnomad_a"):
                (root / suite).mkdir(parents=True)
                (root / suite / "output.txt").write_text("#header\n")
            (root / "provenance.json").write_text(
                json.dumps({"perl_image": "ensemblorg/ensembl-vep:release_115.2", "cache_version": 115})
            )
            self.assertEqual(provenance_file(root / "clinvar_a"), root / "provenance.json")
            self.assertEqual(main(["--provenance-dir", str(root / "clinvar_a")]), 0)
            (root / "gnomad_a" / "output.provenance.json").write_text(
                json.dumps({"perl_image": "ensemblorg/ensembl-vep:release_116.2", "cache_version": 116})
            )
            self.assertEqual(provenance_file(root / "gnomad_a"), root / "gnomad_a" / "output.provenance.json",
                             "a suite's own sidecar wins over the set record")
            bare = Path(d) / "bare_set"
            (bare / "clinvar_a").mkdir(parents=True)
            (bare / "clinvar_a" / "output.txt").write_text("#header\n")
            with self.assertRaisesRegex(ValueError, "has no provenance.json"):
                provenance_file(bare / "clinvar_a")
            with self.assertRaises(ValueError):
                provenance_file(root / "missing")

    def test_a_checkouts_constants_name_its_release(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            pm = Path(d) / "Constants.pm"
            pm.write_text("our $VEP_VERSION = 116;\nour $VEP_SUB_VERSION = 2;\n")
            self.assertEqual(release_from_constants(pm), "116.2")
            pm.write_text("our $VEP_VERSION = 116;\n")
            with self.assertRaises(ValueError):
                release_from_constants(pm)


class CommandLineTests(unittest.TestCase):
    def test_prints_the_release_and_exits_zero(self) -> None:
        out = subprocess.run(
            [sys.executable, str(HERE / "reference_release.py"), "--image", "ensemblorg/ensembl-vep:release_116.2"],
            capture_output=True, text=True,
        )
        self.assertEqual((out.returncode, out.stdout.strip()), (0, "116.2"))

    def test_a_failure_is_one_error_line_on_stderr_and_exit_one(self) -> None:
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "provenance.json"
            p.write_text('{"fork": 4}')
            out = subprocess.run(
                [sys.executable, str(HERE / "reference_release.py"), "--provenance", str(p)],
                capture_output=True, text=True,
            )
        self.assertEqual(out.returncode, 1)
        self.assertEqual(out.stdout, "")
        self.assertTrue(out.stderr.startswith("ERROR: [reference_release] "), out.stderr)

    def test_an_empty_image_and_a_non_object_record_are_one_line_errors(self) -> None:
        self.assertEqual(main(["--image", ""]), 1)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "provenance.json"
            p.write_text("[1, 2]")
            self.assertEqual(main(["--provenance", str(p)]), 1)

    def test_main_takes_exactly_one_source(self) -> None:
        with self.assertRaises(SystemExit):
            main(["--image", "x", "--provenance", "y"])
        with self.assertRaises(SystemExit):
            main([])


if __name__ == "__main__":
    unittest.main(verbosity=2)
