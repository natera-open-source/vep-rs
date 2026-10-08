#!/usr/bin/env python3
"""Derive the comparators' ``--reference-release`` from a reference set's provenance.

The release is read from the Ensembl VEP image tag the reference was produced with
(``ensemblorg/ensembl-vep:release_116.2`` gives ``116.2``), under whichever key the
provenance record carries it; a record with no image tag falls back to its cache
version (115 gives 115.2, 116 gives 116.2, the releases the comparators carry a registry
for). A record with neither, or a tag of a release without a registry, is an error, so a
reference of unknown release is never scored under a guessed one.

usage: reference_release.py --provenance FILE
       reference_release.py --provenance-dir DIR      (provenance.json, the first *.provenance.json,
                                                      or the parent directory's provenance.json)
       reference_release.py --image TAG
       reference_release.py --vep-constants Constants.pm   (a local ensembl-vep checkout)

Prints the release on stdout; exits 1 with ``ERROR: [reference_release] ...`` on failure.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

REFERENCE_RELEASES = ("115.2", "116.2")
CACHE_VERSION_RELEASE = {"115": "115.2", "116": "116.2"}
IMAGE_KEYS = ("perl_docker_image", "perl_image", "perl_vep_docker_image", "vep_image_tag", "image")
_TAG = re.compile(r"release_(\d+\.\d+)$")
# `ensembl-vep          : 116.2`, the version line VEP prints and a golden corpus records.
_REPORTED = re.compile(r"ensembl-vep\s*:\s*(\d+\.\d+)")
_CACHE_VERSION_FLAG = re.compile(r"--cache_version\s+(\d+)")
# `population_r116.2_c116_<date>`, the archived reference set a record names as its prefix.
_PREFIX = re.compile(r"_r(\d+\.\d+)_")


def release_from_image(image: str) -> str:
    """The release a ``release_X.Y`` image tag names; raises ValueError for any other tag
    or for a release the comparators have no registry for."""
    m = _TAG.search(image.strip())
    if not m:
        raise ValueError(f"image {image!r} carries no release_X.Y tag")
    release = m.group(1)
    if release not in REFERENCE_RELEASES:
        raise ValueError(
            f"image {image!r} names release {release}, which the comparators have no registry "
            f"for (one of {', '.join(REFERENCE_RELEASES)}); pass --reference-release to choose one"
        )
    return release


def release_from_provenance(record: dict) -> str:
    """The release of a provenance record or of a ``perl_reference`` block: its image tag
    first (top level or under ``perl_reference``; a digest-only image reference names no
    release and is passed over), else the version line VEP reported, else the release in
    its archived prefix (``population_r116.2_c116_<date>``), else its cache version (the
    ``cache_version`` field, or the ``--cache_version`` of a recorded command);
    ValueError when none is present."""
    candidates = [record.get(k) for k in IMAGE_KEYS]
    nested = record.get("perl_reference")
    if isinstance(nested, dict):
        candidates.append(nested.get("image"))
    for image in candidates:
        if isinstance(image, str) and _TAG.search(image):
            return release_from_image(image)
    reported = record.get("vep_version_reported")
    if isinstance(reported, str) and (m := _REPORTED.search(reported)):
        return release_from_image(f"release_{m.group(1)}")
    for prefix in (record.get("prefix"), record.get("perl_reference_prefix")):
        if isinstance(prefix, str) and (m := _PREFIX.search(prefix)):
            return release_from_image(f"release_{m.group(1)}")
    cache_version = record.get("cache_version")
    if cache_version is None and isinstance(nested, dict):
        cache_version = nested.get("cache_version")
    if cache_version is None:
        command = record.get("vep_command") or record.get("vep_flags")
        if isinstance(command, str) and (m := _CACHE_VERSION_FLAG.search(command)):
            cache_version = m.group(1)
    if cache_version is not None:
        release = CACHE_VERSION_RELEASE.get(str(cache_version).strip())
        if release is None:
            raise ValueError(
                f"cache_version {cache_version!r} maps to no release the comparators have a registry for"
            )
        return release
    raise ValueError(
        f"provenance carries neither an image tag ({', '.join(IMAGE_KEYS)}, perl_reference.image) "
        "nor a cache_version"
    )


def release_from_constants(path: Path) -> str:
    """The release of an ensembl-vep checkout from its Constants.pm."""
    text = path.read_text(encoding="utf-8", errors="replace")
    major = re.search(r"\$VEP_VERSION\s*=\s*(\d+)", text)
    minor = re.search(r"\$VEP_SUB_VERSION\s*=\s*(\d+)", text)
    if not major or not minor:
        raise ValueError(f"{path} carries no $VEP_VERSION / $VEP_SUB_VERSION")
    return release_from_image(f"release_{major.group(1)}.{minor.group(1)}")


def provenance_file(directory: Path) -> Path:
    """The provenance record that names a ground-truth directory's release: its own
    ``provenance.json``, else the first ``*.provenance.json`` sidecar beside an output, else
    the ``provenance.json`` of the parent directory, which is where a reference set written
    as one record for every suite keeps it (``<set>/provenance.json`` over ``<set>/<suite>/``)."""
    if not directory.is_dir():
        raise ValueError(f"{directory} is not a directory")
    direct = directory / "provenance.json"
    if direct.is_file():
        return direct
    sidecars = sorted(directory.glob("*.provenance.json"))
    if sidecars:
        return sidecars[0]
    parent = directory.parent / "provenance.json"
    if parent.is_file():
        return parent
    raise ValueError(
        f"{directory} holds no provenance.json and no *.provenance.json, and {directory.parent} has no provenance.json"
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--provenance", type=Path, help="a provenance.json")
    src.add_argument("--provenance-dir", type=Path, help="a directory holding provenance.json or *.provenance.json")
    src.add_argument("--image", help="an Ensembl VEP image reference, e.g. ensemblorg/ensembl-vep:release_116.2")
    src.add_argument("--vep-constants", type=Path, help="modules/Bio/EnsEMBL/VEP/Constants.pm of a checkout")
    a = ap.parse_args(argv)
    try:
        if a.image is not None:
            release = release_from_image(a.image)
        elif a.vep_constants is not None:
            release = release_from_constants(a.vep_constants)
        else:
            path = a.provenance or provenance_file(a.provenance_dir)
            record = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(record, dict):
                raise ValueError(f"{path} is not a JSON object")
            release = release_from_provenance(record)
    except (ValueError, OSError) as exc:
        print(f"ERROR: [reference_release] {exc}", file=sys.stderr)
        return 1
    print(release)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
