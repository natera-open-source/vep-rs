#!/usr/bin/env python3
"""Read the glibc version each Linux release binary requires and hold it to the floor the release page states.

A dynamically linked executable names the symbol versions it needs in its `.gnu.version_r` section, one
`GLIBC_x.y` entry per glibc version it draws a symbol from (`GLIBC_ABI_DT_RELR`, the packed-relocation
marker, counts as 2.36); the highest is the oldest glibc the binary runs on. Every path is an ELF file or a
`.tar.gz` archive whose ELF members are all read. One line per binary:

    vep-0.5.0-x86_64-unknown-linux-gnu.tar.gz/vep: glibc 2.17

With `--floor X.Y`, a binary needing a newer glibc fails the run and names the symbols that need it.
`--json FILE` writes the measured floor per path (the highest over an archive's binaries), the input of
`render_release_notes.py --glibc-floors`.

usage: check_glibc_floor.py [--floor X.Y] [--json FILE] PATH...

No dependency beyond the standard library.
"""

from __future__ import annotations

import argparse
import json
import re
import struct
import sys
import tarfile
from pathlib import Path

GLIBC_VERSION = re.compile(r"^GLIBC_(\d+(?:\.\d+)*)$")
# A GLIBC version need that names no number: a binary with packed relative relocations loads only on glibc 2.36
# or later, and the dynamic loader checks it under this name.
GLIBC_ABI_VERSIONS = {"GLIBC_ABI_DT_RELR": "2.36"}
SHT_GNU_VERSYM = 0x6FFFFFFF
SHT_GNU_VERNEED = 0x6FFFFFFE
SHT_DYNSYM = 11


def fail(message: str) -> None:
    sys.exit(f"ERROR: [check_glibc_floor] {message}")


def version_key(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


class Elf:
    """The section table of one ELF file, 32- or 64-bit, either byte order."""

    def __init__(self, data: bytes, name: str) -> None:
        if data[:4] != b"\x7fELF":
            raise ValueError(f"{name} is not an ELF file")
        if data[4] not in (1, 2) or data[5] not in (1, 2):
            raise ValueError(f"{name} has an ELF class or byte order this reader does not know")
        self.data = data
        self.name = name
        self.bits = {1: 32, 2: 64}[data[4]]
        self.order = {1: "<", 2: ">"}[data[5]]
        # e_shoff, then past e_flags, e_ehsize, e_phentsize and e_phnum to e_shentsize and e_shnum.
        if self.bits == 64:
            shoff, shentsize, shnum = struct.unpack_from(f"{self.order}Q10xHH", data, 0x28)
        else:
            shoff, shentsize, shnum = struct.unpack_from(f"{self.order}I10xHH", data, 0x20)
        self.sections = [self.section(shoff + i * shentsize) for i in range(shnum)]

    def section(self, at: int) -> dict[str, int]:
        if self.bits == 64:
            fields = struct.unpack_from(f"{self.order}IIQQQQIIQQ", self.data, at)
        else:
            fields = struct.unpack_from(f"{self.order}IIIIIIIIII", self.data, at)
        _, kind, _, _, offset, size, link, info, _, entsize = fields
        return {"type": kind, "offset": offset, "size": size, "link": link, "info": info, "entsize": entsize}

    def u16(self, at: int) -> int:
        return struct.unpack_from(f"{self.order}H", self.data, at)[0]

    def u32(self, at: int) -> int:
        return struct.unpack_from(f"{self.order}I", self.data, at)[0]

    def string(self, table: dict[str, int], index: int) -> str:
        start = table["offset"] + index
        end = self.data.find(b"\0", start)
        if end < 0:
            raise ValueError(f"{self.name} has a string table entry past the end of the file")
        return self.data[start:end].decode("ascii", "replace")

    def find(self, kind: int) -> dict[str, int] | None:
        return next((s for s in self.sections if s["type"] == kind), None)

    def version_needs(self) -> dict[int, str]:
        """Every version this file requires from a shared library, keyed by the index `.gnu.version` uses."""
        verneed = self.find(SHT_GNU_VERNEED)
        if verneed is None:
            return {}
        strings = self.sections[verneed["link"]]
        needs: dict[int, str] = {}
        at = verneed["offset"]
        # Elf_Verneed: vn_cnt at 2, vn_aux at 8, vn_next at 12. Each Elf_Vernaux it points to: vna_other at 6,
        # vna_name at 8, vna_next at 12.
        for libraries_left in range(verneed["info"], 0, -1):
            count, aux, following = self.u16(at + 2), self.u32(at + 8), self.u32(at + 12)
            entry = at + aux
            for versions_left in range(count, 0, -1):
                needs[self.u16(entry + 6)] = self.string(strings, self.u32(entry + 8))
                step = self.u32(entry + 12)
                if step == 0 and versions_left > 1:
                    raise ValueError(f"{self.name} has a version chain that ends before its declared count")
                entry += step
            if following == 0 and libraries_left > 1:
                raise ValueError(f"{self.name} has a library chain that ends before its declared count")
            at += following
        return needs

    def symbols_at(self, wanted: set[int]) -> list[str]:
        """The dynamic symbols bound to any of the version indexes in `wanted`."""
        versym, dynsym = self.find(SHT_GNU_VERSYM), self.find(SHT_DYNSYM)
        if versym is None or dynsym is None or dynsym["entsize"] == 0:
            return []
        strings = self.sections[dynsym["link"]]
        names = []
        for i in range(dynsym["size"] // dynsym["entsize"]):
            # The high bit of a version index marks a hidden symbol and is not part of the index.
            if self.u16(versym["offset"] + 2 * i) & 0x7FFF in wanted:
                names.append(self.string(strings, self.u32(dynsym["offset"] + i * dynsym["entsize"])))
        return sorted(set(names))


def glibc_versions(elf: Elf) -> dict[str, list[str]]:
    """The GLIBC versions `elf` requires, each with the symbols drawn from it, highest first."""
    needs = elf.version_needs()
    versions: dict[str, set[int]] = {}
    for index, name in needs.items():
        if not name.startswith("GLIBC_"):
            continue
        m = GLIBC_VERSION.match(name)
        if m:
            version = m.group(1)
        elif name in GLIBC_ABI_VERSIONS:
            version = GLIBC_ABI_VERSIONS[name]
        else:
            raise ValueError(f"{elf.name} requires {name}, a glibc version this reader cannot place")
        versions.setdefault(version, set()).add(index)
    return {v: elf.symbols_at(versions[v]) for v in sorted(versions, key=version_key, reverse=True)}


def elf_members(path: Path) -> list[tuple[str, bytes]]:
    """(name, bytes) of every ELF file `path` holds: the file itself, or each ELF member of a tar archive."""
    if tarfile.is_tarfile(path):
        members = []
        with tarfile.open(path) as tar:
            for member in tar:
                if member.isfile():
                    data = tar.extractfile(member).read()
                    if data[:4] == b"\x7fELF":
                        members.append((f"{path.name}/{Path(member.name).name}", data))
        if not members:
            fail(f"{path} holds no ELF file")
        return members
    data = path.read_bytes()
    if data[:4] != b"\x7fELF":
        fail(f"{path} is neither an ELF file nor a tar archive")
    return [(path.name, data)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--floor", help="the glibc version the release page states; a binary needing more fails")
    ap.add_argument("--json", type=Path, help="write {path name: measured floor} here")
    ap.add_argument("paths", nargs="+", type=Path)
    a = ap.parse_args()
    if a.floor and not re.fullmatch(r"\d+(\.\d+)+", a.floor):
        fail(f"--floor {a.floor!r} is not a glibc version such as 2.17")

    floors: dict[str, str] = {}
    failures: list[str] = []
    for path in a.paths:
        if not path.is_file():
            fail(f"{path} is not a file")
        for name, data in elf_members(path):
            try:
                versions = glibc_versions(Elf(data, name))
            except ValueError as exc:
                fail(str(exc))
            except (IndexError, struct.error):
                fail(f"{name} is truncated or not a well-formed ELF file")
            if not versions:
                fail(f"{name} requires no GLIBC version, so its floor cannot be read")
            highest = next(iter(versions))
            print(f"{name}: glibc {highest}")
            if a.floor and version_key(highest) > version_key(a.floor):
                above = [f"GLIBC_{v} ({', '.join(syms) or 'no named symbol'})" for v, syms in versions.items()
                         if version_key(v) > version_key(a.floor)]
                failures.append(f"{name} requires {'; '.join(above)}, above the floor {a.floor}")
            current = floors.get(path.name)
            if current is None or version_key(highest) > version_key(current):
                floors[path.name] = highest
    if a.json:
        a.json.write_text(json.dumps(floors, indent=2, sort_keys=True) + "\n")
    if failures:
        fail("; ".join(failures))


if __name__ == "__main__":
    main()
