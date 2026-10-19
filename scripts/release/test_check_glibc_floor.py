"""Tests for check_glibc_floor.py on ELF files the tests assemble: the `.gnu.version_r`, `.gnu.version`
and `.dynsym` sections of a dynamically linked executable, in each ELF class and byte order."""

from __future__ import annotations

import io
import json
import struct
import subprocess
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "check_glibc_floor.py"
SHT_STRTAB, SHT_DYNSYM, SHT_GNU_VERNEED, SHT_GNU_VERSYM = 3, 11, 0x6FFFFFFE, 0x6FFFFFFF


def make_elf(needs: dict[str, dict[str, list[str]]], bits: int = 64, big_endian: bool = False) -> bytes:
    """An ELF executable whose version needs are `needs`: {library: {version name: [symbols bound to it]}}.
    Sections: null, .dynstr, .dynsym, .gnu.version, .gnu.version_r; no program headers."""
    o = ">" if big_endian else "<"
    strings = bytearray(b"\0")

    def intern(s: str) -> int:
        at = len(strings)
        strings.extend(s.encode() + b"\0")
        return at

    # Version indexes start at 2 (0 and 1 mean local and global).
    index = 2
    symbols: list[tuple[int, int]] = []  # (name offset, version index)
    verneed = bytearray()
    libraries = list(needs.items())
    for n, (library, versions) in enumerate(libraries):
        file_at = intern(library)
        aux = bytearray()
        items = list(versions.items())
        for m, (version, syms) in enumerate(items):
            name_at = intern(version)
            for sym in syms:
                symbols.append((intern(sym), index))
            next_aux = 16 if m + 1 < len(items) else 0
            aux += struct.pack(f"{o}IHHII", 0, 0, index, name_at, next_aux)
            index += 1
        entry = struct.pack(f"{o}HHIII", 1, len(items), file_at, 16, 16 + len(aux) if n + 1 < len(libraries) else 0)
        verneed += entry + aux
    sym_size = 24 if bits == 64 else 16
    dynsym = bytearray(b"\0" * sym_size)  # the null symbol
    versym = bytearray(struct.pack(f"{o}H", 0))
    for name_at, version_index in symbols:
        if bits == 64:
            dynsym += struct.pack(f"{o}IBBHQQ", name_at, 0x12, 0, 0, 0, 0)
        else:
            dynsym += struct.pack(f"{o}IIIBBH", name_at, 0, 0, 0x12, 0, 0)
        versym += struct.pack(f"{o}H", version_index)

    ehsize, shentsize = (64, 64) if bits == 64 else (52, 40)
    bodies = [bytes(strings), bytes(dynsym), bytes(versym), bytes(verneed)]
    offsets, at = [], ehsize
    for body in bodies:
        offsets.append(at)
        at += len(body)
    shoff = at

    def section(kind: int, offset: int, size: int, link: int, info: int, entsize: int) -> bytes:
        if bits == 64:
            return struct.pack(f"{o}IIQQQQIIQQ", 0, kind, 0, 0, offset, size, link, info, 0, entsize)
        return struct.pack(f"{o}IIIIIIIIII", 0, kind, 0, 0, offset, size, link, info, 0, entsize)

    headers = b"\0" * shentsize
    headers += section(SHT_STRTAB, offsets[0], len(bodies[0]), 0, 0, 0)
    headers += section(SHT_DYNSYM, offsets[1], len(bodies[1]), 1, 0, sym_size)
    headers += section(SHT_GNU_VERSYM, offsets[2], len(bodies[2]), 2, 0, 2)
    headers += section(SHT_GNU_VERNEED, offsets[3], len(bodies[3]), 1, len(libraries), 0)
    ident = b"\x7fELF" + bytes([2 if bits == 64 else 1, 2 if big_endian else 1, 1, 0]) + b"\0" * 8
    if bits == 64:
        header = ident + struct.pack(f"{o}HHIQQQIHHHHHH", 2, 0, 1, 0, 0, shoff, 0, ehsize, 0, 0, shentsize, 5, 0)
    else:
        header = ident + struct.pack(f"{o}HHIIIIIHHHHHH", 2, 0, 1, 0, 0, shoff, 0, ehsize, 0, 0, shentsize, 5, 0)
    return header + b"".join(bodies) + headers


OLD_VEP = {"libc.so.6": {"GLIBC_2.2.5": ["malloc", "read"], "GLIBC_2.34": ["pthread_create"], "GLIBC_2.39": ["pidfd_spawnp", "pidfd_getpid"]}}
FLOOR_VEP = {"libc.so.6": {"GLIBC_2.17": ["clock_gettime"], "GLIBC_2.2.5": ["malloc"]}, "libssl.so.3": {"OPENSSL_3.0.0": ["SSL_new"]}}
AARCH64_VEP = {"libc.so.6": {"GLIBC_2.17": ["malloc", "read", "pthread_create"]}}


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run([sys.executable, str(SCRIPT), *args], check=False, capture_output=True, text=True)


def write_archive(path: Path, members: dict[str, bytes]) -> None:
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


def test_each_binary_prints_its_highest_required_glibc_version(tmp_path: Path) -> None:
    (tmp_path / "vep").write_bytes(make_elf(OLD_VEP))
    (tmp_path / "builder").write_bytes(make_elf(FLOOR_VEP))
    out = run(str(tmp_path / "vep"), str(tmp_path / "builder"))
    assert out.returncode == 0, out.stderr
    assert out.stdout == "vep: glibc 2.39\nbuilder: glibc 2.17\n"


def test_the_packed_relocation_marker_counts_as_glibc_2_36_and_an_unknown_glibc_name_is_refused(tmp_path: Path) -> None:
    (tmp_path / "relr").write_bytes(make_elf({"libc.so.6": {"GLIBC_2.17": ["malloc"], "GLIBC_ABI_DT_RELR": []}}))
    assert run(str(tmp_path / "relr")).stdout == "relr: glibc 2.36\n"
    out = run("--floor", "2.17", str(tmp_path / "relr"))
    assert out.returncode == 1 and "relr requires GLIBC_2.36 (no named symbol), above the floor 2.17" in out.stderr
    (tmp_path / "odd").write_bytes(make_elf({"libc.so.6": {"GLIBC_2.17": ["malloc"], "GLIBC_PRIVATE": ["_dl_x"]}}))
    out = run(str(tmp_path / "odd"))
    assert out.returncode == 1 and out.stderr.strip() == (
        "ERROR: [check_glibc_floor] odd requires GLIBC_PRIVATE, a glibc version this reader cannot place"
    )


def test_versions_compare_numerically_not_as_text(tmp_path: Path) -> None:
    (tmp_path / "bin").write_bytes(make_elf({"libc.so.6": {"GLIBC_2.9": ["pipe2"], "GLIBC_2.14": ["memcpy"], "GLIBC_2.2.5": ["read"]}}))
    assert run(str(tmp_path / "bin")).stdout == "bin: glibc 2.14\n"
    assert run("--floor", "2.17", str(tmp_path / "bin")).returncode == 0


def test_the_floor_fails_a_binary_above_it_and_names_the_symbols_at_each_version(tmp_path: Path) -> None:
    (tmp_path / "vep").write_bytes(make_elf(OLD_VEP))
    out = run("--floor", "2.17", str(tmp_path / "vep"))
    assert out.returncode == 1
    assert out.stdout == "vep: glibc 2.39\n"
    assert out.stderr.strip() == (
        "ERROR: [check_glibc_floor] vep requires GLIBC_2.39 (pidfd_getpid, pidfd_spawnp); "
        "GLIBC_2.34 (pthread_create), above the floor 2.17"
    )


def test_a_binary_at_or_below_the_floor_passes_and_other_libraries_versions_do_not_count(tmp_path: Path) -> None:
    (tmp_path / "builder").write_bytes(make_elf(FLOOR_VEP))
    assert run("--floor", "2.17", str(tmp_path / "builder")).returncode == 0
    assert run("--floor", "2.17.1", str(tmp_path / "builder")).returncode == 0
    out = run("--floor", "2.16", str(tmp_path / "builder"))
    assert out.returncode == 1 and "builder requires GLIBC_2.17 (clock_gettime), above the floor 2.16" in out.stderr


def test_an_archive_reads_every_elf_member_and_the_json_floor_is_the_highest_of_them(tmp_path: Path) -> None:
    archive = tmp_path / "vep-9.9.9-x86_64-unknown-linux-gnu.tar.gz"
    write_archive(archive, {
        "vep-9.9.9-x86_64-unknown-linux-gnu/vep": make_elf(FLOOR_VEP),
        "vep-9.9.9-x86_64-unknown-linux-gnu/vep-cache-builder": make_elf(OLD_VEP),
        "vep-9.9.9-x86_64-unknown-linux-gnu/README": b"not an ELF file\n",
    })
    arm = tmp_path / "vep-9.9.9-aarch64-unknown-linux-gnu.tar.gz"
    write_archive(arm, {"vep-9.9.9-aarch64-unknown-linux-gnu/vep": make_elf(AARCH64_VEP)})
    floors = tmp_path / "floors.json"
    out = run("--json", str(floors), str(archive), str(arm))
    assert out.returncode == 0, out.stderr
    assert out.stdout == (
        "vep-9.9.9-x86_64-unknown-linux-gnu.tar.gz/vep: glibc 2.17\n"
        "vep-9.9.9-x86_64-unknown-linux-gnu.tar.gz/vep-cache-builder: glibc 2.39\n"
        "vep-9.9.9-aarch64-unknown-linux-gnu.tar.gz/vep: glibc 2.17\n"
    )
    assert json.loads(floors.read_text()) == {
        "vep-9.9.9-aarch64-unknown-linux-gnu.tar.gz": "2.17",
        "vep-9.9.9-x86_64-unknown-linux-gnu.tar.gz": "2.39",
    }
    out = run("--floor", "2.17", "--json", str(floors), str(archive), str(arm))
    assert out.returncode == 1 and "vep-cache-builder requires GLIBC_2.39" in out.stderr
    assert json.loads(floors.read_text())["vep-9.9.9-x86_64-unknown-linux-gnu.tar.gz"] == "2.39"


def test_every_elf_class_and_byte_order_is_read(tmp_path: Path) -> None:
    for bits in (32, 64):
        for big_endian in (False, True):
            path = tmp_path / f"bin{bits}{'be' if big_endian else 'le'}"
            path.write_bytes(make_elf(OLD_VEP, bits=bits, big_endian=big_endian))
            out = run(str(path))
            assert out.returncode == 0, out.stderr
            assert out.stdout == f"{path.name}: glibc 2.39\n"


def test_inputs_the_floor_cannot_be_read_from_are_refused_by_name(tmp_path: Path) -> None:
    (tmp_path / "static").write_bytes(make_elf({"libm.so.6": {"LIBM_1.0": ["sin"]}}))
    out = run(str(tmp_path / "static"))
    assert out.returncode == 1 and out.stderr.strip() == (
        "ERROR: [check_glibc_floor] static requires no GLIBC version, so its floor cannot be read"
    )
    (tmp_path / "notes.txt").write_text("hello\n")
    out = run(str(tmp_path / "notes.txt"))
    assert out.returncode == 1 and "notes.txt is neither an ELF file nor a tar archive" in out.stderr
    write_archive(tmp_path / "docs.tar.gz", {"README": b"text\n"})
    out = run(str(tmp_path / "docs.tar.gz"))
    assert out.returncode == 1 and "docs.tar.gz holds no ELF file" in out.stderr
    out = run(str(tmp_path / "absent"))
    assert out.returncode == 1 and "absent is not a file" in out.stderr
    out = run(str(tmp_path))
    assert out.returncode == 1 and f"{tmp_path} is not a file" in out.stderr
    out = run("--floor", "two", str(tmp_path / "static"))
    assert out.returncode == 1 and "--floor 'two' is not a glibc version such as 2.17" in out.stderr
    (tmp_path / "odd").write_bytes(b"\x7fELF\x03\x01" + b"\0" * 58)
    out = run(str(tmp_path / "odd"))
    assert out.returncode == 1 and out.stderr.strip() == (
        "ERROR: [check_glibc_floor] odd has an ELF class or byte order this reader does not know"
    )
    (tmp_path / "short").write_bytes(b"\x7fELF\x02\x01" + b"\0" * 10)
    out = run(str(tmp_path / "short"))
    assert out.returncode == 1 and out.stderr.strip() == (
        "ERROR: [check_glibc_floor] short is truncated or not a well-formed ELF file"
    )
    (tmp_path / "cut").write_bytes(make_elf(OLD_VEP)[:-40])
    out = run(str(tmp_path / "cut"))
    assert out.returncode == 1 and "cut is truncated or not a well-formed ELF file" in out.stderr
    # e_shoff sits at 0x28 of a 64-bit ELF header; the .dynstr header is the second section header (after the
    # null one), and its sh_offset sits 24 bytes in.
    elf = bytearray(make_elf(OLD_VEP))
    shoff = struct.unpack_from("<Q", elf, 0x28)[0]
    struct.pack_into("<Q", elf, shoff + 64 + 24, len(elf) + 4096)
    (tmp_path / "strtab").write_bytes(elf)
    out = run(str(tmp_path / "strtab"))
    assert out.returncode == 1 and out.stderr.strip() == (
        "ERROR: [check_glibc_floor] strtab has a string table entry past the end of the file"
    )
    # A chain whose first Elf_Vernaux says "no next" while two more are declared would otherwise read the first
    # entry three times and drop the newer versions; the first Vernaux's vna_next is the last field of its 16 bytes.
    elf = bytearray(make_elf(OLD_VEP))
    verneed = elf.rindex(b"\x01\x00\x03\x00")
    struct.pack_into("<I", elf, verneed + 16 + 12, 0)
    (tmp_path / "chain").write_bytes(elf)
    out = run(str(tmp_path / "chain"))
    assert out.returncode == 1 and out.stderr.strip() == (
        "ERROR: [check_glibc_floor] chain has a version chain that ends before its declared count"
    )
    # An Elf_Verneed whose vn_next says "no next" while a second library is declared would otherwise drop that
    # library's versions and pass the floor; vn_next is the last field of the 16-byte Elf_Verneed that opens
    # .gnu.version_r, the fifth section, whose header's sh_offset sits 24 bytes in.
    elf = bytearray(make_elf({"libc.so.6": {"GLIBC_2.17": ["malloc"]}, "libm.so.6": {"GLIBC_2.29": ["exp"]}}))
    shoff = struct.unpack_from("<Q", elf, 0x28)[0]
    verneed = struct.unpack_from("<Q", elf, shoff + 4 * 64 + 24)[0]
    struct.pack_into("<I", elf, verneed + 12, 0)
    (tmp_path / "libraries").write_bytes(elf)
    out = run("--floor", "2.17", str(tmp_path / "libraries"))
    assert out.returncode == 1 and out.stderr.strip() == (
        "ERROR: [check_glibc_floor] libraries has a library chain that ends before its declared count"
    )
