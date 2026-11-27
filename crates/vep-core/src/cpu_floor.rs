// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The startup check that refuses, in one line, an x86_64 Linux CPU below the
//! floor a binary was compiled for, where the binary would otherwise die of an
//! illegal instruction.
//!
//! A build with `-C target-cpu` above `x86-64` (`x86-64-v3`, `native`) may run
//! instructions of that floor from its first Rust constructor on, so the check
//! is assembly alone: [`cpu_floor_check!`](crate::cpu_floor_check) puts a
//! constructor into the binary's `.init_array.00000`, which runs before every
//! other constructor of the binary, and the constructor calls the CPUID test
//! that [`missing`] also runs from Rust. A default build requires no bit, and
//! its check passes on every CPU.

/// The bit `bit` when `on`.
const fn bit(on: bool, bit: u32) -> u32 {
    if on {
        1 << bit
    } else {
        0
    }
}

/// The CPUID leaf 1 ECX bits of the features the build assumes; AVX adds
/// OSXSAVE (bit 27), which the XCR0 read needs.
pub const LEAF1_ECX: u32 = bit(cfg!(target_feature = "sse3"), 0)
    | bit(cfg!(target_feature = "ssse3"), 9)
    | bit(cfg!(target_feature = "fma"), 12)
    | bit(cfg!(target_feature = "cmpxchg16b"), 13)
    | bit(cfg!(target_feature = "sse4.1"), 19)
    | bit(cfg!(target_feature = "sse4.2"), 20)
    | bit(cfg!(target_feature = "movbe"), 22)
    | bit(cfg!(target_feature = "popcnt"), 23)
    | bit(cfg!(target_feature = "xsave"), 26)
    | bit(cfg!(target_feature = "avx"), 27)
    | bit(cfg!(target_feature = "avx"), 28)
    | bit(cfg!(target_feature = "f16c"), 29);

/// The CPUID leaf 7 subleaf 0 EBX bits of the features the build assumes.
pub const LEAF7_EBX: u32 = bit(cfg!(target_feature = "bmi1"), 3)
    | bit(cfg!(target_feature = "avx2"), 5)
    | bit(cfg!(target_feature = "bmi2"), 8)
    | bit(cfg!(target_feature = "avx512f"), 16)
    | bit(cfg!(target_feature = "avx512dq"), 17)
    | bit(cfg!(target_feature = "avx512cd"), 28)
    | bit(cfg!(target_feature = "avx512bw"), 30)
    | bit(cfg!(target_feature = "avx512vl"), 31);

/// The CPUID leaf 0x8000_0001 ECX bits of the features the build assumes.
pub const EXT1_ECX: u32 = bit(cfg!(target_feature = "lzcnt"), 5);

/// The XCR0 bits the build needs the operating system to save: SSE and AVX
/// state for AVX, the three AVX-512 states for AVX-512.
pub const XCR0: u32 = bit(cfg!(target_feature = "avx"), 1)
    | bit(cfg!(target_feature = "avx"), 2)
    | bit(cfg!(target_feature = "avx512f"), 5)
    | bit(cfg!(target_feature = "avx512f"), 6)
    | bit(cfg!(target_feature = "avx512f"), 7);

// The masks cover the x86-64-v2, v3 and v4 features; a `native` build that
// also enables later ones (GFNI, VAES, AVX-512 VNNI) can still fault.
/// The refusal after the `ERROR: [<component>] ` prefix: the floor the build
/// assumes, and the builds that run on the CPU.
pub const MESSAGE: &str = if cfg!(target_feature = "avx512f") {
    "this CPU lacks x86-64-v4 instructions (AVX-512 F, BW, CD, DQ, VL) this build needs; a build without -C target-cpu runs on any x86_64 CPU\n"
} else if cfg!(target_feature = "avx2") {
    "this CPU lacks x86-64-v3 instructions (AVX2, BMI1, BMI2, FMA, F16C, LZCNT, MOVBE) this build needs; a build without -C target-cpu runs on any x86_64 CPU, and the x86_64-unknown-linux-gnu-v2 release archive on one with SSE4.2 and POPCNT\n"
} else {
    "this CPU lacks x86-64-v2 instructions (SSE4.2, SSSE3, POPCNT, CMPXCHG16B) this build needs; a build without -C target-cpu runs on any x86_64 CPU\n"
};

/// `prefix` followed by `message`, as the bytes the constructor writes.
pub const fn join<const N: usize>(prefix: &str, message: &str) -> [u8; N] {
    let (a, b) = (prefix.as_bytes(), message.as_bytes());
    assert!(a.len() + b.len() == N);
    let mut out = [0u8; N];
    let mut i = 0;
    while i < N {
        out[i] = if i < a.len() { a[i] } else { b[i - a.len()] };
        i += 1;
    }
    out
}

// System V call: the four masks in edi, esi, edx and ecx; eax is 0 when the
// CPU has every bit, 1 otherwise. A leaf the CPU lacks reads as no bit, and
// XCR0 is read only after the leaf 1 mask, which then holds OSXSAVE, passed.
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
core::arch::global_asm!(
    ".pushsection .text.vep_core_cpu_floor_missing,\"ax\",@progbits",
    ".globl vep_core_cpu_floor_missing",
    ".hidden vep_core_cpu_floor_missing",
    ".p2align 4",
    "vep_core_cpu_floor_missing:",
    "push rbx",
    "mov r8d, edi",
    "mov r9d, esi",
    "mov r10d, edx",
    "mov r11d, ecx",
    "xor eax, eax",
    "cpuid",
    "mov esi, eax",
    "mov eax, 1",
    "xor ecx, ecx",
    "cpuid",
    "and ecx, r8d",
    "cmp ecx, r8d",
    "jne 9f",
    "xor ebx, ebx",
    "cmp esi, 7",
    "jb 4f",
    "mov eax, 7",
    "xor ecx, ecx",
    "cpuid",
    "4:",
    "and ebx, r9d",
    "cmp ebx, r9d",
    "jne 9f",
    "mov eax, 0x80000000",
    "cpuid",
    "mov esi, eax",
    "xor ecx, ecx",
    "cmp esi, 0x80000001",
    "jb 2f",
    "mov eax, 0x80000001",
    "cpuid",
    "2:",
    "and ecx, r10d",
    "cmp ecx, r10d",
    "jne 9f",
    "test r11d, r11d",
    "jz 3f",
    "xor ecx, ecx",
    "xgetbv",
    "and eax, r11d",
    "cmp eax, r11d",
    "jne 9f",
    "3:",
    "xor eax, eax",
    "pop rbx",
    "ret",
    "9:",
    "mov eax, 1",
    "pop rbx",
    "ret",
    ".popsection",
);

#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
unsafe extern "sysv64" {
    fn vep_core_cpu_floor_missing(leaf1_ecx: u32, leaf7_ebx: u32, ext1_ecx: u32, xcr0: u32) -> u32;
}

/// Whether this CPU lacks a bit of the masks, read as [`LEAF1_ECX`],
/// [`LEAF7_EBX`], [`EXT1_ECX`] and [`XCR0`] are; an XCR0 mask adds OSXSAVE.
#[cfg(all(target_arch = "x86_64", target_os = "linux"))]
pub fn missing(leaf1_ecx: u32, leaf7_ebx: u32, ext1_ecx: u32, xcr0: u32) -> bool {
    let leaf1_ecx = leaf1_ecx | bit(xcr0 != 0, 27);
    // SAFETY: the routine reads CPUID, and XCR0 only once OSXSAVE is present;
    // it keeps rbx and writes no memory.
    unsafe { vep_core_cpu_floor_missing(leaf1_ecx, leaf7_ebx, ext1_ecx, xcr0) != 0 }
}

/// Puts the startup check into the binary that names it, with `component` in
/// its `ERROR: [<component>]` line; the refusal exits with status 1 before
/// `main`. Elsewhere than x86_64 Linux it adds nothing.
#[doc(hidden)]
#[macro_export]
macro_rules! cpu_floor_check {
    ($component:literal) => {
        #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
        mod vep_cpu_floor_check {
            const PREFIX: &str = concat!("ERROR: [", $component, "] ");
            #[used]
            static MESSAGE: [u8; PREFIX.len() + $crate::cpu_floor::MESSAGE.len()] =
                $crate::cpu_floor::join(PREFIX, $crate::cpu_floor::MESSAGE);
            // write(2, MESSAGE) and exit_group(1) by system call when the test fails.
            ::core::arch::global_asm!(
                ".pushsection .text.vep_cpu_floor_check,\"ax\",@progbits",
                ".p2align 4",
                "vep_cpu_floor_check:",
                "sub rsp, 8",
                "mov edi, {leaf1}",
                "mov esi, {leaf7}",
                "mov edx, {ext1}",
                "mov ecx, {xcr0}",
                "call vep_core_cpu_floor_missing",
                "add rsp, 8",
                "test eax, eax",
                "jnz 5f",
                "ret",
                "5:",
                "mov eax, 1",
                "mov edi, 2",
                "lea rsi, [rip + {message}]",
                "mov edx, {len}",
                "syscall",
                "mov eax, 231",
                "mov edi, 1",
                "syscall",
                ".popsection",
                ".pushsection .init_array.00000,\"aw\"",
                ".p2align 3",
                ".quad vep_cpu_floor_check",
                ".popsection",
                leaf1 = const $crate::cpu_floor::LEAF1_ECX,
                leaf7 = const $crate::cpu_floor::LEAF7_EBX,
                ext1 = const $crate::cpu_floor::EXT1_ECX,
                xcr0 = const $crate::cpu_floor::XCR0,
                message = sym MESSAGE,
                len = const PREFIX.len() + $crate::cpu_floor::MESSAGE.len(),
            );
        }
    };
}

#[cfg(all(test, target_arch = "x86_64", target_os = "linux"))]
mod tests {
    use super::*;

    #[test]
    fn the_build_floor_and_no_floor_pass() {
        assert!(!missing(0, 0, 0, 0));
        assert!(!missing(LEAF1_ECX, LEAF7_EBX, EXT1_ECX, XCR0));
    }

    #[test]
    fn a_feature_no_cpu_has_is_missing() {
        assert!(missing(u32::MAX, 0, 0, 0));
        assert!(missing(0, u32::MAX, 0, 0));
        assert!(missing(0, 0, u32::MAX, 0));
        assert!(missing(0, 0, 0, u32::MAX));
    }

    #[test]
    fn the_refusal_is_one_error_line() {
        const PREFIX: &str = "ERROR: [vep] ";
        let line: [u8; PREFIX.len() + MESSAGE.len()] = join(PREFIX, MESSAGE);
        let line = std::str::from_utf8(&line).unwrap();
        assert!(line.starts_with("ERROR: [vep] this CPU lacks x86-64-v"));
        assert_eq!(line.matches('\n').count(), 1);
        assert!(line.ends_with('\n'));
    }
}
