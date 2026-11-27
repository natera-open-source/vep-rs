// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The startup check in a binary of its own: this test binary carries the
//! constructor `vep_core::cpu_floor_check!` adds, which returns on the CPU the
//! binary was built for.

vep_core::cpu_floor_check!("cpu-floor-test");

#[test]
fn the_check_passes_on_the_cpu_the_binary_was_built_for() {
    // The constructor ran before the harness: a refusal would have ended the
    // process with status 1 before this test.
    #[cfg(all(target_arch = "x86_64", target_os = "linux"))]
    assert!(!vep_core::cpu_floor::missing(
        vep_core::cpu_floor::LEAF1_ECX,
        vep_core::cpu_floor::LEAF7_EBX,
        vep_core::cpu_floor::EXT1_ECX,
        vep_core::cpu_floor::XCR0,
    ));
}
