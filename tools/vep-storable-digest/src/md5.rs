// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! MD5 (RFC 1321) of a byte string, as the lower-case hexadecimal text
//! `Digest::MD5::md5_hex` returns.

const SHIFTS: [u32; 64] = [
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9,
    14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10, 15,
    21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
];

/// `floor(2^32 * abs(sin(i + 1)))` for i in 0..64.
const SINES: [u32; 64] = [
    0xd76a_a478,
    0xe8c7_b756,
    0x2420_70db,
    0xc1bd_ceee,
    0xf57c_0faf,
    0x4787_c62a,
    0xa830_4613,
    0xfd46_9501,
    0x6980_98d8,
    0x8b44_f7af,
    0xffff_5bb1,
    0x895c_d7be,
    0x6b90_1122,
    0xfd98_7193,
    0xa679_438e,
    0x49b4_0821,
    0xf61e_2562,
    0xc040_b340,
    0x265e_5a51,
    0xe9b6_c7aa,
    0xd62f_105d,
    0x0244_1453,
    0xd8a1_e681,
    0xe7d3_fbc8,
    0x21e1_cde6,
    0xc337_07d6,
    0xf4d5_0d87,
    0x455a_14ed,
    0xa9e3_e905,
    0xfcef_a3f8,
    0x676f_02d9,
    0x8d2a_4c8a,
    0xfffa_3942,
    0x8771_f681,
    0x6d9d_6122,
    0xfde5_380c,
    0xa4be_ea44,
    0x4bde_cfa9,
    0xf6bb_4b60,
    0xbebf_bc70,
    0x289b_7ec6,
    0xeaa1_27fa,
    0xd4ef_3085,
    0x0488_1d05,
    0xd9d4_d039,
    0xe6db_99e5,
    0x1fa2_7cf8,
    0xc4ac_5665,
    0xf429_2244,
    0x432a_ff97,
    0xab94_23a7,
    0xfc93_a039,
    0x655b_59c3,
    0x8f0c_cc92,
    0xffef_f47d,
    0x8584_5dd1,
    0x6fa8_7e4f,
    0xfe2c_e6e0,
    0xa301_4314,
    0x4e08_11a1,
    0xf753_7e82,
    0xbd3a_f235,
    0x2ad7_d2bb,
    0xeb86_d391,
];

fn block(state: &mut [u32; 4], chunk: &[u8]) {
    let mut words = [0u32; 16];
    for (word, bytes) in words.iter_mut().zip(chunk.chunks_exact(4)) {
        *word = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
    }
    let [mut a, mut b, mut c, mut d] = *state;
    for round in 0..64 {
        let (mix, index) = match round / 16 {
            0 => ((b & c) | (!b & d), round),
            1 => ((d & b) | (!d & c), (5 * round + 1) % 16),
            2 => (b ^ c ^ d, (3 * round + 5) % 16),
            _ => (c ^ (b | !d), (7 * round) % 16),
        };
        let sum = a
            .wrapping_add(mix)
            .wrapping_add(SINES[round])
            .wrapping_add(words[index]);
        a = d;
        d = c;
        c = b;
        b = b.wrapping_add(sum.rotate_left(SHIFTS[round]));
    }
    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
}

/// The sixteen bytes of the digest.
pub fn md5(data: &[u8]) -> [u8; 16] {
    let mut state = [0x6745_2301u32, 0xefcd_ab89, 0x98ba_dcfe, 0x1032_5476];
    let mut chunks = data.chunks_exact(64);
    for chunk in &mut chunks {
        block(&mut state, chunk);
    }
    // The tail: the bytes left, one 0x80, zeros, and the length in bits.
    let rest = chunks.remainder();
    let mut tail = [0u8; 128];
    tail[..rest.len()].copy_from_slice(rest);
    tail[rest.len()] = 0x80;
    let tail_len = if rest.len() < 56 { 64 } else { 128 };
    let bits = (data.len() as u64).wrapping_mul(8);
    tail[tail_len - 8..tail_len].copy_from_slice(&bits.to_le_bytes());
    for chunk in tail[..tail_len].chunks_exact(64) {
        block(&mut state, chunk);
    }
    let mut digest = [0u8; 16];
    for (bytes, word) in digest.chunks_exact_mut(4).zip(state) {
        bytes.copy_from_slice(&word.to_le_bytes());
    }
    digest
}

/// The digest as 32 lower-case hexadecimal digits.
pub fn md5_hex(data: &[u8]) -> [u8; 32] {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut hex = [0u8; 32];
    for (pair, byte) in hex.chunks_exact_mut(2).zip(md5(data)) {
        pair[0] = DIGITS[usize::from(byte >> 4)];
        pair[1] = DIGITS[usize::from(byte & 0x0f)];
    }
    hex
}

#[cfg(test)]
mod tests {
    use super::md5_hex;

    fn hex(data: &[u8]) -> String {
        String::from_utf8(md5_hex(data).to_vec()).expect("hexadecimal digits")
    }

    /// The test vectors of RFC 1321, appendix A.5.
    #[test]
    fn the_test_vectors_of_rfc_1321() {
        assert_eq!(hex(b""), "d41d8cd98f00b204e9800998ecf8427e");
        assert_eq!(hex(b"a"), "0cc175b9c0f1b6a831c399e269772661");
        assert_eq!(hex(b"abc"), "900150983cd24fb0d6963f7d28e17f72");
        assert_eq!(hex(b"message digest"), "f96b697d7cb7938d525a2f31aaf161d0");
        assert_eq!(
            hex(b"abcdefghijklmnopqrstuvwxyz"),
            "c3fcd3d76192e4007dfb496cca67e13b"
        );
        assert_eq!(
            hex(b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"),
            "d174ab98d277d9f5a5611c2c9f419d9f"
        );
        assert_eq!(
            hex(
                b"12345678901234567890123456789012345678901234567890123456789012345678901234567890"
            ),
            "57edf4a22be3c955ac49da2e2107b67a"
        );
    }

    /// Lengths on both sides of the padding's two block boundaries.
    #[test]
    fn lengths_around_a_block_boundary() {
        for (length, expected) in [
            (55, "ef1772b6dff9a122358552954ad0df65"),
            (56, "3b0c8ac703f828b04c6c197006d17218"),
            (63, "b06521f39153d618550606be297466d5"),
            (64, "014842d480b571495a4a0363793f7367"),
            (65, "c743a45e0d2e6a95cb859adae0248435"),
            (119, "8a7bd0732ed6a28ce75f6dabc90e1613"),
            (120, "5f61c0ccad4cac44c75ff505e1f1e537"),
            (128, "e510683b3f5ffe4093d021808bc6ff70"),
        ] {
            assert_eq!(hex(&vec![b'a'; length]), expected, "{length} bytes");
        }
    }
}
