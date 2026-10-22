// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! SHA256 (FIPS 180-4, section 6.2), the digest the cache checksum manifest
//! records, computed incrementally over a byte stream.

/// The initial hash value (FIPS 180-4, section 5.3.3).
const H0: [u32; 8] = [
    0x6a09_e667,
    0xbb67_ae85,
    0x3c6e_f372,
    0xa54f_f53a,
    0x510e_527f,
    0x9b05_688c,
    0x1f83_d9ab,
    0x5be0_cd19,
];

/// The round constants (FIPS 180-4, section 4.2.2).
const K: [u32; 64] = [
    0x428a_2f98,
    0x7137_4491,
    0xb5c0_fbcf,
    0xe9b5_dba5,
    0x3956_c25b,
    0x59f1_11f1,
    0x923f_82a4,
    0xab1c_5ed5,
    0xd807_aa98,
    0x1283_5b01,
    0x2431_85be,
    0x550c_7dc3,
    0x72be_5d74,
    0x80de_b1fe,
    0x9bdc_06a7,
    0xc19b_f174,
    0xe49b_69c1,
    0xefbe_4786,
    0x0fc1_9dc6,
    0x240c_a1cc,
    0x2de9_2c6f,
    0x4a74_84aa,
    0x5cb0_a9dc,
    0x76f9_88da,
    0x983e_5152,
    0xa831_c66d,
    0xb003_27c8,
    0xbf59_7fc7,
    0xc6e0_0bf3,
    0xd5a7_9147,
    0x06ca_6351,
    0x1429_2967,
    0x27b7_0a85,
    0x2e1b_2138,
    0x4d2c_6dfc,
    0x5338_0d13,
    0x650a_7354,
    0x766a_0abb,
    0x81c2_c92e,
    0x9272_2c85,
    0xa2bf_e8a1,
    0xa81a_664b,
    0xc24b_8b70,
    0xc76c_51a3,
    0xd192_e819,
    0xd699_0624,
    0xf40e_3585,
    0x106a_a070,
    0x19a4_c116,
    0x1e37_6c08,
    0x2748_774c,
    0x34b0_bcb5,
    0x391c_0cb3,
    0x4ed8_aa4a,
    0x5b9c_ca4f,
    0x682e_6ff3,
    0x748f_82ee,
    0x78a5_636f,
    0x84c8_7814,
    0x8cc7_0208,
    0x90be_fffa,
    0xa450_6ceb,
    0xbef9_a3f7,
    0xc671_78f2,
];

const BLOCK_LEN: usize = 64;

/// A 32-byte SHA256 digest.
pub type Digest = [u8; 32];

/// A SHA256 computation in progress: bytes go in through [`Sha256::update`]
/// in any number of calls and the digest comes out of [`Sha256::finalize`].
#[derive(Clone)]
pub struct Sha256 {
    state: [u32; 8],
    block: [u8; BLOCK_LEN],
    /// Bytes of `block` filled by the current partial block.
    pending: usize,
    /// Bytes hashed so far, the message length the padding records.
    length: u64,
}

impl Default for Sha256 {
    fn default() -> Self {
        Self::new()
    }
}

impl Sha256 {
    /// A computation over no bytes yet.
    pub fn new() -> Self {
        Sha256 {
            state: H0,
            block: [0; BLOCK_LEN],
            pending: 0,
            length: 0,
        }
    }

    /// The digest of `data` in one call.
    pub fn digest(data: &[u8]) -> Digest {
        let mut hasher = Sha256::new();
        hasher.update(data);
        hasher.finalize()
    }

    /// Appends `data` to the message.
    pub fn update(&mut self, mut data: &[u8]) {
        self.length = self.length.wrapping_add(data.len() as u64);
        if self.pending > 0 {
            let take = (BLOCK_LEN - self.pending).min(data.len());
            self.block[self.pending..self.pending + take].copy_from_slice(&data[..take]);
            self.pending += take;
            data = &data[take..];
            if self.pending < BLOCK_LEN {
                return;
            }
            compress(&mut self.state, &self.block);
            self.pending = 0;
        }
        let mut blocks = data.chunks_exact(BLOCK_LEN);
        for block in &mut blocks {
            compress(&mut self.state, block.try_into().expect("a 64-byte chunk"));
        }
        let rest = blocks.remainder();
        self.block[..rest.len()].copy_from_slice(rest);
        self.pending = rest.len();
    }

    /// Pads the message (FIPS 180-4, section 5.1.1) and returns its digest.
    pub fn finalize(mut self) -> Digest {
        let bit_length = self.length.wrapping_mul(8);
        self.block[self.pending] = 0x80;
        self.pending += 1;
        if self.pending > BLOCK_LEN - 8 {
            self.block[self.pending..].fill(0);
            compress(&mut self.state, &self.block);
            self.pending = 0;
        }
        self.block[self.pending..BLOCK_LEN - 8].fill(0);
        self.block[BLOCK_LEN - 8..].copy_from_slice(&bit_length.to_be_bytes());
        compress(&mut self.state, &self.block);
        let mut digest = [0u8; 32];
        for (chunk, word) in digest.chunks_exact_mut(4).zip(self.state) {
            chunk.copy_from_slice(&word.to_be_bytes());
        }
        digest
    }
}

/// One compression of a 64-byte block into `state` (FIPS 180-4, section 6.2.2).
fn compress(state: &mut [u32; 8], block: &[u8; BLOCK_LEN]) {
    let mut w = [0u32; 64];
    for (t, chunk) in block.chunks_exact(4).enumerate() {
        w[t] = u32::from_be_bytes(chunk.try_into().expect("a 4-byte chunk"));
    }
    for t in 16..64 {
        let s0 = w[t - 15].rotate_right(7) ^ w[t - 15].rotate_right(18) ^ (w[t - 15] >> 3);
        let s1 = w[t - 2].rotate_right(17) ^ w[t - 2].rotate_right(19) ^ (w[t - 2] >> 10);
        w[t] = w[t - 16]
            .wrapping_add(s0)
            .wrapping_add(w[t - 7])
            .wrapping_add(s1);
    }
    let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = *state;
    for t in 0..64 {
        let big_s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
        let ch = (e & f) ^ (!e & g);
        let t1 = h
            .wrapping_add(big_s1)
            .wrapping_add(ch)
            .wrapping_add(K[t])
            .wrapping_add(w[t]);
        let big_s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
        let maj = (a & b) ^ (a & c) ^ (b & c);
        let t2 = big_s0.wrapping_add(maj);
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(t1);
        d = c;
        c = b;
        b = a;
        a = t1.wrapping_add(t2);
    }
    for (word, value) in state.iter_mut().zip([a, b, c, d, e, f, g, h]) {
        *word = word.wrapping_add(value);
    }
}

/// The lower-case hexadecimal spelling of a digest, as `sha256sum` prints it.
pub fn to_hex(digest: &Digest) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut text = String::with_capacity(64);
    for byte in digest {
        text.push(HEX[usize::from(byte >> 4)] as char);
        text.push(HEX[usize::from(byte & 0x0f)] as char);
    }
    text
}

/// The digest a 64-character hexadecimal spelling denotes, in either case;
/// `None` for any other text.
pub fn from_hex(text: &str) -> Option<Digest> {
    let bytes = text.as_bytes();
    if bytes.len() != 64 {
        return None;
    }
    let mut digest = [0u8; 32];
    for (byte, pair) in digest.iter_mut().zip(bytes.chunks_exact(2)) {
        let high = (pair[0] as char).to_digit(16)?;
        let low = (pair[1] as char).to_digit(16)?;
        *byte = (high as u8) << 4 | low as u8;
    }
    Some(digest)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The digests FIPS 180-4's example vectors and `sha256sum` give.
    #[test]
    fn known_vectors() {
        assert_eq!(
            to_hex(&Sha256::digest(b"")),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        assert_eq!(
            to_hex(&Sha256::digest(b"abc")),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        assert_eq!(
            to_hex(&Sha256::digest(
                b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"
            )),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
        assert_eq!(
            to_hex(&Sha256::digest(&[b'a'; 1_000_000])),
            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
        );
    }

    /// Messages whose padding crosses a block boundary: 55 bytes fit the
    /// length in their own block, 56 to 64 need a second one.
    #[test]
    fn padding_boundaries_match_one_shot_digests() {
        for len in [55usize, 56, 63, 64, 65, 119, 120, 127, 128] {
            let data: Vec<u8> = (0..len).map(|i| i as u8).collect();
            let mut split = Sha256::new();
            split.update(&data[..len / 3]);
            split.update(&data[len / 3..]);
            assert_eq!(split.finalize(), Sha256::digest(&data), "length {len}");
        }
    }

    /// Byte-by-byte updates and one update give the same digest.
    #[test]
    fn incremental_updates_match_one_shot_digest() {
        let data: Vec<u8> = (0..1000u32).map(|i| (i * 7 % 251) as u8).collect();
        let mut hasher = Sha256::new();
        for byte in &data {
            hasher.update(std::slice::from_ref(byte));
        }
        assert_eq!(hasher.finalize(), Sha256::digest(&data));
    }

    #[test]
    fn hex_round_trips() {
        let digest = Sha256::digest(b"abc");
        assert_eq!(from_hex(&to_hex(&digest)), Some(digest));
        assert_eq!(from_hex(&to_hex(&digest).to_uppercase()), Some(digest));
        assert_eq!(from_hex("abc"), None);
        assert_eq!(from_hex(&"g".repeat(64)), None);
    }
}
