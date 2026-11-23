// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The generator's random stream: SplitMix64 and a keyed hash, both written here so that
//! no dependency's release can change the records a cache and a seed give.

/// The SplitMix64 generator: a 64-bit state stepped by the golden-ratio increment and
/// finished by two multiply-and-shift rounds.
#[derive(Clone, Debug)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    /// A stream that starts from `state`.
    pub fn new(state: u64) -> SplitMix64 {
        SplitMix64 { state }
    }

    /// The stream's next 64 bits.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        finish(self.state)
    }

    /// A number below `bound`, every value equally likely: a draw from the uneven tail of
    /// the 64-bit range is thrown away and drawn again. `bound` is above zero.
    pub fn below(&mut self, bound: u64) -> u64 {
        assert!(bound > 0, "a draw needs a bound above zero");
        let tail = (u64::MAX - bound + 1) % bound;
        loop {
            let draw = self.next_u64();
            if draw <= u64::MAX - tail {
                return draw % bound;
            }
        }
    }

    /// A number from `low` to `high`, both included. `low` is at most `high`.
    pub fn between(&mut self, low: u64, high: u64) -> u64 {
        assert!(low <= high, "a draw needs a range");
        match (high - low).checked_add(1) {
            Some(width) => low + self.below(width),
            None => self.next_u64(),
        }
    }
}

/// SplitMix64's output function.
fn finish(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The hash of a seed and a list of byte strings: FNV-1a over the seed's eight
/// little-endian bytes and then over each part with its length before it, finished by
/// SplitMix64's output function. The length keeps `["ab", "c"]` apart from `["a", "bc"]`.
pub fn keyed_hash(seed: u64, parts: &[&[u8]]) -> u64 {
    const OFFSET: u64 = 0xCBF2_9CE4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01B3;
    let mut hash = OFFSET;
    let mut absorb = |bytes: &[u8]| {
        for &byte in bytes {
            hash ^= u64::from(byte);
            hash = hash.wrapping_mul(PRIME);
        }
    };
    absorb(&seed.to_le_bytes());
    for part in parts {
        absorb(&(part.len() as u64).to_le_bytes());
        absorb(part);
    }
    finish(hash)
}

/// The stream a seed and a list of byte strings name: equal names give equal streams on
/// every machine, and the order in which other streams are drawn from changes nothing.
pub fn keyed_stream(seed: u64, parts: &[&[u8]]) -> SplitMix64 {
    SplitMix64::new(keyed_hash(seed, parts))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The first outputs of SplitMix64 from state 0, as its reference implementation
    /// gives them.
    #[test]
    fn the_stream_is_splitmix64() {
        let mut stream = SplitMix64::new(0);
        assert_eq!(stream.next_u64(), 0xE220_A839_7B1D_CDAF);
        assert_eq!(stream.next_u64(), 0x6E78_9E6A_A1B9_65F4);
        assert_eq!(stream.next_u64(), 0x06C4_5D18_8009_454F);
    }

    #[test]
    fn a_bounded_draw_stays_below_its_bound() {
        let mut stream = SplitMix64::new(7);
        for bound in [1, 2, 3, 7, 64, 5000, u64::MAX] {
            for _ in 0..200 {
                assert!(stream.below(bound) < bound);
            }
        }
        for _ in 0..200 {
            let draw = stream.between(10, 12);
            assert!((10..=12).contains(&draw));
        }
        assert_eq!(stream.between(4, 4), 4);
    }

    #[test]
    fn the_keyed_hash_tells_seeds_and_part_boundaries_apart() {
        let one = keyed_hash(1, &[b"ab", b"c"]);
        assert_eq!(one, keyed_hash(1, &[b"ab", b"c"]));
        assert_ne!(one, keyed_hash(2, &[b"ab", b"c"]));
        assert_ne!(one, keyed_hash(1, &[b"a", b"bc"]));
        assert_ne!(one, keyed_hash(1, &[b"abc"]));
    }

    /// Pins the hash itself: a change of these bytes changes every genome's records.
    #[test]
    fn the_keyed_hash_is_pinned() {
        assert_eq!(keyed_hash(1, &[b"missense_forward", b"19"]), PINNED);
    }

    const PINNED: u64 = 0x48AD_BC65_2056_40AE;
}
