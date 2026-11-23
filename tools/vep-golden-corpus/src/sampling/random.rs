// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! CPython's `random.Random(seed)` for an integer seed: the Mersenne Twister MT19937 seeded
//! as `init_by_array` seeds it from the seed's magnitude, and the draws `randrange` and
//! `shuffle` make from it through `getrandbits`.

const N: usize = 624;
const M: usize = 397;
const UPPER: u32 = 0x8000_0000;
const LOWER: u32 = 0x7fff_ffff;

/// The generator `random.Random(seed)` returns for an integer seed.
#[derive(Clone)]
pub struct Random {
    state: [u32; N],
    /// The place of the next word in `state`; `N` when every word is taken.
    index: usize,
}

impl Random {
    /// `random.Random(seed)`: the generator seeded from `abs(seed)`, whose 32-bit digits
    /// from the least significant up are the key of `init_by_array`; a seed of zero is the
    /// key of one zero.
    pub fn new(seed: i64) -> Random {
        let magnitude = seed.unsigned_abs();
        let low = magnitude as u32;
        let high = (magnitude >> 32) as u32;
        if high == 0 {
            Random::from_key(&[low])
        } else {
            Random::from_key(&[low, high])
        }
    }

    /// `init_genrand(19650218)` and `init_by_array(key)` of the reference implementation.
    fn from_key(key: &[u32]) -> Random {
        let mut state = [0u32; N];
        state[0] = 19_650_218;
        for i in 1..N {
            let previous = state[i - 1];
            state[i] = 1_812_433_253u32
                .wrapping_mul(previous ^ (previous >> 30))
                .wrapping_add(i as u32);
        }
        let (mut i, mut j) = (1, 0);
        for _ in 0..N.max(key.len()) {
            let previous = state[i - 1];
            state[i] = (state[i] ^ (previous ^ (previous >> 30)).wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                state[0] = state[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..N - 1 {
            let previous = state[i - 1];
            state[i] = (state[i] ^ (previous ^ (previous >> 30)).wrapping_mul(1_566_083_941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                state[0] = state[N - 1];
                i = 1;
            }
        }
        state[0] = UPPER;
        Random { state, index: N }
    }

    /// The next 32-bit word (`genrand_uint32`).
    fn word(&mut self) -> u32 {
        if self.index >= N {
            let twist = |upper: u32, lower: u32, far: u32| {
                let mixed = (upper & UPPER) | (lower & LOWER);
                far ^ (mixed >> 1) ^ if mixed & 1 == 0 { 0 } else { 0x9908_b0df }
            };
            for k in 0..N - M {
                self.state[k] = twist(self.state[k], self.state[k + 1], self.state[k + M]);
            }
            for k in N - M..N - 1 {
                self.state[k] = twist(self.state[k], self.state[k + 1], self.state[k + M - N]);
            }
            self.state[N - 1] = twist(self.state[N - 1], self.state[0], self.state[M - 1]);
            self.index = 0;
        }
        let mut word = self.state[self.index];
        self.index += 1;
        word ^= word >> 11;
        word ^= (word << 7) & 0x9d2c_5680;
        word ^= (word << 15) & 0xefc6_0000;
        word ^ (word >> 18)
    }

    /// `getrandbits(bits)` for 1 to 64 bits: the top bits of one word, or for more than 32
    /// bits a whole word as the low half and the top bits of a second word above it.
    pub fn getrandbits(&mut self, bits: u32) -> u64 {
        assert!((1..=64).contains(&bits), "getrandbits of 1 to 64 bits");
        if bits <= 32 {
            return u64::from(self.word() >> (32 - bits));
        }
        let low = u64::from(self.word());
        let high = u64::from(self.word() >> (64 - bits));
        low | (high << 32)
    }

    /// `randrange(n)` for a count of one or more: a whole number below `n`, each as likely
    /// as the others. As many bits as `n` has are drawn again until they spell a number
    /// below `n` (`_randbelow_with_getrandbits`).
    pub fn randrange(&mut self, n: u64) -> u64 {
        assert!(n > 0, "randrange of an empty range");
        let bits = u64::BITS - n.leading_zeros();
        loop {
            let drawn = self.getrandbits(bits);
            if drawn < n {
                return drawn;
            }
        }
    }

    /// `shuffle(items)`: from the last place down to the second, each item changes places
    /// with one drawn among the places up to its own.
    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.randrange(i as u64 + 1) as usize;
            items.swap(i, j);
        }
    }
}

#[cfg(test)]
#[path = "random_tests.rs"]
mod tests;
