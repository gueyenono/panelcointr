//! Random numbers for the null-distribution simulations: xoshiro256++
//! (Blackman and Vigna 2019) seeded through SplitMix64, and standard normal
//! draws by Marsaglia's polar method. Each block of replicates gets its own
//! stream, derived from a seed and the block index, so the simulations do not
//! depend on the number of threads.

/// SplitMix64 step: advances `state` and returns the next output.
#[inline]
pub fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
    spare: Option<f64>,
}

impl Rng {
    /// Generator with the given internal state (for the reference test vectors).
    pub fn from_state(s: [u64; 4]) -> Self {
        Rng { s, spare: None }
    }

    /// Stream `stream` of `seed`: the state is four SplitMix64 outputs from a
    /// starting value that mixes both.
    pub fn stream(seed: u64, stream: u64) -> Self {
        let mut mix = stream;
        let mut sm = seed ^ splitmix64(&mut mix);
        let s = [splitmix64(&mut sm), splitmix64(&mut sm), splitmix64(&mut sm), splitmix64(&mut sm)];
        Rng { s, spare: None }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Uniform on `[0, 1)` with 53 random bits.
    #[inline]
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Standard normal by Marsaglia's polar method; the second draw of each
    /// pair is kept.
    #[inline]
    pub fn normal(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        loop {
            let u = 2.0 * self.uniform() - 1.0;
            let v = 2.0 * self.uniform() - 1.0;
            let s = u * u + v * v;
            if s > 0.0 && s < 1.0 {
                let m = (-2.0 * s.ln() / s).sqrt();
                self.spare = Some(v * m);
                return u * m;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xoshiro_reference_outputs() {
        // Reference implementation (prng.di.unimi.it) from state {1, 2, 3, 4}.
        let mut g = Rng::from_state([1, 2, 3, 4]);
        assert_eq!(g.next_u64(), 41_943_041);
        assert_eq!(g.next_u64(), 58_720_359);
        assert_eq!(g.next_u64(), 3_588_806_011_781_223);
        assert_eq!(g.next_u64(), 3_591_011_842_654_386);
    }

    #[test]
    fn splitmix_reference_output() {
        // First output of SplitMix64 seeded with 0 (Vigna's splitmix64.c).
        let mut s = 0u64;
        assert_eq!(splitmix64(&mut s), 0xE220_A839_7B1D_CDAF);
    }

    #[test]
    fn normal_moments() {
        let mut g = Rng::stream(1, 0);
        let n = 200_000;
        let (mut s1, mut s2) = (0.0, 0.0);
        for _ in 0..n {
            let z = g.normal();
            s1 += z;
            s2 += z * z;
        }
        let m = s1 / n as f64;
        let v = s2 / n as f64 - m * m;
        assert!(m.abs() < 0.01, "{m}");
        assert!((v - 1.0).abs() < 0.01, "{v}");
    }

    #[test]
    fn streams_differ() {
        let a = Rng::stream(7, 0).next_u64();
        let b = Rng::stream(7, 1).next_u64();
        let c = Rng::stream(8, 0).next_u64();
        assert!(a != b && a != c && b != c);
    }
}
