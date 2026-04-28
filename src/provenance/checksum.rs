// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Build-time integrity verification for embedded data tables.
//!
//! Datasets that ship in the crate via [`include_str!`] / [`include_bytes!`]
//! should be pinned to a SHA-256 hash so that an accidental edit to the
//! bundled file (or a corrupted vendored copy) is rejected at compile time
//! rather than producing silently wrong scientific results.
//!
//! ## Usage
//!
//! ```ignore
//! use siderust::assert_data_checksum;
//!
//! const RAW: &str = include_str!("../../data/o3trans.dat");
//! assert_data_checksum!(
//!     "o3trans.dat",
//!     RAW.as_bytes(),
//!     "cb06c173f393d6d55e3c39551665abb8f5d6c1a846cd0fd739a15d0155f94502"
//! );
//! ```
//!
//! The macro expands to a `const _: () = …` block that calls
//! [`assert_sha256_eq`] in const context. SHA-256 is implemented as a
//! `const fn` here ([`sha256`]), so the check is fully evaluated by the
//! compiler — no new runtime dependency, no extra crates.
//!
//! ## Updating a pinned hash
//!
//! When a bundled file legitimately changes (e.g. a new release of an
//! upstream catalog) recompute the hash with:
//!
//! ```text
//! sha256sum siderust/data/<file>
//! ```
//!
//! or run the helper test:
//!
//! ```text
//! cargo test -p siderust provenance::checksum::dev::print_hashes -- --nocapture --ignored
//! ```
//!
//! Then update the pinned literal at the call site **and** record the
//! reason for the update in the commit message / file-level provenance
//! note. Hashes must never be updated without verifying that the new
//! bytes are intentional.
//!
//! ## Defense in depth
//!
//! Each protected file is also covered by a `#[test]` that recomputes
//! the SHA-256 at runtime and asserts it matches the pinned value, so
//! any future refactor that bypasses the const path still fails
//! `cargo test`.

// ── const-evaluable SHA-256 ──────────────────────────────────────────────────

const K: [u32; 64] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
];

const H0: [u32; 8] = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
];

const fn rotr(x: u32, n: u32) -> u32 {
    (x >> n) | (x << (32 - n))
}

const fn compress(state: &mut [u32; 8], block: &[u8; 64]) {
    let mut w = [0u32; 64];
    let mut i = 0;
    while i < 16 {
        w[i] = ((block[i * 4] as u32) << 24)
            | ((block[i * 4 + 1] as u32) << 16)
            | ((block[i * 4 + 2] as u32) << 8)
            | (block[i * 4 + 3] as u32);
        i += 1;
    }
    let mut i = 16;
    while i < 64 {
        let s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        let s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16]
            .wrapping_add(s0)
            .wrapping_add(w[i - 7])
            .wrapping_add(s1);
        i += 1;
    }

    let mut a = state[0];
    let mut b = state[1];
    let mut c = state[2];
    let mut d = state[3];
    let mut e = state[4];
    let mut f = state[5];
    let mut g = state[6];
    let mut h = state[7];

    let mut i = 0;
    while i < 64 {
        let s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        let ch = (e & f) ^ ((!e) & g);
        let t1 = h
            .wrapping_add(s1)
            .wrapping_add(ch)
            .wrapping_add(K[i])
            .wrapping_add(w[i]);
        let s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        let mj = (a & b) ^ (a & c) ^ (b & c);
        let t2 = s0.wrapping_add(mj);
        h = g;
        g = f;
        f = e;
        e = d.wrapping_add(t1);
        d = c;
        c = b;
        b = a;
        a = t1.wrapping_add(t2);
        i += 1;
    }

    state[0] = state[0].wrapping_add(a);
    state[1] = state[1].wrapping_add(b);
    state[2] = state[2].wrapping_add(c);
    state[3] = state[3].wrapping_add(d);
    state[4] = state[4].wrapping_add(e);
    state[5] = state[5].wrapping_add(f);
    state[6] = state[6].wrapping_add(g);
    state[7] = state[7].wrapping_add(h);
}

/// Compute the SHA-256 digest of `input` at compile time or runtime.
///
/// Returns the 32-byte raw digest (big-endian word order, MSB-first within
/// each word — i.e. the canonical SHA-256 byte sequence).
pub const fn sha256(input: &[u8]) -> [u8; 32] {
    let mut state = H0;
    let len = input.len();
    let bit_len = (len as u64).wrapping_mul(8);

    // Process every full 64-byte block straight from the input.
    let mut off = 0;
    while off + 64 <= len {
        let mut block = [0u8; 64];
        let mut k = 0;
        while k < 64 {
            block[k] = input[off + k];
            k += 1;
        }
        compress(&mut state, &block);
        off += 64;
    }

    // The final 1 or 2 padded blocks are assembled in a 128-byte tail
    // buffer (worst case = 55-byte remainder needs one extra block, so
    // 64 < remainder + 1 + len(8) ≤ 128).
    let rem = len - off;
    let mut tail = [0u8; 128];
    let mut k = 0;
    while k < rem {
        tail[k] = input[off + k];
        k += 1;
    }
    tail[rem] = 0x80;

    let pad_zeros = if rem % 64 < 56 {
        55 - (rem % 64)
    } else {
        55 + 64 - (rem % 64)
    };
    let final_len = rem + 1 + pad_zeros + 8;

    let len_off = final_len - 8;
    tail[len_off]     = (bit_len >> 56) as u8;
    tail[len_off + 1] = (bit_len >> 48) as u8;
    tail[len_off + 2] = (bit_len >> 40) as u8;
    tail[len_off + 3] = (bit_len >> 32) as u8;
    tail[len_off + 4] = (bit_len >> 24) as u8;
    tail[len_off + 5] = (bit_len >> 16) as u8;
    tail[len_off + 6] = (bit_len >> 8)  as u8;
    tail[len_off + 7] =  bit_len        as u8;

    let mut bi = 0;
    while bi < final_len {
        let mut block = [0u8; 64];
        let mut m = 0;
        while m < 64 {
            block[m] = tail[bi + m];
            m += 1;
        }
        compress(&mut state, &block);
        bi += 64;
    }

    let mut out = [0u8; 32];
    let mut o = 0;
    while o < 8 {
        let v = state[o];
        out[o * 4]     = (v >> 24) as u8;
        out[o * 4 + 1] = (v >> 16) as u8;
        out[o * 4 + 2] = (v >> 8)  as u8;
        out[o * 4 + 3] =  v        as u8;
        o += 1;
    }
    out
}

// ── hex helpers ──────────────────────────────────────────────────────────────

const fn hex_nibble(b: u8) -> u8 {
    match b {
        b'0'..=b'9' => b - b'0',
        b'a'..=b'f' => b - b'a' + 10,
        b'A'..=b'F' => b - b'A' + 10,
        _ => panic!("assert_data_checksum!: expected SHA-256 hex digit"),
    }
}

/// Decode a 64-character lowercase/uppercase hex string into a 32-byte SHA-256.
pub const fn hex32(s: &str) -> [u8; 32] {
    let bytes = s.as_bytes();
    if bytes.len() != 64 {
        panic!("assert_data_checksum!: expected a 64-character SHA-256 hex string");
    }
    let mut out = [0u8; 32];
    let mut i = 0;
    while i < 32 {
        out[i] = (hex_nibble(bytes[2 * i]) << 4) | hex_nibble(bytes[2 * i + 1]);
        i += 1;
    }
    out
}

/// Render a 32-byte digest as a lowercase 64-char hex string.
pub fn to_hex(digest: &[u8; 32]) -> String {
    const TABLE: &[u8; 16] = b"0123456789abcdef";
    let mut s = String::with_capacity(64);
    for &b in digest {
        s.push(TABLE[(b >> 4) as usize] as char);
        s.push(TABLE[(b & 0x0f) as usize] as char);
    }
    s
}

// ── compile-time assertion machinery ─────────────────────────────────────────

/// Compare a SHA-256 of `data` against the hex-encoded `expected` value in
/// const context. Triggers a `panic!` (which becomes a `compile_error!` when
/// invoked from a `const _: () = …` binding) on mismatch.
///
/// The `name` argument is currently informational; it is retained in the
/// signature so that future revisions can surface it via a richer
/// diagnostic without breaking the macro contract.
pub const fn assert_sha256_eq(_name: &'static str, data: &[u8], expected: &str) {
    let actual = sha256(data);
    let expected = hex32(expected);
    let mut i = 0;
    while i < 32 {
        if actual[i] != expected[i] {
            panic!(
                "siderust: SHA-256 mismatch for an embedded data table — \
                 the bundled file no longer matches the pinned hash. If \
                 the data change is intentional, recompute the hash with \
                 `sha256sum` (or the `print_hashes` test) and update the \
                 pinned literal."
            );
        }
        i += 1;
    }
}

/// Pin the SHA-256 of an embedded data blob at compile time.
///
/// Expands to a `const _: () = …` binding that const-evaluates
/// [`assert_sha256_eq`] against the byte slice `$data`. A mismatch becomes
/// a hard compile error.
///
/// `$data` must be `&[u8]`. For `&'static str` data (`include_str!`) call
/// `.as_bytes()` at the macro site:
///
/// ```ignore
/// const RAW: &str = include_str!("../../data/o3trans.dat");
/// siderust::assert_data_checksum!(
///     "siderust/data/o3trans.dat",
///     RAW.as_bytes(),
///     "cb06c173f393d6d55e3c39551665abb8f5d6c1a846cd0fd739a15d0155f94502"
/// );
/// ```
#[macro_export]
macro_rules! assert_data_checksum {
    ($name:literal, $data:expr, $expected_hex:literal) => {
        #[allow(long_running_const_eval)]
        const _: () = $crate::provenance::checksum::assert_sha256_eq(
            $name,
            $data,
            $expected_hex,
        );
    };
}

// ── developer helpers ────────────────────────────────────────────────────────

#[cfg(test)]
pub(crate) mod dev {
    //! Helpers for refreshing pinned hashes during dataset updates.

    use super::{sha256, to_hex};

    /// Print the SHA-256 of every dataset currently protected by
    /// [`assert_data_checksum!`]. Run with:
    ///
    /// ```text
    /// cargo test -p siderust provenance::checksum::dev::print_hashes \
    ///     -- --nocapture --ignored
    /// ```
    #[test]
    #[ignore = "developer helper — prints pinned-dataset hashes"]
    fn print_hashes() {
        let datasets: &[(&str, &[u8])] = &[
            (
                "siderust/data/o3trans.dat",
                include_bytes!("../../data/o3trans.dat"),
            ),
            (
                "siderust/data/passbands/bessell1990/U.dat",
                include_bytes!("../../data/passbands/bessell1990/U.dat"),
            ),
            (
                "siderust/data/passbands/bessell1990/B.dat",
                include_bytes!("../../data/passbands/bessell1990/B.dat"),
            ),
            (
                "siderust/data/passbands/bessell1990/V.dat",
                include_bytes!("../../data/passbands/bessell1990/V.dat"),
            ),
            (
                "siderust/data/passbands/bessell1990/R.dat",
                include_bytes!("../../data/passbands/bessell1990/R.dat"),
            ),
            (
                "siderust/data/passbands/bessell1990/I.dat",
                include_bytes!("../../data/passbands/bessell1990/I.dat"),
            ),
        ];
        for (name, data) in datasets {
            println!("{}  {}", to_hex(&sha256(data)), name);
        }
    }
}

// ── self-tests ───────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Known-answer tests for the const-fn SHA-256 implementation.
    /// Vectors taken from FIPS 180-2 / NIST CAVS.
    #[test]
    fn sha256_known_answers() {
        // empty string
        assert_eq!(
            to_hex(&sha256(b"")),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        // "abc"
        assert_eq!(
            to_hex(&sha256(b"abc")),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        // 56-byte boundary case (exercises the two-block padding path).
        assert_eq!(
            to_hex(&sha256(
                b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"
            )),
            "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
        );
        // 1 000 000 'a' characters — exercises multi-block streaming.
        let mut million_a = Vec::with_capacity(1_000_000);
        million_a.resize(1_000_000, b'a');
        assert_eq!(
            to_hex(&sha256(&million_a)),
            "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
        );
    }

    #[test]
    fn const_eval_path_matches_runtime_path() {
        const DIGEST: [u8; 32] = sha256(b"abc");
        assert_eq!(
            DIGEST,
            hex32("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad")
        );
    }

    #[test]
    fn hex32_round_trip() {
        let want = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";
        let bytes = hex32(want);
        assert_eq!(to_hex(&bytes), want);
    }

    /// Sanity: `assert_data_checksum!` accepts a matching pinned hash.
    /// (A mismatched hash would fail to compile, which we cannot test
    /// without a separate compile-fail harness.)
    #[test]
    fn macro_accepts_matching_hash() {
        const PAYLOAD: &[u8] = b"abc";
        crate::assert_data_checksum!(
            "test/abc",
            PAYLOAD,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }
}
