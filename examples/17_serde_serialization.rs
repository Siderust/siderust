// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Serde Serialization Example (prefixed)
//!
//! Run with: `cargo run --example 17_serde_serialization --features serde`

use qtty::*;
use serde::{Deserialize, Serialize};
use serde_json;
use siderust::coordinates::{cartesian, frames, spherical};
use siderust::time::JulianDate;
use std::fs;

fn main() {
    println!("Serde serialization (prefixed)");
}
