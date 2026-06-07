// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # POD process-noise models
//!
//! ## Scientific scope
//! Q-matrix construction for EKF and related sequential estimators.

pub mod noise;

pub use noise::{
    GaussMarkovParams, PiecewiseSegment, ProcessNoise, ProcessNoiseModel, WhiteAccelPsd,
};
