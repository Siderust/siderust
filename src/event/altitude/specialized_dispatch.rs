// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Crate-private fast paths for in-crate targets with specialized search.

use super::provider::AltitudeProvider;
use super::search::{InternalSearchConfig, SearchOpts};
use super::types::CrossingEvent;
use crate::bodies::{solar_system, Star};
use crate::coordinates::centers::Geodetic;
use crate::coordinates::frames::ECEF;
use crate::coordinates::spherical::direction;
use crate::qtty::*;
use crate::time::{Interval, ModifiedJulianDate};
use std::any::Any;

fn internal_config(opts: SearchOpts) -> InternalSearchConfig {
    InternalSearchConfig::from_public_opts(opts)
}

/// Attempt a specialized above-threshold search for known in-crate target types.
pub(crate) fn above_threshold_search<T: AltitudeProvider + Any + 'static>(
    target: &T,
    observer: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Option<Vec<Interval<ModifiedJulianDate>>> {
    let any = target as &dyn Any;
    if let Some(_sun) = any.downcast_ref::<solar_system::Sun>() {
        return Some(crate::event::solar::solar_above_threshold_impl(
            observer,
            window,
            threshold,
            internal_config(opts),
        ));
    }
    if let Some(_moon) = any.downcast_ref::<solar_system::Moon>() {
        return Some(crate::event::lunar::lunar_above_threshold_impl(
            observer,
            window,
            threshold,
            internal_config(opts),
        ));
    }
    if let Some(dir) = any.downcast_ref::<direction::ICRS>() {
        return Some(crate::event::stellar::find_star_above_periods(
            dir.ra(),
            dir.dec(),
            observer,
            window,
            threshold,
        ));
    }
    if let Some(star) = any.downcast_ref::<Star<'static>>() {
        let dir = direction::ICRS::from(star);
        return Some(crate::event::stellar::find_star_above_periods(
            dir.ra(),
            dir.dec(),
            observer,
            window,
            threshold,
        ));
    }
    None
}

/// Attempt a specialized below-threshold search for known in-crate target types.
pub(crate) fn below_threshold_search<T: AltitudeProvider + Any + 'static>(
    target: &T,
    observer: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Option<Vec<Interval<ModifiedJulianDate>>> {
    let any = target as &dyn Any;
    if let Some(_sun) = any.downcast_ref::<solar_system::Sun>() {
        return Some(crate::event::solar::solar_below_threshold_impl(
            observer,
            window,
            threshold,
            internal_config(opts),
        ));
    }
    if let Some(_moon) = any.downcast_ref::<solar_system::Moon>() {
        return Some(crate::event::lunar::lunar_below_threshold_impl(
            observer,
            window,
            threshold,
            internal_config(opts),
        ));
    }
    if let Some(dir) = any.downcast_ref::<direction::ICRS>() {
        return Some(crate::event::stellar::find_star_below_periods(
            dir.ra(),
            dir.dec(),
            observer,
            window,
            threshold,
        ));
    }
    if let Some(star) = any.downcast_ref::<Star<'static>>() {
        let dir = direction::ICRS::from(star);
        return Some(crate::event::stellar::find_star_below_periods(
            dir.ra(),
            dir.dec(),
            observer,
            window,
            threshold,
        ));
    }
    None
}

/// Attempt a specialized altitude-range search for known in-crate target types.
pub(crate) fn altitude_range_search<T: AltitudeProvider + Any + 'static>(
    target: &T,
    observer: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    min_altitude: Degrees,
    max_altitude: Degrees,
    opts: SearchOpts,
) -> Option<Vec<Interval<ModifiedJulianDate>>> {
    let any = target as &dyn Any;
    if let Some(_sun) = any.downcast_ref::<solar_system::Sun>() {
        return Some(crate::event::solar::solar_altitude_ranges_impl(
            observer,
            window,
            (min_altitude, max_altitude),
            internal_config(opts),
        ));
    }
    if let Some(_moon) = any.downcast_ref::<solar_system::Moon>() {
        return Some(crate::event::lunar::lunar_altitude_ranges_impl(
            observer,
            window,
            (min_altitude, max_altitude),
            internal_config(opts),
        ));
    }
    if let Some(dir) = any.downcast_ref::<direction::ICRS>() {
        return Some(crate::event::stellar::find_star_range_periods(
            dir.ra(),
            dir.dec(),
            observer,
            window,
            (min_altitude, max_altitude),
        ));
    }
    if let Some(star) = any.downcast_ref::<Star<'static>>() {
        let dir = direction::ICRS::from(star);
        return Some(crate::event::stellar::find_star_range_periods(
            dir.ra(),
            dir.dec(),
            observer,
            window,
            (min_altitude, max_altitude),
        ));
    }
    None
}

/// Attempt a specialized crossing search for known in-crate target types.
pub(crate) fn crossings_search<T: AltitudeProvider + Any + 'static>(
    target: &T,
    observer: Geodetic<ECEF>,
    window: Interval<ModifiedJulianDate>,
    threshold: Degrees,
    opts: SearchOpts,
) -> Option<Vec<CrossingEvent>> {
    let any = target as &dyn Any;
    if let Some(_sun) = any.downcast_ref::<solar_system::Sun>() {
        return Some(crate::event::solar::solar_crossings_impl(
            observer,
            window,
            threshold,
            internal_config(opts),
        ));
    }
    if let Some(_moon) = any.downcast_ref::<solar_system::Moon>() {
        return Some(crate::event::lunar::lunar_crossings_impl(
            observer,
            window,
            threshold,
            internal_config(opts),
        ));
    }
    None
}
