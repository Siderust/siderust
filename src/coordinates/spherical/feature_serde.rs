//! Serde implementations for spherical coordinate types.

use super::direction_core::Direction;
use super::position_core::Position;
use crate::coordinates::{centers, frames};
use qtty::LengthUnit;
use serde::{Deserialize, Serialize};

// =============================================================================
// Direction Serde
// =============================================================================

impl<F: frames::SphericalNaming> Serialize for Direction<F> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

impl<'de, F: frames::SphericalNaming> Deserialize<'de> for Direction<F> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        affn::spherical::Direction::<F>::deserialize(deserializer).map(Self::from)
    }
}

// =============================================================================
// Position Serde
// =============================================================================

impl<C, F, U> Serialize for Position<C, F, U>
where
    C: centers::ReferenceCenter,
    C::Params: Serialize,
    F: frames::SphericalNaming,
    U: LengthUnit,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.inner.serialize(serializer)
    }
}

impl<'de, C, F, U> Deserialize<'de> for Position<C, F, U>
where
    C: centers::ReferenceCenter,
    C::Params: Deserialize<'de> + Default,
    F: frames::SphericalNaming,
    U: LengthUnit,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        affn::spherical::Position::<C, F, U>::deserialize(deserializer).map(Self::from)
    }
}
