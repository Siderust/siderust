mod brightness;
mod builder;
pub mod csv;
mod error;
mod magnitude;
mod map;
mod photometry;
mod provenance;
mod record;
mod validation;

pub use brightness::StellarSurfaceBrightness;
pub use builder::StellarSurfaceBrightnessMapBuilder;
pub use error::{Result, StellarMapError};
pub use magnitude::ApparentMagnitude;
pub use map::StellarSurfaceBrightnessMap;
pub use photometry::flux_10mag_units;
pub use provenance::StellarMapProvenance;
pub use record::StellarCatalogueRecord;
pub use validation::{
    validate_flux_conservation, validate_no_longitude_wrap_artifact,
    validate_plane_pole_contrast, validate_stellar_map_values,
};
