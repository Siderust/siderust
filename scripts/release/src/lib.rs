//! Release metadata helpers used by tag-triggered publication.

mod changelog;
mod publish;
mod validate;

pub use changelog::{
    extract_version_section, format_release_notes, ChangelogError, ChangelogSection,
};
pub use publish::{crates_io_has_version, publish};
pub use validate::{
    default_root, explicit_version, package_meta, parse_tag, run_validation, PackageMeta,
    ValidationError,
};
