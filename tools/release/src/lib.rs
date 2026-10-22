//! Release metadata helpers used by tag-triggered publication.

mod changelog;
mod validate;

pub use changelog::{
    extract_version_section, format_release_notes, ChangelogError, ChangelogSection,
};
pub use validate::{
    crates_io_has_version, default_root, package_meta, parse_tag, run_validation, PackageMeta,
    ValidationError,
};
