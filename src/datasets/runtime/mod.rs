// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Runtime acquisition of `RuntimeDownload` datasets.
//!
//! This module is only compiled when the `runtime-data` feature is enabled.
//! It provides [`DatasetManager`], which downloads and caches
//! [`super::Acquisition::RuntimeDownload`] datasets. Calling acquire methods
//! on `Embedded` or `ExternalProvider` datasets returns
//! [`super::DatasetError::NotDownloadable`].

mod cache;
mod download;
mod manager;

pub use download::ProgressCallback;
pub use manager::DatasetManager;
