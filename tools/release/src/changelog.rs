//! Keep a Changelog section extraction for GitHub Releases.

use std::fmt;

/// One `## [version]` block from CHANGELOG.md.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChangelogSection {
    /// Label inside the brackets (`0.15.0`, `Unreleased`, …).
    pub label: String,
    /// Optional date after ` - `.
    pub date: Option<String>,
    /// Section body excluding the header line.
    pub body: String,
}

/// Invalid or missing changelog content.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChangelogError(pub String);

impl fmt::Display for ChangelogError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ChangelogError {}

fn is_section_header(line: &str) -> Option<(&str, Option<&str>)> {
    let trimmed = line.trim_end();
    let rest = trimmed.strip_prefix("## [")?;
    let close = rest.find(']')?;
    let label = rest[..close].trim();
    let after = rest[close + 1..].trim();
    let date = if let Some(date) = after.strip_prefix('-') {
        let date = date.trim();
        if date.is_empty() {
            None
        } else {
            Some(date)
        }
    } else if after.is_empty() {
        None
    } else {
        return None;
    };
    Some((label, date))
}

/// Parse all top-level version sections in changelog order.
pub(crate) fn iter_sections(text: &str) -> Vec<ChangelogSection> {
    let lines: Vec<&str> = text.lines().collect();
    let mut headers: Vec<(usize, String, Option<String>)> = Vec::new();
    for (index, line) in lines.iter().enumerate() {
        if let Some((label, date)) = is_section_header(line) {
            headers.push((index, label.to_owned(), date.map(str::to_owned)));
        }
    }

    let mut sections = Vec::with_capacity(headers.len());
    for (position, (start_line, label, date)) in headers.iter().enumerate() {
        let body_start = start_line + 1;
        let body_end = headers
            .get(position + 1)
            .map(|(next, _, _)| *next)
            .unwrap_or(lines.len());
        let body = lines[body_start..body_end].join("\n");
        sections.push(ChangelogSection {
            label: label.clone(),
            date: date.clone(),
            body: body.trim_matches('\n').to_owned(),
        });
    }
    sections
}

/// Return `X.Y.Z` when *label* is a version heading, else `None`.
pub(crate) fn normalize_version(label: &str) -> Option<String> {
    if label.eq_ignore_ascii_case("unreleased") {
        return None;
    }
    let candidate = label
        .strip_prefix('v')
        .or_else(|| label.strip_prefix('V'))
        .unwrap_or(label);
    if is_semver_like(candidate) {
        Some(candidate.to_owned())
    } else {
        None
    }
}

fn is_semver_like(value: &str) -> bool {
    let mut parts = value.splitn(2, ['-', '+']);
    let core = parts.next().unwrap_or("");
    let mut nums = core.split('.');
    let (Some(major), Some(minor), Some(patch)) = (nums.next(), nums.next(), nums.next()) else {
        return false;
    };
    if nums.next().is_some() {
        return false;
    }
    major.chars().all(|c| c.is_ascii_digit())
        && !major.is_empty()
        && minor.chars().all(|c| c.is_ascii_digit())
        && !minor.is_empty()
        && patch.chars().all(|c| c.is_ascii_digit())
        && !patch.is_empty()
}

/// Return the unique nonempty section for *version* (with or without a leading `v`).
pub fn extract_version_section(
    text: &str,
    version: &str,
) -> Result<ChangelogSection, ChangelogError> {
    let wanted = version
        .strip_prefix('v')
        .or_else(|| version.strip_prefix('V'))
        .unwrap_or(version);
    let matches: Vec<_> = iter_sections(text)
        .into_iter()
        .filter(|section| normalize_version(&section.label).as_deref() == Some(wanted))
        .collect();
    match matches.len() {
        0 => Err(ChangelogError(format!(
            "CHANGELOG.md has no section for version {wanted}"
        ))),
        1 => {
            let section = matches.into_iter().next().expect("len checked");
            if section.body.trim().is_empty() {
                Err(ChangelogError(format!(
                    "CHANGELOG.md section for {wanted} is empty"
                )))
            } else {
                Ok(section)
            }
        }
        n => Err(ChangelogError(format!(
            "CHANGELOG.md has {n} sections for version {wanted}"
        ))),
    }
}

/// Build GitHub Release notes from a changelog section body.
pub fn format_release_notes(section: &ChangelogSection, version: &str) -> String {
    let heading = match &section.date {
        Some(date) => format!("## [{version}] - {date}"),
        None => format!("## [{version}]"),
    };
    format!("{heading}\n\n{}\n", section.body.trim())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_matching_section() {
        let text = "\
# Changelog

## [Unreleased]

### Added
- pending

## [0.15.0]

### Added
- feature one

## [0.14.1] - 2026-10-14

### Fixed
- bug
";
        let section = extract_version_section(text, "0.15.0").unwrap();
        assert!(section.body.contains("feature one"));
        assert!(!section.body.contains("bug"));
        assert!(section.date.is_none());
        let notes = format_release_notes(&section, "0.15.0");
        assert!(notes.starts_with("## [0.15.0]\n"));
        assert!(notes.contains("feature one"));
    }

    #[test]
    fn extracts_dated_section() {
        let text = "## [0.14.1] - 2026-10-14\n\n### Fixed\n- bug\n";
        let section = extract_version_section(text, "v0.14.1").unwrap();
        assert_eq!(section.date.as_deref(), Some("2026-10-14"));
        let notes = format_release_notes(&section, "0.14.1");
        assert!(notes.starts_with("## [0.14.1] - 2026-10-14\n"));
    }

    #[test]
    fn missing_section() {
        let err = extract_version_section("## [Unreleased]\n\n- x\n", "0.15.0").unwrap_err();
        assert!(err.0.contains("no section"));
    }

    #[test]
    fn empty_section() {
        let text = "## [0.15.0]\n\n## [0.14.1]\n\n- older\n";
        let err = extract_version_section(text, "0.15.0").unwrap_err();
        assert!(err.0.contains("empty"));
    }

    #[test]
    fn duplicate_sections() {
        let text = "## [0.15.0]\n\n- a\n\n## [0.15.0]\n\n- b\n";
        let err = extract_version_section(text, "0.15.0").unwrap_err();
        assert!(err.0.contains("2 sections"));
    }
}
