//! Check relative inline links and heading fragments in repository documentation.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::LazyLock;

use percent_encoding::percent_decode_str;
use regex::Regex;

use siderust_release::default_root;

static FENCE: LazyLock<Regex> = LazyLock::new(|| regex(r"^ {0,3}(`{3,}|~{3,})(.*)$"));
static LINK: LazyLock<Regex> =
    LazyLock::new(|| regex(r#"!?(?:\[[^\]\n]*\])\(<?([^\s)>]+)>?(?: +"[^"]*")?\)"#));
static HEADING: LazyLock<Regex> = LazyLock::new(|| regex(r"^ {0,3}#{1,6}\s+(.+?)\s*#*\s*$"));
static HEADING_LINK: LazyLock<Regex> = LazyLock::new(|| regex(r"\[([^\]]+)\]\([^)]+\)"));
static ANCHOR: LazyLock<Regex> = LazyLock::new(|| regex(r#"<a\s+(?:id|name)="([^"]+)""#));
static SCHEME: LazyLock<Regex> = LazyLock::new(|| regex(r"^[A-Za-z][A-Za-z0-9+.-]*:"));
static BACKTICKS: LazyLock<Regex> = LazyLock::new(|| regex(r"`+"));

fn regex(pattern: &str) -> Regex {
    Regex::new(pattern).expect("valid built-in pattern")
}

fn markdown_lines(text: &str) -> Result<Vec<(usize, &str)>, usize> {
    let mut lines = Vec::new();
    let mut fence: Option<(&str, usize)> = None;
    for (index, line) in text.lines().enumerate() {
        let number = index + 1;
        let marker = FENCE.captures(line);
        match (fence, marker) {
            (Some((open, _)), Some(close))
                if close[1].starts_with(open) && close[2].trim().is_empty() =>
            {
                fence = None;
            }
            (Some(_), _) => {}
            (None, Some(marker)) => fence = Some((marker.get(1).unwrap().as_str(), number)),
            (None, None) => lines.push((number, line)),
        }
    }
    match fence {
        Some((_, number)) => Err(number),
        None => Ok(lines),
    }
}

fn anchors(lines: &[(usize, &str)]) -> HashSet<String> {
    let mut found = HashSet::new();
    let mut counts = HashMap::new();
    for (_, line) in lines {
        found.extend(
            ANCHOR
                .captures_iter(line)
                .map(|anchor| anchor[1].to_owned()),
        );
        if let Some(heading) = HEADING.captures(line) {
            let label = HEADING_LINK.replace_all(&heading[1], "$1");
            let slug: String = label
                .to_lowercase()
                .chars()
                .filter(|ch| ch.is_alphanumeric() || matches!(ch, '_' | '-' | ' '))
                .map(|ch| if ch == ' ' { '-' } else { ch })
                .collect();
            let mut anchor = slug.clone();
            let count = counts.entry(slug.clone()).or_insert(0);
            while found.contains(&anchor) {
                *count += 1;
                anchor = format!("{slug}-{count}");
            }
            found.insert(anchor);
        }
    }
    found
}

fn inline_code(line: &str) -> Vec<std::ops::Range<usize>> {
    let runs: Vec<_> = BACKTICKS.find_iter(line).collect();
    let mut spans = Vec::new();
    let mut index = 0;
    while let Some(open) = runs.get(index) {
        if let Some(close) = runs[index + 1..]
            .iter()
            .position(|run| run.len() == open.len())
        {
            index += close + 1;
            spans.push(open.start()..runs[index].end());
        }
        index += 1;
    }
    spans
}

fn local_destination(destination: &str) -> Option<(&str, &str)> {
    if destination.starts_with("//") || SCHEME.is_match(destination) {
        return None;
    }
    let (path, fragment) = destination.split_once('#').unwrap_or((destination, ""));
    let path = path.split_once('?').map_or(path, |(path, _)| path);
    Some((path, fragment))
}

fn is_markdown(path: &Path) -> bool {
    path.extension().is_some_and(|ext| ext == "md")
}

fn collect_markdown(dir: &Path, pages: &mut Vec<PathBuf>) -> io::Result<()> {
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        if entry.file_type()?.is_dir() {
            collect_markdown(&entry.path(), pages)?;
        } else if is_markdown(&entry.path()) {
            pages.push(entry.path());
        }
    }
    Ok(())
}

fn pages(root: &Path) -> io::Result<Vec<PathBuf>> {
    let mut pages = Vec::new();
    for dir in ["doc", "benches"] {
        collect_markdown(&root.join(dir), &mut pages)?;
    }
    for entry in fs::read_dir(root.join("scripts"))? {
        let path = entry?.path().join("README.md");
        if path.is_file() {
            pages.push(path);
        }
    }
    pages.extend(
        [
            "README.md",
            "CHANGELOG.md",
            "examples/README.md",
            "scripts/README.md",
            "siderust-ffi/README.md",
            "tests/fixtures/lunar_de440/README.md",
            "tests/fixtures/sgp4/README.md",
        ]
        .map(|path| root.join(path)),
    );
    pages.sort();
    Ok(pages)
}

fn document(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|error| format!("{}: {error}", path.display()))
}

fn checked_lines<'a>(path: &Path, text: &'a str) -> Result<Vec<(usize, &'a str)>, String> {
    markdown_lines(text)
        .map_err(|number| format!("{}:{number}: unclosed code fence", path.display()))
}

fn check_links(root: &Path, pages: &[PathBuf]) -> Result<Vec<String>, String> {
    let mut failures = Vec::new();
    let mut headings = HashMap::new();
    for page in pages {
        let relative = page.strip_prefix(root).unwrap_or(page);
        let text = document(page)?;
        for (number, line) in checked_lines(relative, &text)? {
            let code = inline_code(line);
            for link in LINK.captures_iter(line) {
                if code
                    .iter()
                    .any(|span| span.contains(&link.get(0).unwrap().start()))
                {
                    continue;
                }
                let destination = &link[1];
                let Some((path, fragment)) = local_destination(destination) else {
                    continue;
                };
                let target = if path.is_empty() {
                    page.clone()
                } else {
                    let path = percent_decode_str(path).decode_utf8_lossy();
                    page.parent().unwrap_or(root).join(path.as_ref())
                };
                let mut target = match target.canonicalize() {
                    Ok(path) => path,
                    Err(error) if error.kind() == io::ErrorKind::NotFound => {
                        failures.push(format!(
                            "{}:{number}: missing target {destination}",
                            relative.display()
                        ));
                        continue;
                    }
                    Err(error) => {
                        return Err(format!(
                            "{}:{number}: {destination}: {error}",
                            relative.display()
                        ))
                    }
                };
                if target.is_dir() && target.join("README.md").is_file() {
                    target.push("README.md");
                }
                if !fragment.is_empty() && is_markdown(&target) {
                    if !headings.contains_key(&target) {
                        let text = document(&target)?;
                        headings.insert(target.clone(), anchors(&checked_lines(&target, &text)?));
                    }
                    let fragment = percent_decode_str(fragment).decode_utf8_lossy();
                    if !headings[&target].contains(fragment.as_ref()) {
                        failures.push(format!(
                            "{}:{number}: missing heading {destination}",
                            relative.display()
                        ));
                    }
                }
            }
        }
    }
    Ok(failures)
}

fn main() -> ExitCode {
    let root = default_root();
    let result = pages(&root)
        .map_err(|error| format!("could not list Markdown files: {error}"))
        .and_then(|pages| check_links(&root, &pages).map(|failures| (pages.len(), failures)));
    match result {
        Ok((count, failures)) if failures.is_empty() => {
            println!("Checked relative links and code fences in {count} Markdown files.");
            ExitCode::SUCCESS
        }
        Ok((_, failures)) => {
            eprintln!("{}", failures.join("\n"));
            ExitCode::FAILURE
        }
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checks_local_destinations_and_ignores_code_and_external_links() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path();
        let page = root.join("README.md");
        fs::create_dir(root.join("with space")).unwrap();
        fs::write(
            root.join("with space/README.md"),
            "# Café\n# Café\n<a name=\"old\"></a>",
        )
        .unwrap();
        fs::write(&page, concat!(
            "# Home\n",
            "[`valid`](with%20space?query#caf%C3%A9-1) [legacy](with%20space/#old)\n",
            "[self](#home) [external](https://example.com/missing) [remote](//example.com/missing)\n",
            "``[example](missing) ` code``\n",
            "~~~~rust\n[example](missing)\n~~~\n[example](missing)\n~~~~\n",
            "[bad](missing) [bad heading](#missing)\n",
        )).unwrap();
        assert_eq!(
            check_links(root, &[page]).unwrap(),
            [
                "README.md:10: missing target missing",
                "README.md:10: missing heading #missing",
            ]
        );
    }

    #[test]
    fn reports_unclosed_fence_at_its_opening() {
        assert_eq!(markdown_lines("# Title\n   ````rust\n```\n"), Err(2));
    }
}
