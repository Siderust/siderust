#!/usr/bin/env python3
"""Check relative inline links and GitHub heading fragments in the manual and benchmark reports."""

import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
LINK = re.compile(r"!?\[[^\]\n]*\]\(<?([^\s)>]+)>?(?: +\"[^\"]*\")?\)")
INLINE_CODE = re.compile(r"(`+).*?\1")
HEADING = re.compile(r"^ {0,3}#{1,6}\s+(.+?)\s*#*\s*$")


def markdown_lines(path):
    fence = None
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        marker = FENCE.match(line)
        if fence:
            if marker:
                closing_run, suffix = marker.groups()
                if closing_run[0] == fence[0] and len(closing_run) >= len(fence) and not suffix.strip():
                    fence = None
        elif marker:
            fence = marker[1]
        else:
            yield number, line
    if fence:
        raise ValueError(f"{path.relative_to(ROOT)}: unclosed code fence")


def anchors(path):
    counts = defaultdict(int)
    found = set()
    for _, line in markdown_lines(path):
        found.update(re.findall(r'<a\s+(?:id|name)="([^"]+)"', line))
        heading = HEADING.match(line)
        if heading:
            label = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", heading[1]).lower()
            slug = re.sub(r"[^\w -]", "", label).replace(" ", "-")
            anchor = slug
            while anchor in found:
                counts[slug] += 1
                anchor = f"{slug}-{counts[slug]}"
            found.add(anchor)
    return found


def main():
    pages = sorted((ROOT / "doc").rglob("*.md"))
    pages.extend(sorted((ROOT / "benches").rglob("*.md")))
    pages.extend(sorted((ROOT / "scripts").glob("README.md")))
    pages.extend(sorted((ROOT / "scripts").glob("*/README.md")))
    pages.extend(ROOT / name for name in (
        "README.md",
        "CHANGELOG.md",
        "examples/README.md",
        "siderust-ffi/README.md",
        "tests/fixtures/lunar_de440/README.md",
        "tests/fixtures/sgp4/README.md",
    ))
    failures = []
    for page in pages:
        for number, line in markdown_lines(page):
            # Keep link labels containing code spans, but ignore links inside code.
            code_spans = [code.span() for code in INLINE_CODE.finditer(line)]
            for match in LINK.finditer(line):
                if any(start <= match.start() < end for start, end in code_spans):
                    continue
                url = urlsplit(match[1])
                if url.scheme or url.netloc:
                    continue
                target = (page.parent / unquote(url.path)).resolve() if url.path else page
                if target.is_dir() and (target / "README.md").is_file():
                    target = target / "README.md"
                if not target.exists():
                    failures.append(f"{page.relative_to(ROOT)}:{number}: missing target {match[1]}")
                elif url.fragment and target.suffix == ".md" and unquote(url.fragment) not in anchors(target):
                    failures.append(f"{page.relative_to(ROOT)}:{number}: missing heading {match[1]}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"Checked relative links and code fences in {len(pages)} Markdown files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
