from pathlib import Path

HEADER = """// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

"""

# directories to skip
SKIP_DIRS = {".git", "target", "vendor", "node_modules", "qtty"}

def should_skip(path: Path) -> bool:
    return any(part in SKIP_DIRS for part in path.parts)

changed = 0
for p in Path(".").rglob("*.rs"):
    if should_skip(p):
        continue
    text = p.read_text(encoding="utf-8")

    # Don't double-insert
    if "SPDX-License-Identifier:" in text or "GNU Affero General Public License" in text:
        continue

    p.write_text(HEADER + text, encoding="utf-8")
    changed += 1

print(f"Updated {changed} file(s).")
