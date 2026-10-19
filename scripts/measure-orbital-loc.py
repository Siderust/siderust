# /// script
# dependencies = ["tree-sitter==0.26.0", "tree-sitter-rust==0.24.2"]
# ///
"""Count physical Rust code lines using syntax-aware comment/test exclusion.

Usage: uv run measure.py /path/to/repository REF
Counts tracked src/**/*.rs and */src/**/*.rs, excluding test paths, cfg(test)
items/statements, comments and blank lines. Test/fixture/generated/Markdown
figures are raw physical lines, reported separately. No source is modified.
"""
import subprocess, re, json, pathlib, sys
from tree_sitter import Parser, Language
import tree_sitter_rust
parser = Parser(Language(tree_sitter_rust.language()))

def measure(repo, ref):
    def git(*args):
        return subprocess.check_output(['git','-C',str(repo),*args])
    counts = dict(production=0, tests=0, fixtures=0, generated=0, markdown=0, rust_comments=0)
    for name in git('ls-tree','-r','--name-only',ref).decode().splitlines():
        fixture = any(part in name for part in ('fixtures/', 'test-data/', 'reference_data/'))
        if not fixture and not name.endswith(('.rs','.md','.h')): continue
        data = git('show',f'{ref}:{name}')
        raw_lines = len(data.splitlines())
        if fixture:
            if b'\0' not in data: counts['fixtures'] += raw_lines
            continue
        if name.endswith('.md'):
            counts['markdown'] += raw_lines; continue
        if name.endswith('.h') or '/generated/' in name or b'@generated' in data[:500]:
            counts['generated'] += raw_lines; continue
        if not name.endswith('.rs'): continue
        if name.startswith('tests/') or '/tests/' in name or name.endswith('/tests.rs'):
            counts['tests'] += raw_lines; continue
        if not (name.startswith('src/') or '/src/' in name): continue
        tree = parser.parse(data)
        tests, comments = [], []
        def visit(node):
            pending = False
            for child in node.children:
                if child.type == 'attribute_item' and (pending or re.fullmatch(rb'#\[\s*cfg\(\s*test\s*\)\s*\]', child.text)):
                    tests.append((child.start_byte,child.end_byte)); pending=True
                elif pending:
                    tests.append((child.start_byte,child.end_byte))
                    if child.type not in ('line_comment','block_comment'): pending=False
                elif child.type in ('line_comment','block_comment'):
                    comments.append((child.start_byte,child.end_byte))
                else: visit(child)
        visit(tree.root_node)
        def occupied(ranges):
            return {data.count(b'\n',0,a)+i for a,b in ranges for i in range(len(data[a:b].splitlines()))}
        counts['tests'] += len(occupied(tests))
        counts['rust_comments'] += len(occupied(comments))
        stripped=bytearray(data)
        for a,b in tests+comments:
            for i in range(a,b):
                if stripped[i] not in (10,13): stripped[i]=32
        counts['production'] += sum(bool(line.strip()) for line in stripped.splitlines())
    return counts
if __name__ == '__main__':
    repo, ref = sys.argv[1:]
    print(json.dumps({'ref':ref, **measure(pathlib.Path(repo),ref)},indent=2))
