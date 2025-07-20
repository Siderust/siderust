import re
from pathlib import Path
from typing import List, Dict

# === PARSER ===

def parse_line_exact(line: str, n_ints: int, n_floats: int) -> Dict:
    tokens = re.findall(r"[-+]?\d+\.\d+|[-+]?\d+", line.strip())
    total = n_ints + n_floats
    if len(tokens) != total:
        raise ValueError(
            f"Expected {total} elements, found {len(tokens)} in:\n{line}"
        )
    return {
        "i": [int(tok) for tok in tokens[:n_ints]],
        "coef": tokens[n_ints:]
    }

def parse_file(path: Path, n_ints: int, n_floats: int) -> List[Dict]:
    data = []
    with path.open("r") as f:
        for line in f:
            if not re.match(r"^\s*[-+]?\d", line):
                continue
            data.append(parse_line_exact(line, n_ints, n_floats))
    return data

# Formatos por archivo
FILE_FORMATS = {
    **{f"ELP{i}": (4, 7) for i in [1, 2, 3]},
    **{f"ELP{i}": (11, 3) for i in range(10, 22)},
    **{f"ELP{i}": (5, 3) for i in list(range(4, 10)) + list(range(22, 37))},
}

def parse_all_elps(directory: Path) -> Dict[str, List[Dict]]:
    parsed_data = {}
    for path in sorted(directory.glob("ELP*")):
        key = path.stem.upper()
        if key not in FILE_FORMATS:
            print(f"⚠️  {key}: without defined format, skiping.")
            continue
        n_ints, n_floats = FILE_FORMATS[key]
        print(f"🔍 Processing {key} → {n_ints} integers, {n_floats} floats")
        parsed_data[key] = parse_file(path, n_ints, n_floats)
    return parsed_data


def generate_rust(parsed_data: Dict[str, List[Dict]], output_file="elp_data.rs"):
    with open(output_file, "w") as f:

        for name, entries in parsed_data.items():
            idx = int(name.replace("ELP", ""))
            if 1 <= idx <= 3:
                typ = "MainProblem"
            elif 4 <= idx <= 9 or 22 <= idx <= 36:
                typ = "EarthPert"
            elif 10 <= idx <= 21:
                typ = "PlanetPert"
            else:
                print(f"⚠️  {name}: index {idx} out of range 1–36, skiping.")
                continue

            f.write(f"pub static {name}: &[{typ}] = &[\n")
            for e in entries:
                coef = e["coef"]
                if typ == "MainProblem":
                    ilu = ", ".join(str(x) for x in e["i"])
                    A = coef[0]
                    B = ", ".join(coef[1:])
                    f.write(f"    {typ} {{ ilu: [{ilu}], a: {A}, b: [{B}] }},\n")

                elif typ == "EarthPert":
                    iz = float(e["i"][0])
                    ilu = ", ".join(str(x) for x in e["i"][1:5])
                    O, A, P = coef
                    f.write(
                        f"    {typ} {{ iz: {iz}, ilu: [{ilu}], o: {O}, a: {A}, p: {P} }},\n"
                    )

                elif typ == "PlanetPert":
                    ipla = ", ".join(str(x) for x in e["i"])
                    theta, O, P = coef
                    f.write(
                        f"    {typ} {{ ipla: [{ipla}], theta: {theta}, o: {O}, p: {P} }},\n"
                    )

            f.write("];\n\n")

    print(f"✅ Rust file generated with up to ELP36 structs: {output_file}")



# === MAIN ===

if __name__ == "__main__":
    import sys
    base_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    data = parse_all_elps(base_dir)
    total = sum(len(v) for v in data.values())
    print(f"\n✅ {total} obtained terms from {len(data)} files.")
    generate_rust(data, "elp_data.rs")
