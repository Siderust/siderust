import re
from pathlib import Path
from typing import List, Dict

# === PARSER ===

def parse_line_exact(line: str, n_ints: int, n_floats: int) -> Dict:
    tokens = re.findall(r"[-+]?\d+\.\d+|[-+]?\d+", line.strip())
    total = n_ints + n_floats
    if len(tokens) != total:
        raise ValueError(
            f"Esperados {total} elementos, encontrados {len(tokens)} en:\n{line}"
        )
    return {
        "i": [int(tok) for tok in tokens[:n_ints]],
        "coef": tokens[n_ints:]  # mantener precisión exacta
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
            print(f"⚠️  {key}: sin formato definido, se omite.")
            continue
        n_ints, n_floats = FILE_FORMATS[key]
        print(f"🔍 Procesando {key} → {n_ints} enteros, {n_floats} flotantes")
        parsed_data[key] = parse_file(path, n_ints, n_floats)
    return parsed_data

# === GENERADOR RUST ===

def generate_rust(parsed_data: Dict[str, List[Dict]], output_file="elp_data.rs"):
    with open(output_file, "w") as f:
        # … definición de los tres structs igual que antes …

        for name, entries in parsed_data.items():
            idx = int(name.replace("ELP", ""))
            if 1 <= idx <= 3:
                typ = "MainProblem"
            # ahora earth_pert para 4–9 **y** 22–36
            elif 4 <= idx <= 9 or 22 <= idx <= 36:
                typ = "EarthPert"
            elif 10 <= idx <= 21:
                typ = "PlanetPert"
            else:
                print(f"⚠️  {name}: índice {idx} fuera de rango 1–36, se omite.")
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
                    # primer entero → iz, siguientes 4 → ilu
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

    print(f"✅ Archivo Rust generado con structs hasta el ELP36: {output_file}")



# === MAIN ===

if __name__ == "__main__":
    import sys
    base_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
    data = parse_all_elps(base_dir)
    total = sum(len(v) for v in data.values())
    print(f"\n✅ Se obtuvieron {total} términos de {len(data)} archivos.")
    generate_rust(data, "elp_data.rs")
