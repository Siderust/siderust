use std::{
    collections::BTreeMap,
    env,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::PathBuf,
};

use anyhow::Context;
use regex::Regex;
use walkdir::WalkDir;

/// Representa un término Vsop87 final a·cos(b + c·T)
#[derive(Clone, Copy)]
struct Term {
    a: f64,
    b: f64,
    c: f64,
}

// Tipos auxiliares para el árbol planet → coord → t_power
//
// PlanetMap        = BTreeMap< String  , CoordMap   >
// CoordMap         = BTreeMap<    u8    , TPowerMap >
// TPowerMap        = BTreeMap<    u8    , Vec<Term>>
//
// En esta versión mantenemos un nivel adicional para agrupar por versión
// (VSOP87A, VSOP87E, ...).
type CoordMap = BTreeMap<u8, TPowerMap>;
type TPowerMap = BTreeMap<u8, Vec<Term>>;
type PlanetMap = BTreeMap<String, CoordMap>;
type VersionMap = BTreeMap<char, PlanetMap>; // ‘A’, ‘E’, …

/// Formatea flotantes de modo que los enteros salgan como “0.0”
fn fmt_f(v: f64) -> String {
    if (v - v.round()).abs() < 1e-15 {
        format!("{v:.1}")
    } else {
        // 14 dígitos significativos como en el script original
        format!("{v:.14}")
    }
}

/// Convierte S·sin() + K·cos() en un único término coseno,
/// a menos que se haya dado A directamente.
fn sk_to_term(s: f64, k: f64, a: f64, b: f64, c: f64) -> Option<Term> {
    if a.abs() > 1e-15 {
        Some(Term { a, b, c })
    } else if s.abs() > 1e-15 || k.abs() > 1e-15 {
        let r = (s * s + k * k).sqrt();
        let alpha = s.atan2(k); // atan2(s, k)
        Some(Term {
            a: r,
            b: b - alpha,
            c,
        })
    } else {
        None
    }
}

/// Parsea una línea de datos y devuelve (S, K, A, B, C)
fn parse_data_line(line: &str) -> Option<(f64, f64, f64, f64, f64)> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 5 {
        return None;
    }
    let len = parts.len();
    let c = parts[len - 1].parse().ok()?;
    let b = parts[len - 2].parse().ok()?;
    let a = parts[len - 3].parse().ok()?;
    let k = parts[len - 4].parse().ok()?;
    let s = parts[len - 5].parse().ok()?;
    Some((s, k, a, b, c))
}

// ---------------------------------------------------------------------
// Configuración
// ---------------------------------------------------------------------
const DATA_DIR: &str = "dataset";
const REGEX_FILE: &str = r"^VSOP87([A-Z])\.[A-Za-z]{3}$"; // capturamos la versión
const HEADER_REGEX: &str =
    r"VSOP87 VERSION [A-Z]?\d+\s+(\S+)\s+VARIABLE\s+(\d+)\s+\(XYZ\)\s+\*T\*\*(\d+)\s+(\d+)\s+TERMS";

fn main() -> anyhow::Result<()> {
    // ---------- Configuración de expresiones regulares ----------
    let vsop87_path = PathBuf::from(DATA_DIR);
    let file_re = Regex::new(REGEX_FILE)?;
    let header_re = Regex::new(HEADER_REGEX)?;

    // Vuelve a compilar si cambian los datos
    println!("cargo:rerun-if-changed={}", vsop87_path.display());

    // Directorio de salida del build script
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Tabla agregada por versión → planeta → coordenada → potencia
    let mut versions: VersionMap = BTreeMap::new();

    // -----------------------------------------------------------------
    // 1. Recorremos todos los ficheros y acumulamos los términos
    // -----------------------------------------------------------------
    for entry in WalkDir::new(&vsop87_path).min_depth(1).into_iter().filter_map(Result::ok) {
        let path = entry.path();
        let fname = path.file_name().unwrap().to_string_lossy();
        
        // Comprobamos si es un fichero VSOP87X.xxx
        let caps = match file_re.captures(&fname) {
            Some(c) => c,
            None => continue, // no coincide
        };

        // Carácter de versión (A, E, ...)
        let v_char = caps[1].chars().next().unwrap();

        // Abrimos fichero
        let file = File::open(path)
            .with_context(|| format!("No se pudo abrir {path:?}"))?;
        let reader = BufReader::new(file);

        // Estado de parseo
        let mut cur_planet = String::new();
        let mut cur_coord: u8 = 0;
        let mut cur_t: u8 = 0;
        let mut remaining = 0u32;

        for line in reader.lines() {
            let line = line.expect("Error leyendo línea");

            // ¿Nueva cabecera?
            if let Some(cap) = header_re.captures(&line) {
                cur_planet = cap[1].to_uppercase();
                cur_coord = cap[2].parse::<u8>().unwrap(); // 1 = X, 2 = Y, 3 = Z
                cur_t = cap[3].parse::<u8>().unwrap(); // potencia de T
                remaining = cap[4].parse::<u32>().unwrap();
                continue;
            }

            // Dentro de un bloque, parsear términos
            if remaining > 0 {
                if let Some((s, k, a, b, c)) = parse_data_line(&line) {
                    if let Some(term) = sk_to_term(s, k, a, b, c) {
                        versions
                            .entry(v_char)
                            .or_default() // PlanetMap
                            .entry(cur_planet.clone())
                            .or_default() // CoordMap
                            .entry(cur_coord)
                            .or_default() // TPowerMap
                            .entry(cur_t)
                            .or_default()
                            .push(term);
                    }
                    remaining -= 1;
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // 2. Generamos un módulo .rs por cada versión
    // -----------------------------------------------------------------
    let coord_letter = |c| match c {
        1 => 'X',
        2 => 'Y',
        3 => 'Z',
        _ => '?',
    };

    for (version, planets) in versions {
        let mut code = String::new();
        code.push_str("// -----------------------------------------------\n");
        code.push_str("// Generado automáticamente por build.rs\n");
        code.push_str("// -----------------------------------------------\n");
        code.push_str("use crate::calculus::vsop87::Vsop87;\n\n");

        for (planet, coords) in &planets {
            let planet_up = planet.to_uppercase();
            for (coord, t_powers) in coords {
                let letter = coord_letter(*coord);
                for (t, terms) in t_powers {
                    let array_name = format!("{planet_up}_{letter}{t}");
                    code.push_str(&format!(
                        "pub static {array_name}: [Vsop87; {}] = [\n",
                        terms.len()
                    ));
                    for Term { a, b, c } in terms {
                        code.push_str(&format!(
                            "    Vsop87 {{ a: {}, b: {}, c: {} }},\n",
                            fmt_f(*a),
                            fmt_f(*b),
                            fmt_f(*c)
                        ));
                    }
                    code.push_str("];\n\n");
                }
            }
        }

        // Nombre del fichero de salida: vsop87a.rs, vsop87e.rs, ...
        let file_name = format!("vsop87{}.rs", version.to_ascii_lowercase());
        let out_file = out_dir.join(file_name);

        let mut f = File::create(&out_file)
            .unwrap_or_else(|e| panic!("No se pudo crear {out_file:?}: {e}"));
        f.write_all(code.as_bytes())
            .unwrap_or_else(|e| panic!("Error escribiendo {out_file:?}: {e}"));
    }

    Ok(())
}
