#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

#[path = "scripts/elp2000/mod.rs"]
mod elp2000_build;

use std::{env, fs, path::{PathBuf, Path}};

fn write_stub_files(out_dir: &Path) {
    // Minimal VSOP87 stub with all required constants but empty data
    fn stub_vsop() -> String {
        use std::fmt::Write;
        let mut s = String::new();
        writeln!(s, "use crate::calculus::vsop87::vsop87_impl::Vsop87;").unwrap();
        let bodies = ["MERCURY","VENUS","EARTH","MARS","JUPITER","SATURN","URANUS","NEPTUNE","EMB","SUN"];
        for body in bodies.iter() {
            for coord in ["X","Y","Z"].iter() {
                for power in 0..=5 {
                    writeln!(s, "pub static {body}_{coord}{power}: [Vsop87; 0] = [];").unwrap();
                }
            }
        }
        s
    }

    // Minimal ELP2000 stub with empty series arrays
    fn stub_elp() -> String {
        use std::fmt::Write;
        let mut s = String::new();
        writeln!(s, "use crate::calculus::elp2000::elp_structs::{{MainProblem, EarthPert, PlanetPert}};").unwrap();
        for idx in 1..=36 {
            let typ = if (1..=3).contains(&idx) {
                "MainProblem"
            } else if (10..=21).contains(&idx) {
                "PlanetPert"
            } else {
                "EarthPert"
            };
            writeln!(s, "pub static ELP{idx}: &[{typ}] = &[];").unwrap();
        }
        s
    }

    fs::write(out_dir.join("vsop87a.rs"), stub_vsop()).expect("write vsop87a stub");
    fs::write(out_dir.join("vsop87e.rs"), stub_vsop()).expect("write vsop87e stub");
    fs::write(out_dir.join("elp_data.rs"), stub_elp()).expect("write elp stub");
}

fn main() {
    let out_dir = PathBuf::from(
        env::var("OUT_DIR").expect("OUT_DIR not set by Cargo"),
    );

    if env::var("SIDERUST_STUBS").is_ok() {
        write_stub_files(&out_dir);
        return;
    }

    // VSOP87
    let data_dir = out_dir.join("vsop87_dataset");
    vsop87_build::run(data_dir.as_path())
        .expect("VSOP87 codegen failed");

    // ELP2000
    let elp_dir = out_dir.join("elp2000_dataset");
    elp2000_build::run(elp_dir.as_path())
        .expect("ELP2000 codegen failed");
}
