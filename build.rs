#[path = "scripts/vsop87/mod.rs"]
mod vsop87_build;

fn main() {
    vsop87_build::run("dataset")
        .expect("VSOP87 codegen failed");
}
