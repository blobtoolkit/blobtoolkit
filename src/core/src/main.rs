use std::process;

use blobtoolkit_core;
use blobtoolkit_core::cli;

fn main() {
    let options = cli::parse();
    if let Err(e) = blobtoolkit_core::run(options) {
        eprintln!("Application error: {e}");
        process::exit(1);
    }
}
