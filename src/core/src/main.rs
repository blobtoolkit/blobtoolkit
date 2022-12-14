use std::process;

use blobtoolkit_core::cli;
use blobtoolkit_core::run;

fn main() {
    let args = cli::parse();
    if let Err(e) = run::cmd(args) {
        eprintln!("Application error: {e}");
        process::exit(1);
    }
}
