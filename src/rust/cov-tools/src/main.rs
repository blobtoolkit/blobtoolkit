use std::process;

use cov_tools;
use cov_tools::cli;

fn main() {
    let options = cli::parse();
    if let Err(e) = cov_tools::run(options) {
        println!("Application error: {e}");
        process::exit(1);
    }
}
