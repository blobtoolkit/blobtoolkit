extern crate atty;
use std::io::{self, BufRead, BufReader};
use std::path::PathBuf;

use std::fs::File;

fn read_stdin() -> Vec<String> {
    let stdin = io::stdin();
    let mut list: Vec<String> = vec![];
    if atty::is(atty::Stream::Stdin) {
        println!("No input on STDIN!");
        return list;
    }
    for line in stdin.lock().lines() {
        list.extend(line);
    }
    list
}

fn read_file(file_path: &PathBuf) -> Vec<String> {
    let file = File::open(file_path).expect("no such file");
    let buf = BufReader::new(file);
    buf.lines()
        .map(|line| line.expect("Could not parse line"))
        .collect()
}

pub fn get_list(file_path: &Option<PathBuf>) -> Vec<String> {
    let list = match file_path {
        None => read_stdin(),
        &Some(_) => read_file(file_path.as_ref().unwrap()),
    };
    list
}
