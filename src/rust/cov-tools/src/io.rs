extern crate atty;
use std::collections::HashSet;
use std::io::{self, BufRead, BufReader, BufWriter, Result, Write};
use std::path::{Path, PathBuf};

use std::fs::File;

use flate2::write;
use flate2::Compression;
use std::ffi::OsStr;

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

fn write_stdout(entries: &HashSet<Vec<u8>>) -> Result<()> {
    let stdout = io::stdout();
    let lock = stdout.lock();
    let mut w = BufWriter::new(lock);
    for line in entries.into_iter() {
        writeln!(&mut w, "{}", String::from_utf8(line.to_vec()).unwrap()).unwrap();
    }
    Ok(())
}

fn write_file(entries: &HashSet<Vec<u8>>, file_path: &PathBuf) -> Result<()> {
    let mut w = File::create(file_path).unwrap();
    for line in entries.into_iter() {
        writeln!(&mut w, "{}", String::from_utf8(line.to_vec()).unwrap()).unwrap();
    }
    Ok(())
}

pub fn get_list(file_path: &Option<PathBuf>) -> Vec<String> {
    let list = match file_path {
        None => read_stdin(),
        &Some(_) => read_file(file_path.as_ref().unwrap()),
    };
    list
}

pub fn write_list(entries: &HashSet<Vec<u8>>, file_path: &Option<PathBuf>) -> Result<()> {
    match &file_path {
        None => return Ok(()),
        &Some(p) if p == Path::new("-") => write_stdout(&entries),
        &Some(_) => write_file(&entries, file_path.as_ref().unwrap()),
    }?;
    Ok(())
}

pub fn get_writer(file_path: PathBuf) -> Box<dyn Write> {
    let file = match File::create(&file_path) {
        Err(why) => panic!("couldn't open {}: {}", file_path.display(), why),
        Ok(file) => file,
    };

    let writer: Box<dyn Write> = if file_path.extension() == Some(OsStr::new("gz")) {
        // Error is here: Created file isn't gzip-compressed
        Box::new(BufWriter::with_capacity(
            128 * 1024,
            write::GzEncoder::new(file, Compression::default()),
        ))
    } else {
        Box::new(BufWriter::with_capacity(128 * 1024, file))
    };
    writer
}
