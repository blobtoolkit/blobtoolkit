use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};

extern crate needletail;
use needletail::parser::{write_fastq, LineEnding};
use needletail::{parse_fastx_file, FastxReader};

use crate::io::get_writer;

fn open_fastx(fastx_path: &Option<PathBuf>) -> Option<Box<dyn FastxReader>> {
    let reader = match fastx_path {
        None => None,
        &Some(_) => Some(parse_fastx_file(&fastx_path.as_ref().unwrap()).expect("valid path/file")),
    };
    return reader;
}

fn trim_read_id(input: &[u8]) -> Vec<u8> {
    input
        .iter()
        .copied()
        .by_ref()
        .take_while(|&x| x != b' ' && x != b'/')
        .collect()
}

fn subsample_paired(
    read_names: &HashSet<Vec<u8>>,
    mut reader: Box<dyn FastxReader>,
    mut paired_reader: Box<dyn FastxReader>,
) {
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let paired_record = paired_reader.next().unwrap();
        let paired_seqrec = paired_record.expect("invalid paired record");
        let seq_id: Vec<u8> = trim_read_id(seqrec.id());
        let paired_id: Vec<u8> = trim_read_id(paired_seqrec.id());
        // println!("{:?}", String::from_utf8(&seq_id));
        // println!("{:?}", String::from_utf8(&paired_id));
        if seq_id != paired_id {
            // println!(
            //     "Fasta files not sorted consistently {:?} ne {:?}",
            //     String::from_utf8(seq_id),
            //     String::from_utf8(paired_id)
            // );
        }
        if read_names.contains(&seq_id) || read_names.contains(&paired_id) {
            println!("{:?}", String::from_utf8(seq_id));
        }
    }
}

fn subsample_single(
    read_names: &HashSet<Vec<u8>>,
    mut reader: Box<dyn FastxReader>,
    writer: &mut dyn Write,
) {
    while let Some(record) = reader.next() {
        let seqrec = record.as_ref().expect("invalid record");
        let seq_id: Vec<u8> = trim_read_id(seqrec.id());
        if read_names.contains(&seq_id) {
            write_fastq(
                &seq_id as &[u8],
                &seqrec.seq(),
                seqrec.qual(),
                writer,
                LineEnding::Unix,
            )
            .expect("Unable to write FASTQ")
        }
    }
}

fn suffix_file_name(path: impl AsRef<Path>, suffix: &String) -> PathBuf {
    let path = path.as_ref();
    let mut result = path.to_owned();
    let mut new_name = String::from("_tmp");
    let mut extension = String::from("fastq");
    if let Some(ext) = path.extension() {
        extension = String::from(ext.to_str().unwrap());
        if extension == "gz" {
            let file_path = PathBuf::from(path.file_stem().unwrap().to_str().unwrap());
            if let Some(pre_ext) = file_path.extension() {
                extension = format!("{}.{}", pre_ext.to_str().unwrap(), extension);

                new_name = format!(
                    "{}.{}",
                    file_path.file_stem().unwrap().to_str().unwrap(),
                    suffix
                );
            } else {
                new_name = format!("{}.{}", path.file_stem().unwrap().to_str().unwrap(), suffix);
            }
        } else {
            new_name = format!("{}.{}", path.file_stem().unwrap().to_str().unwrap(), suffix);
        }
    }
    result.set_extension(&extension);
    result.set_file_name(format!("{}.{}", &new_name, &extension));
    result
}

pub fn subsample(
    read_names: &HashSet<Vec<u8>>,
    fastq_path_1: &Option<PathBuf>,
    fastq_path_2: &Option<PathBuf>,
    suffix: &String,
) -> () {
    let reader = open_fastx(fastq_path_1);
    let paired_reader = open_fastx(fastq_path_2);

    let out_path = suffix_file_name(fastq_path_1.as_ref().unwrap(), &suffix);
    let mut writer = get_writer(out_path);
    // let mut writer: Box<&mut dyn Write> = get_writer(PathBuf::from("filtered.reads_1.fq"));
    if let Some(_) = paired_reader {
        subsample_paired(read_names, reader.unwrap(), paired_reader.unwrap());
    } else if let Some(_) = reader {
        subsample_single(read_names, reader.unwrap(), &mut *writer);
    }

    // let mut paired_reader = paired_reader.unwrap();
}
