use std::{collections::HashSet, path::PathBuf};

extern crate needletail;
use needletail::{parse_fastx_file, FastxReader, Sequence};

fn open_fastx(fastx_path: &Option<PathBuf>) -> Option<Box<dyn FastxReader>> {
    let reader = match fastx_path {
        None => None,
        &Some(_) => Some(parse_fastx_file(&fastx_path.as_ref().unwrap()).expect("valid path/file")),
    };
    return reader;
}

fn trim_to_character(input: &[u8], character: u8) -> Vec<u8> {
    input
        .iter()
        .copied()
        .by_ref()
        .take_while(|&x| x != character)
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
        let seq_id: Vec<u8> = trim_to_character(seqrec.id(), b' ');
        let paired_id: Vec<u8> = trim_to_character(paired_seqrec.id(), b' ');
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

fn subsample_single(read_names: &HashSet<Vec<u8>>, mut reader: Box<dyn FastxReader>) {
    while let Some(record) = reader.next() {
        let seqrec = record.as_ref().expect("invalid record");
        println!("{:?}", String::from_utf8(seqrec.id().to_vec()));
    }
}

pub fn subsample(
    read_names: &HashSet<Vec<u8>>,
    fastq_path_1: &Option<PathBuf>,
    fastq_path_2: &Option<PathBuf>,
) -> () {
    let mut reader = open_fastx(fastq_path_1);
    let mut paired_reader = open_fastx(fastq_path_2);
    if let Some(_) = paired_reader {
        subsample_paired(read_names, reader.unwrap(), paired_reader.unwrap());
    } else if let Some(_) = reader {
        subsample_single(read_names, reader.unwrap());
    }

    // let mut paired_reader = paired_reader.unwrap();
}
