use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};

extern crate needletail;
use needletail::parser::{write_fastq, LineEnding};
use needletail::{parse_fastx_file, FastxReader};

use crate::io::get_writer;
use crate::utils::styled_progress_bar;

pub fn open_fastx(fastx_path: &Option<PathBuf>) -> Option<Box<dyn FastxReader>> {
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
    writer: &mut dyn Write,
    paired_writer: &mut dyn Write,
    read_suffix: &[Vec<u8>; 2],
) {
    let total = read_names.len() as u64;
    let progress_bar = styled_progress_bar(total, "Subsampling FASTQ");

    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        let paired_record = paired_reader.next().unwrap();
        let paired_seqrec = paired_record.expect("invalid paired record");
        let mut seq_id: Vec<u8> = trim_read_id(seqrec.id());
        let mut paired_id: Vec<u8> = trim_read_id(paired_seqrec.id());
        if seq_id != paired_id {
            panic!(
                "Fasta files not sorted consistently {:?} ne {:?}",
                String::from_utf8(seq_id),
                String::from_utf8(paired_id)
            );
        }
        if read_names.contains(&seq_id) || read_names.contains(&paired_id) {
            seq_id.extend(&read_suffix[0]);
            write_fastq(
                &seqrec.id(),
                &seqrec.seq(),
                seqrec.qual(),
                writer,
                LineEnding::Unix,
            )
            .expect("Unable to write FASTQ");
            paired_id.extend(&read_suffix[1]);
            write_fastq(
                &paired_seqrec.id(),
                &paired_seqrec.seq(),
                paired_seqrec.qual(),
                paired_writer,
                LineEnding::Unix,
            )
            .expect("Unable to write FASTQ");
            progress_bar.inc(1);
            if progress_bar.position() == total {
                break;
            }
        }
    }
    progress_bar.finish();
}

fn subsample_single(
    read_names: &HashSet<Vec<u8>>,
    mut reader: Box<dyn FastxReader>,
    writer: &mut dyn Write,
    read_suffix: &[Vec<u8>; 2],
) {
    let total = read_names.len() as u64;
    let progress_bar = styled_progress_bar(total, "Subsampling FASTQ");

    while let Some(record) = reader.next() {
        let seqrec = record.as_ref().expect("invalid record");
        let mut seq_id: Vec<u8> = trim_read_id(seqrec.id());
        if read_names.contains(&seq_id) {
            seq_id.extend(&read_suffix[0]);
            write_fastq(
                &seqrec.id(),
                &seqrec.seq(),
                seqrec.qual(),
                writer,
                LineEnding::Unix,
            )
            .expect("Unable to write FASTQ");
            progress_bar.inc(1);
            if progress_bar.position() == total {
                break;
            }
        }
    }
    progress_bar.finish();
}

pub fn suffix_file_name(path: impl AsRef<Path>, suffix: &String) -> PathBuf {
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

fn set_read_suffix(read_names: &HashSet<Vec<u8>>) -> [Vec<u8>; 2] {
    let first_name = read_names.iter().next().unwrap();
    if first_name.contains(&b'/') {
        return [vec![] as Vec<u8>, vec![] as Vec<u8>];
    }
    [vec![b'/', b'1'], vec![b'/', b'2']]
}

pub fn subsample(
    read_names: &HashSet<Vec<u8>>,
    fastq_path_1: &Option<PathBuf>,
    fastq_path_2: &Option<PathBuf>,
    fastq_out: &bool,
    suffix: &String,
) -> () {
    if let None = fastq_path_1 {
        return;
    }
    if !fastq_out {
        return;
    }
    let reader = open_fastx(fastq_path_1);
    let paired_reader = open_fastx(fastq_path_2);
    let read_suffix = set_read_suffix(read_names);
    let out_path = suffix_file_name(fastq_path_1.as_ref().unwrap(), &suffix);
    let mut writer = get_writer(out_path);
    if let Some(_) = paired_reader {
        let paired_out_path = suffix_file_name(fastq_path_2.as_ref().unwrap(), &suffix);
        let mut paired_writer = get_writer(paired_out_path);
        subsample_paired(
            read_names,
            reader.unwrap(),
            paired_reader.unwrap(),
            &mut *writer,
            &mut *paired_writer,
            &read_suffix,
        );
    } else if let Some(_) = reader {
        subsample_single(read_names, reader.unwrap(), &mut *writer, &read_suffix);
    }
}
