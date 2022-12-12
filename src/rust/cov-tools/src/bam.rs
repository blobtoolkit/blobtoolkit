use crate::utils::styled_progress_bar;
use rust_htslib::bam::{index, IndexedReader, Read};
use rust_htslib::htslib;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

fn add_extension(path: &mut PathBuf, extension: impl AsRef<Path>) {
    match path.extension() {
        Some(ext) => {
            let mut ext = ext.to_os_string();
            ext.push(".");
            ext.push(extension.as_ref());
            path.set_extension(ext)
        }
        None => path.set_extension(extension.as_ref()),
    };
}

pub fn create_index(bam_path: &PathBuf) -> () {
    let mut csi = PathBuf::from(bam_path);
    add_extension(&mut csi, "csi");
    if Path::new(&csi).exists() {
        return;
    }
    let mut bai = PathBuf::from(bam_path);
    add_extension(&mut bai, "bai");
    if Path::new(&bai).exists() {
        return;
    }
    match index::build(bam_path, None, index::Type::Csi(14), 1) {
        Err(e) => eprintln!("Error writing BAM index: {e:?}"),
        Ok(_) => eprintln!("Successfully created BAM index"),
    }
}

pub fn open_bam(
    bam_path: &Option<PathBuf>,
    cram_path: &Option<PathBuf>,
    _fasta_path: &Option<PathBuf>,
) -> IndexedReader {
    let bam_cram_path = match bam_path {
        None => cram_path.as_ref().unwrap(),
        &Some(_) => bam_path.as_ref().unwrap(),
    };
    create_index(&bam_cram_path);
    let bam = IndexedReader::from_path(&bam_cram_path).unwrap();
    bam
}

pub fn reads_from_bam(seq_names: &HashSet<Vec<u8>>, mut bam: IndexedReader) -> HashSet<Vec<u8>> {
    let mut wanted_reads = HashSet::new();
    let total = seq_names.len() as u64;
    let progress_bar = styled_progress_bar(total, "Locating alignments");

    for seq_name in seq_names {
        match bam.fetch(seq_name) {
            Err(_) => eprintln!("Sequence {:?} not found in BAM file", seq_name),
            Ok(_) => (),
        }

        for read in bam
            .rc_records()
            .map(|x| x.expect("Failure parsing Bam file"))
            // TODO: include filter options in config
            .filter(|read| {
                read.flags()
                    & (htslib::BAM_FUNMAP
                        | htslib::BAM_FSECONDARY
                        | htslib::BAM_FQCFAIL
                        | htslib::BAM_FDUP) as u16
                    == 0
            })
        {
            wanted_reads.insert(read.qname().to_vec());
            // TODO: add option to print info from matching records file, e.g.
            // println!(
            //     "{:?}: {:?}",
            //     String::from_utf8(read.qname().to_vec()).unwrap(),
            //     read.cigar().to_string()
            // );
        }
        progress_bar.inc(1);
    }
    progress_bar.finish();
    wanted_reads
}
