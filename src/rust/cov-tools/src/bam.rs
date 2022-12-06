use rust_htslib::bam::{index, IndexedReader, Read};
use rust_htslib::htslib;
use std::collections::HashSet;
use std::path::PathBuf;

pub fn create_index(bam_path: &PathBuf) -> () {
    match index::build(bam_path, None, index::Type::Bai, 1) {
        Err(e) => println!("Error writing BAM index: {e:?}"),
        Ok(_) => println!("Successfully created BAM index"),
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

pub fn reads_from_bam(seq_names: &Vec<String>, mut bam: IndexedReader) -> HashSet<Vec<u8>> {
    let mut wanted_reads = HashSet::new();

    for seq_name in seq_names {
        match bam.fetch(seq_name) {
            Err(_) => println!("Sequence name not found in BAM file: {:?}", seq_name),
            Ok(_) => println!("Found {:?}", seq_name),
        }

        for read in bam
            .rc_records()
            .map(|x| x.expect("Failure parsing Bam file"))
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
            // let read_name = str::from_utf8(read.qname()).unwrap();
            // println!("Found a forward read: {:?}", read_name);
            // println!("Mapped to contig: {:?}", read.tid());
        }
    }
    wanted_reads
    // .into_iter()
    // .map(|x| String::from_utf8(x).unwrap())
    // .collect()
}
