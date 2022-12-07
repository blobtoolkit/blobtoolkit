use rust_htslib::bam::{index, Header, IndexedReader, Read};
use rust_htslib::htslib; // for BAM_F*
use std::collections::HashSet;
// use std::str;
fn main() {
    let bam_path = "test/test.bam";
    match index::build(bam_path, None, index::Type::Bai, 1) {
        Err(e) => println!("Error writing BAM index: {e:?}"),
        Ok(_) => println!("Successfully created BAM index"),
    }

    let mut bam = IndexedReader::from_path(&bam_path).unwrap();

    let header = Header::from_template(bam.header());

    // let mut all_seqs = HashSet::new();

    // let mut seq_id_map = HashSet::new();

    for (key, records) in header.to_hashmap() {
        // for record in records.iter().filter(|r| r.contains_key("SN")) {
        for record in records {
            if record.contains_key("SN") {
                println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
                // all_seqs.insert(record["SN"].to_string());
                // seq_id_map.insert(record["SN"].to_string(), record.)
            }
        }
    }
    // println!("List of sequence names: {:?}", all_seqs);

    let mut wanted_seqs = HashSet::new();
    wanted_seqs.insert("FJNM01000076.1");
    wanted_seqs.insert("FJNM01002842.1");

    // println!("List of wanted names: {:?}", wanted_seqs);

    // let retained_seqs = all_seqs.intersection(&wanted_seqs);

    // println!("List of sequences to keep: {:?}", retained_seqs);

    let mut wanted_reads = HashSet::new();

    for seq_name in wanted_seqs {
        match bam.fetch(seq_name) {
            Err(e) => println!("Sequence name not found in BAM file: {:?}", seq_name),
            Ok(_) => println!("Found {:?}", seq_name),
        }
        // for read in bam.records() {
        //     println!("read name: {:?}", read.unwrap().qname());
        // }

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

    println!("{:#?}", wanted_reads);
}
