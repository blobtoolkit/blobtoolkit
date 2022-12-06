use std::error::Error;

use cli::Config;

pub mod bam;
pub mod cli;
pub mod fastq;
pub mod io;

pub fn run(options: Config) -> Result<(), Box<dyn Error>> {
    let seq_names = io::get_list(&options.list);
    let bam = bam::open_bam(&options.bam, &options.cram, &options.fasta);
    let read_names = bam::reads_from_bam(&seq_names, bam);
    println!("{:#?}", read_names);
    fastq::subsample(&read_names, &options.fastq1, &options.fastq2);
    // let entries = io::process_fasta(&options);
    // let assembly_stats = assembly::assembly_stats(entries);
    // let snail_stats = snail::snail_stats(&options, assembly_stats);
    // io::write_snail_json(&options, snail_stats)?;
    Ok(())
}
