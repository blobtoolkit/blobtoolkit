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
    io::write_list(&read_names, &options.read_list)?;
    fastq::subsample(
        &read_names,
        &options.fastq1,
        &options.fastq2,
        &options.suffix,
    );
    // fasta::subsample;
    Ok(())
}
