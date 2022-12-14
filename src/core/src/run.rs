use std::error::Error;

use crate::bam;
use crate::cli;
use crate::fasta;
use crate::fastq;
use crate::io;

pub fn filter(options: &cli::FilterOptions) -> Result<(), Box<dyn Error>> {
    let seq_names = io::get_list(&options.list_file);
    if seq_names.len() == 0 {
        return Ok(());
    }
    fasta::subsample(
        &seq_names,
        &options.fasta,
        &options.fasta_out,
        &options.suffix,
    );
    let bam = bam::open_bam(&options.bam, &options.cram, &options.fasta);
    let read_names = bam::reads_from_bam(&seq_names, bam);
    io::write_list(&read_names, &options.read_list)?;
    fastq::subsample(
        &read_names,
        &options.fastq1,
        &options.fastq2,
        &options.fastq_out,
        &options.suffix,
    );
    Ok(())
}

pub fn cmd(args: cli::Arguments) -> Result<(), Box<dyn Error>> {
    match args.cmd {
        cli::SubCommand::Filter(options) => filter(&options),
    }
}
