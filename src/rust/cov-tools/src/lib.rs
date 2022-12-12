use std::error::Error;

use cli::{Config, Options};

pub mod bam;
pub mod cli;
pub mod fasta;
pub mod fastq;
pub mod io;
pub mod utils;

use pyo3::prelude::*;

/// runs python cov_tools.
#[pyfunction]
fn cov_filter(options: &Options) -> PyResult<usize> {
    let seq_names = io::get_list(&options.list);
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
    Ok(read_names.len())
}

#[pymodule]
fn cov_tools(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cov_filter, m)?)?;
    m.add_class::<Options>()?;
    Ok(())
}

pub fn run(options: Config) -> Result<(), Box<dyn Error>> {
    let seq_names = io::get_list(&options.list);
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
