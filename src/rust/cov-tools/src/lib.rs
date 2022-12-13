use std::collections::HashMap;
use std::error::Error;
use std::path::PathBuf;

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

fn extract_to_option_pathbuf(
    py: Python<'_>,
    map: &HashMap<String, PyObject>,
    key: &str,
) -> Option<PathBuf> {
    let hash_key = String::from(key);
    let option: Option<PathBuf> = match map.get(&hash_key) {
        Some(value) => Some(value.extract::<PathBuf>(py).unwrap()),
        _ => None,
    };
    option
}

fn extract_to_default_string(
    py: Python<'_>,
    map: &HashMap<String, PyObject>,
    key: &str,
    default: &str,
) -> String {
    let hash_key = String::from(key);
    let value = match map.get(&hash_key) {
        Some(value) => value.extract::<String>(py).unwrap(),
        _ => String::from(default),
    };
    value
}

fn extract_to_bool(py: Python<'_>, map: &HashMap<String, PyObject>, key: &str) -> bool {
    let hash_key = String::from(key);
    let value = match map.get(&hash_key) {
        Some(value) => value.extract::<bool>(py).unwrap(),
        _ => false,
    };
    value
}

fn convert_hashmap_to_options(py: Python<'_>, map: HashMap<String, PyObject>) -> Options {
    let list = extract_to_option_pathbuf(py, &map, "list");
    let bam = extract_to_option_pathbuf(py, &map, "bam");
    let cram = extract_to_option_pathbuf(py, &map, "cram");
    let fasta = extract_to_option_pathbuf(py, &map, "fasta");
    let fastq1 = extract_to_option_pathbuf(py, &map, "fastq1");
    let fastq2 = extract_to_option_pathbuf(py, &map, "fastq2");
    let read_list = extract_to_option_pathbuf(py, &map, "read_list");
    let suffix = extract_to_default_string(py, &map, "suffix", "filtered");
    let fasta_out = extract_to_bool(py, &map, "fasta_out");
    let fastq_out = extract_to_bool(py, &map, "fastq_out");
    Options {
        list,
        bam,
        cram,
        fasta,
        fastq1,
        fastq2,
        suffix,
        fasta_out,
        fastq_out,
        read_list,
    }
}

#[pyfunction]
fn cov_filter_dict(py: Python<'_>, map: HashMap<String, PyObject>) -> PyResult<usize> {
    let options = &convert_hashmap_to_options(py, map);
    cov_filter(options)
}

#[pymodule]
fn cov_tools(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(cov_filter, m)?)?;
    m.add_function(wrap_pyfunction!(cov_filter_dict, m)?)?;
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
