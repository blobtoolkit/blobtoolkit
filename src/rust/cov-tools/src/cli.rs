use clap::{ArgGroup, Parser};
use pyo3::prelude::*;
use std::path::PathBuf;

#[derive(Debug)]
#[pyclass]
pub struct Options {
    /// File containing a list of sequence IDs
    pub list: Option<PathBuf>,
    /// Path to BAM file
    pub bam: Option<PathBuf>,
    /// Path to CRAM file
    pub cram: Option<PathBuf>,
    /// Path to assembly FASTA input file (required for CRAM)
    pub fasta: Option<PathBuf>,
    /// Path to FASTQ file to filter (forward or single reads)
    pub fastq1: Option<PathBuf>,
    /// Path to paired FASTQ file to filter (reverse reads)
    pub fastq2: Option<PathBuf>,
    /// Suffix to use for output filtered files
    pub suffix: String,
    /// Flag to output a filtered FASTA file
    pub fasta_out: bool,
    /// Flag to output filtered FASTQ files
    pub fastq_out: bool,
    /// Path to output list of read IDs
    pub read_list: Option<PathBuf>,
}

#[pymethods]
impl Options {
    #[new]
    fn new(
        list: Option<PathBuf>,
        bam: Option<PathBuf>,
        cram: Option<PathBuf>,
        fasta: Option<PathBuf>,
        fastq1: Option<PathBuf>,
        fastq2: Option<PathBuf>,
        suffix: String,
        fasta_out: bool,
        fastq_out: bool,
        read_list: Option<PathBuf>,
    ) -> Self {
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
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(group(
    ArgGroup::new("alignment")
        .required(true)
        .args(["bam", "cram"]),
))]
pub struct Config {
    /// File containing a list of sequence IDs
    // TODO: add option to invert list (use BAM header)
    #[arg(long, short = 'l', value_name = "TXT")]
    pub list: Option<PathBuf>,
    /// Path to BAM file
    #[arg(long, short = 'b')]
    pub bam: Option<PathBuf>,
    /// Path to CRAM file
    #[arg(long, short = 'c', requires = "fasta")]
    pub cram: Option<PathBuf>,
    /// Path to assembly FASTA input file (required for CRAM)
    #[arg(long, short = 'a')]
    pub fasta: Option<PathBuf>,
    /// Path to FASTQ file to filter (forward or single reads)
    #[arg(long = "fastq", short = 'f', value_name = "FASTQ")]
    pub fastq1: Option<PathBuf>,
    /// Path to paired FASTQ file to filter (reverse reads)
    #[arg(
        long = "fastq2",
        short = 'r',
        value_name = "FASTQ",
        requires = "fastq1"
    )]
    pub fastq2: Option<PathBuf>,
    /// Suffix to use for output filtered files
    #[arg(long, short = 'S', value_name = "SUFFIX", default_value_t = String::from("filtered"))]
    pub suffix: String,
    /// Flag to output a filtered FASTA file
    #[arg(
        long = "fasta-out",
        short = 'A',
        requires = "fasta",
        default_value_t = false
    )]
    pub fasta_out: bool,
    /// Flag to output filtered FASTQ files
    #[arg(
        long = "fastq-out",
        short = 'F',
        requires = "fastq1",
        default_value_t = false
    )]
    pub fastq_out: bool,
    /// Path to output list of read IDs
    #[arg(long = "read-list", short = 'O', value_name = "TXT")]
    pub read_list: Option<PathBuf>,
}

pub fn parse() -> Config {
    Config::parse()
}
