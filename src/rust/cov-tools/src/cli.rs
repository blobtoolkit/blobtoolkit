use clap::{ArgGroup, Parser};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(group(
    ArgGroup::new("alignment")
        .required(true)
        .args(["bam", "cram"]),
))]
pub struct Config {
    /// File containing a list of sequence IDs
    #[arg(long, short = 'l', value_name = "TXT")]
    pub list: Option<PathBuf>,
    /// Path to BAM file
    #[arg(long, short = 'b')]
    pub bam: Option<PathBuf>,
    /// Path to CRAM file
    #[arg(long, short = 'c', requires = "fasta")]
    pub cram: Option<PathBuf>,
    /// Assembly FASTA file (required for CRAM)
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
    /// Output list of read IDs
    #[arg(long = "read-list", short = 'o', value_name = "TXT")]
    pub read_list: Option<PathBuf>,
}

pub fn parse() -> Config {
    Config::parse()
}
