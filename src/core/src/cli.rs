use clap::{ArgGroup, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]

pub struct Arguments {
    #[clap(subcommand)]
    pub cmd: SubCommand,
}

#[derive(Subcommand, Debug)]
pub enum SubCommand {
    /// Filter files based on list of sequence names
    Filter(FilterOptions),
}

#[derive(Parser, Debug)]
#[command(group(
    ArgGroup::new("alignment")
        .required(true)
        .args(["bam", "cram"]),
))]
pub struct FilterOptions {
    /// Path to input file containing a list of sequence IDs
    // TODO: add option to invert list (use BAM header)
    #[arg(long = "list", short = 'i', value_name = "TXT")]
    pub list_file: Option<PathBuf>,
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

pub fn parse() -> Arguments {
    Arguments::parse()
}
