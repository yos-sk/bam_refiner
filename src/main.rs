use clap::Parser;
use std::process;

mod filter;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.3.1", about = "Refine alignments by unique kmers.", long_about = None)]
struct Arguments {
    #[arg(short = 'i', long)]
    input_bam: String,

    #[arg(short = 'o', long)]
    output_bam: String,

    #[arg(short = 't', long)]
    hap1_tabix: String,

    #[arg(short = 'u', long)]
    hap2_tabix: String,

    #[arg(short = 'k', long)]
    kmer_size: u32,
}

fn main() {
    let arguments = Arguments::parse();
    if let Err(error) = filter::run(
        &arguments.input_bam,
        &arguments.output_bam,
        &arguments.hap1_tabix,
        &arguments.hap2_tabix,
        arguments.kmer_size,
    ) {
        eprintln!("{}", error);
        process::exit(1);
    }
}

