use clap::Parser;
use std::process;

mod filter;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
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

//TODO: check output_bam, generate test dataset
// 1. Each read has only one primary alignment
// 2. Each read deosn't have any secondary alignment
