use clap::{Parser, Subcommand};
use std::process;

mod filter;
mod single;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.3.1", about = "Refine alignments by unique kmers.", long_about = None)]

struct Arguments {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Filter {
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
    },
    Single {
        #[arg(short = 'i', long)]
        input_bam: String,

        #[arg(short = 'o', long)]
        output_bam: String,

        #[arg(short = 't', long)]
        ref_tabix: String,

        #[arg(short = 'k', long)]
        kmer_size: u32,
    },
}


fn main() {
    let arguments = Arguments::parse();

    match &arguments.command {
        Commands::Filter {
            input_bam, 
            output_bam,
            hap1_tabix,
            hap2_tabix,
            kmer_size, 
        } => {
            if let Err(error) = filter::run(
                input_bam,
                output_bam,
                hap1_tabix,
                hap2_tabix,
                *kmer_size,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        },

        Commands::Single {
            input_bam, 
            output_bam,
            ref_tabix,
            kmer_size,
        } => {
            if let Err(error) = single::run(
                input_bam,
                output_bam,
                ref_tabix,
                *kmer_size,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        },
    }
    
}

//TODO: check output_bam, generate test dataset
// 1. Each read has only one primary alignment
// 2. Each read deosn't have any secondary alignment
