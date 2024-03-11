use clap::{Parser, Subcommand};
use std::process;

mod kmer_locator;
mod filter;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.3.3", about = "Refine alignments by unique kmers.", long_about = None)]
struct Arguments {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    LocateKmers {
        #[arg(short = 'i', long)]
        kmer_file: String,

        #[arg(short = 'f', long)]
        input_fasta: String,

        #[arg(short = 'k', long)]
        kmer_size: usize,
    },

    Refine {
        #[arg(short = 'i', long)]
        input_bam: String,

        #[arg(short = 'o', long)]
        output_bam: String,

        #[arg(short = 't', long)]
        hap1_tabix: String,

        #[arg(short = 'u', long)]
        hap2_tabix: String,

        #[arg(short = 'l', long)]
        hap1_list: String,

        #[arg(short = 'm', long)]
        hap2_list: String,

        #[arg(short = 'k', long)]
        kmer_size: u32,
    },
}

fn main() {
    let arguments = Arguments::parse();

    match &arguments.command {
        Commands::LocateKmers {
            kmer_file,
            input_fasta,
            kmer_size,
        } => {
            if let Err(error) = kmer_locator::run(
                kmer_file,
                input_fasta,
                *kmer_size,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        },

        Commands::Refine {
            input_bam,
            output_bam,
            hap1_tabix,
            hap2_tabix,
            hap1_list,
            hap2_list,
            kmer_size,
        } => {
            if let Err(error) = filter::run(
                input_bam,
                output_bam,
                hap1_tabix,
                hap2_tabix,
                hap1_list,
                hap2_list,
                *kmer_size,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        },
    }    
}

