mod kmer_locator;
mod kmer_ratio;
mod refine;
mod single;

use clap::{Parser, Subcommand};
use std::process;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.3.6", about = "Refine alignments by unique kmers.", long_about = None)]
struct Arguments {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    KmerRatio {
        #[clap(value_parser, default_value = "-")]
        input_bam: String,

        #[clap(short, long, value_parser, default_value_t = 4)]
        threads: usize,
    },

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

        #[clap(short = 'p', long, value_parser, default_value_t = 4)]
        threads: usize,
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

        #[clap(short = 'p', long, value_parser, default_value_t = 4)]
        threads: usize,
    },
}

fn main() {
    let arguments = Arguments::parse();

    match &arguments.command {
        Commands::KmerRatio {
            input_bam,
            threads,
        } => {
            if let Err(error) = kmer_ratio::run(
                input_bam,
                *threads,
            ) {
                eprintln!("{}", error);
                process::exit(1); 
            }
        },
        
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
            threads,
        } => {
            if let Err(error) = refine::run(
                input_bam,
                output_bam,
                hap1_tabix,
                hap2_tabix,
                hap1_list,
                hap2_list,
                *kmer_size,
                *threads,
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
            threads,
        } => {
            if let Err(error) = single::run(
                input_bam,
                output_bam,
                ref_tabix,
                *kmer_size,
                *threads,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        },

    }    
}

