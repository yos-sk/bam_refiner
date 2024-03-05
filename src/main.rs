use clap::{Parser, Subcommand};
use std::process;

mod kmer_locator;
mod collect_alignments;
mod kmer_count;
mod refine;
mod write_bam;

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

        #[arg(short = 't', long)]
        threads: usize,

        #[arg(short = 'o', long)]
        output_bam: String,

        #[arg(short = 'a', long)]
        hap1_tabix: String,

        #[arg(short = 'b', long)]
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
            threads,
            output_bam,
            hap1_tabix,
            hap2_tabix,
            hap1_list,
            hap2_list,
            kmer_size,
        } => {
            let (alignments, sequences, qualities, indices) =
                collect_alignments::run(input_bam, *threads)
                    .expect("Error: collect_alignments step");
            

            let kmer_counted_data = kmer_count::run(
                alignments,
                sequences.clone(),
                hap1_tabix,
                hap2_tabix,
                hap1_list,
                hap2_list,
                *kmer_size,
                *threads,
            )
            .expect("Error: kmer_count step");
            let refined_alignments = refine::run(kmer_counted_data, indices).expect("Error: refine step");
            let _ = write_bam::run(
                input_bam,
                output_bam,
                &refined_alignments,
                &sequences,
                &qualities,
                hap1_list,
                hap2_list,
            );
        }
    }
}
