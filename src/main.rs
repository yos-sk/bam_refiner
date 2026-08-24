mod kmer_locator;
mod kmer_ratio;
mod refine;
mod single;

use clap::{Parser, Subcommand};
use std::process;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.4.0", about = "Refine alignments by unique kmers.", long_about = None)]
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

        /// Ratio reported for a read that spans no haplotype-specific locus,
        /// and the value small denominators are pulled toward.
        /// Suggested: 0.8 for HiFi, 0.6 for ONT.
        #[clap(short = 'm', long, value_parser, default_value_t = 0.5)]
        prior_mean: f64,

        /// Strength of that prior, in pseudo-loci. 0 restores the raw quotient.
        #[clap(short = 'w', long, value_parser, default_value_t = 1.0)]
        prior_weight: f64,
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

        /// Adopt the best placement only if
        /// max_kmer / (max_kmer + second_max_kmer) >= this value.
        /// 0.5 restores the pre-0.4.0 behaviour (any margin wins).
        #[clap(short = 'r', long, value_parser, default_value_t = 0.8)]
        ratio_threshold: f64,

        /// Leave a placement undetermined when it matched fewer than this many
        /// haplotype-specific loci AND left some of the loci available to it
        /// unmatched. 1 disables the rule.
        #[clap(short = 'n', long, value_parser, default_value_t = 3)]
        min_markers: usize,
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

        /// Adopt the best placement only if
        /// max_kmer / (max_kmer + second_max_kmer) >= this value.
        /// 0.5 restores the pre-0.4.0 behaviour (any margin wins).
        #[clap(short = 'r', long, value_parser, default_value_t = 0.8)]
        ratio_threshold: f64,

        /// Leave a placement undetermined when it matched fewer than this many
        /// haplotype-specific loci AND left some of the loci available to it
        /// unmatched. 1 disables the rule.
        #[clap(short = 'n', long, value_parser, default_value_t = 3)]
        min_markers: usize,
    },
}

fn main() {
    let arguments = Arguments::parse();

    match &arguments.command {
        Commands::KmerRatio {
            input_bam,
            threads,
            prior_mean,
            prior_weight,
        } => {
            if let Err(error) = kmer_ratio::run(
                input_bam,
                *threads,
                *prior_mean,
                *prior_weight,
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
            ratio_threshold,
            min_markers,
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
                *ratio_threshold,
                *min_markers,
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
            ratio_threshold,
            min_markers,
        } => {
            if let Err(error) = single::run(
                input_bam,
                output_bam,
                ref_tabix,
                *kmer_size,
                *threads,
                *ratio_threshold,
                *min_markers,
            ) {
                eprintln!("{}", error);
                process::exit(1);
            }
        },

    }    
}

