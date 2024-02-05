use clap::Parser;
//use std::process;

mod collect_alignments;
mod filter;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.3.3", about = "Refine alignments by unique kmers.", long_about = None)]
struct Arguments {
    #[arg(short = 'i', long)]
    input_bam: String,

    #[arg(short = 't', long)]
    threads: usize,

    /*
    #[arg(short = 'o', long)]
    output_bam: String,
    */
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
}

fn main() {
    let arguments = Arguments::parse();
    let (alignments, sequences) = collect_alignments::run(
        &arguments.input_bam,
        arguments.threads,
    ).expect("Error");

    let _ = filter::run(alignments, sequences, &arguments.hap1_tabix, &arguments.hap2_tabix, &arguments.hap1_list, &arguments.hap2_list, arguments.kmer_size, arguments.threads).expect("Error");
}



