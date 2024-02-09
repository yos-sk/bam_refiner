use clap::Parser;
//use std::process;

mod collect_alignments;
mod kmer_count;
mod refine;
mod write_bam;

#[derive(Parser)]
#[command(author = "Yoshitaka Sakamoto", version = "0.3.3", about = "Refine alignments by unique kmers.", long_about = None)]
struct Arguments {
    #[arg(short = 'i', long)]
    input_bam: String,

    #[arg(short = 't', long)]
    threads: usize,

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
}

fn main() {
    let arguments = Arguments::parse();
    let (alignments, sequences, qualities, indices) =
        collect_alignments::run(&arguments.input_bam, arguments.threads)
            .expect("Error: collect_alignments step");

    let kmer_counted_data = kmer_count::run(
        alignments,
        sequences.clone(),
        &arguments.hap1_tabix,
        &arguments.hap2_tabix,
        &arguments.hap1_list,
        &arguments.hap2_list,
        arguments.kmer_size,
        arguments.threads,
    )
    .expect("Error: kmer_count step");
    let refined_alignments = refine::run(kmer_counted_data, indices).expect("Error: refine step");
    let _ = write_bam::run(
        &arguments.input_bam,
        &arguments.output_bam,
        &refined_alignments,
        &sequences,
        &qualities,
        &arguments.hap1_list,
        &arguments.hap2_list,
    );
}
