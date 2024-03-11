use std::collections::HashMap;
use std::error::Error;
use std::io::BufRead;

use bam_refiner::open_file;

pub fn run(kmer_file: &str, input_fasta: &str, kmer_size: usize) -> Result<(), Box<dyn Error>> {
    let kmer_db = make_kmer_db(kmer_file).expect("Failed to make kmer db");
    let fasta_reader = open_file(input_fasta).expect(&format!("Failed to open file {}", input_fasta));

    let mut sequence = String::new();
    let mut header = String::new();
    for line in fasta_reader.lines() {
        let line = line.unwrap();
        if line.starts_with(">") {
            if !sequence.is_empty() {
                search_kmers(&header, &sequence.to_uppercase(), &kmer_db, kmer_size);
                sequence.clear();
            }
            let split_cname: Vec<&str> = line.split(' ').collect();
            header = split_cname[0].to_string();
        } else {
            sequence.push_str(&line);
        }
    }

    if !sequence.is_empty() {
        search_kmers(&header, &sequence.to_uppercase(), &kmer_db, kmer_size);
    }
    Ok(())
}

fn search_kmers(header: &str, sequence: &str, kmer_db: &HashMap<String, String>, k: usize) {
    for i in 0..(sequence.len() - k + 1) {
        let slice = &sequence[i..(i + k)];
        if let Some(value) = kmer_db.get(slice) {
            println!(
                "{}\t{}\t{}\t{}\t+\t{}",
                &header[1..],
                i,
                i + k,
                value,
                slice
            );
        }
    }

    let rev_comp_sequence = reverse_complement(sequence);
    for i in 0..(rev_comp_sequence.len() - k + 1) {
        let slice = &rev_comp_sequence[i..(i + k)];
        if let Some(value) = kmer_db.get(slice) {
            println!(
                "{}\t{}\t{}\t{}\t-\t{}",
                &header[1..],
                rev_comp_sequence.len() - i - k,
                rev_comp_sequence.len() - i,
                value,
                slice
            );
        }
    }
}

fn reverse_complement(sequence: &str) -> String {
    // complement
    let complement = sequence
        .chars()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => c,
        })
        .collect::<String>();

    // reverse
    let rev_comp = complement.chars().rev().collect::<String>();

    rev_comp
}

fn make_kmer_db(kmer_file: &str) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let reader = open_file(kmer_file).expect("Failed to open file");

    let mut kmer_set: HashMap<String, String> = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let split_line: Vec<&str> = line.split('\t').collect();
        let kmer = split_line[0].to_string();
        let kmer_id = "kmer".to_string() + &i.to_string();
        kmer_set.insert(kmer, kmer_id);
    }

    Ok(kmer_set)
}