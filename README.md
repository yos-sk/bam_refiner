# Bam_refiner
Refine alignments by unique kmers

## Dependencies
- clap >=4.3.9
- rust_htslib 0.43.0

I have not tried other versions of rust_htslib.

## Install 
```
git clone https://github.com/yos-sk/bam_refiner.git
cd bam_refiner
cargo build --release
```

## Usage
```
./target/release/bam_refiner -input-bam ${INPUT_BAM} --output-bam ${OUTPUT_BAM} --hap1-tabix ${hap1_tabix} --hap2-tabix ${hap2_tabix} --kmer-size ${kmer_size} 1>output.tsv 2>log
```