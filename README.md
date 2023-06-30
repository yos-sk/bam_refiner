# Marker_filter
Filter alignments by unique kmers

## Install 
```
git clone https://github.com/yos-sk/marker_filter.git
cd marker_filter
cargo build --release
```

## Usage
```
./target/release/marker_filter -input-bam ${INPUT_BAM} --output-bam ${OUTPUT_BAM} --hap1-tabix ${hap1_tabix} --hap2-tabix ${hap2_tabix} --kmer-size ${kmer_size} 1>output.tsv 2>log
```
