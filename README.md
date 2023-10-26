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

## Preparation
### Step 1: Diploid genome assembly by Hifiasm (Verkko)

### Step2: Extract haplotype-specific unique k-mer
You can use [meryl](https://github.com/marbl/meryl.git) to count kmers.\
You can also use [kmer_locate](https://github.com/yos-sk/kmer_locate.git) to search kmer positions. 

```
for hap in hap1 hap2
do
    meryl count k=21 threads=16 ${hap}.contig.fa output ${hap}.meryl
done

meryl difference hap1.meryl hap2.meryl output hap1.uniq.meryl
meryl difference hap2.meryl hap1.meryl output hap2.uniq.meryl

meryl print threads=16 hap1.uniq.meryl > hap1.uniq.tsv
meryl print threads=16 hap2.uniq.meryl > hap2.uniq.tsv

for hap in hap1 hap2
do
    awk '{if ($2 <= 10) print}' ${hap}.uniq.tsv > ${hap}.cnt10.uniq.tsv
    gzip ${hap}.cnt10.uniq.tsv
    gzip ${hap}.uniq.tsv
done

for hap in hap1 hap2
do
    python ${path-to-kmer_locate}/script/kmercounts2fasta.py ${hap}.cnt10.uniq.tsv.gz > unique_kmerCounts_cnt10_${hap}.fa ${hap}
    kmer_locate --kmer-path unique_kmerCounts_cnt10_${hap}.fa --input-file ${hap}.contig.fa --kmer-size 21 | sort -k 1,1 -k 2,2n > ${hap}_cnt10_kmerposition.bed
    bgzip -f ${hap}_cnt10_kmerposition.bed
    tabix -p bed ${hap}_cnt10_kmerposition.bed.gz
done
```
### Step3: Align sequencing reads to the diploid genome assembly constructed in Step1
You can use [minimap2](https://github.com/lh3/minimap2.git) and [samtools](http://www.htslib.org) for alignment.\
You should sort the bam file by read name for [bam_refiner](https://github.com/yos-sk/bam_refiner.git).


## Usage
Please split the BAM file using [split_bam](https://github.com/yos-sk/split_bam.git) if necessary.\

```
./target/release/bam_refiner \
    --input-bam ${INPUT_BAM} \
    --output-bam ${OUTPUT_BAM} \
    --hap1-tabix hap1_cnt10_kmerposition.bed.gz \
    --hap2-tabix hap2_cnt10_kmerposition.bed.gz \
    --kmer-size 21 \
    1>output.tsv 2>log
```

If you want to use bam_refiner for the alignment data to the exisitng reference genome (e.g. GRCh38 or CHM13), plase try single_mode branch and see [document](https://github.com/yos-sk/bam_refiner/blob/master/document/single_mode.md).
