# How to use bam_refiner single reference mode
Alignments of reads to reference genomes should be done before running this tool.

## Preparation
### Step 1: Extract sequences of target regions

```
# Please specify the ${reference}.
samtools faidx ${reference}
samtools faidx ${reference} \
    -o target.fa \
    -r target.txt
```
You should prepare the fasta file of target regions whose format is chr:from-to using [samtools](http://www.htslib.org).

### Step 2: Extract unique k-mer
You can use [meryl](https://github.com/marbl/meryl.git) to count kmers. \
You can also use [kmer_locate](https://github.com/yos-sk/kmer_locate.git) to search kmer positions. 

```
meryl count \
    k=21 \
    threads=16 \
    target.fa \
    output target.meryl
meryl print \
    threads=16 \
    target.meryl > target.tsv
awk '{if ($2 <= 10) print}' target.tsv > target.cnt10.tsv
gzip target.cnt10.tsv

python3 ${path-to-kmer_locate}/script/kmercounts2fasta.py \
    target.cnt10.tsv.gz target > kmerCounts_cnt10.fa
kmer_locate \
    --kmer-path kmerCounts_cnt10.fa \
    --input-file ${reference} \
    --kmer-size 21 | sort -k 1,1 -k 2,2n > kmerposition.bed
bgzip -f kmerposition.bed
tabix -p bed kmerposition.bed.gz
```
You should use the same ${reference} as Step 1.

### Step 3: Align sequencing reads to thd diploid genome assembly constructed in Step 1
You can use [minimap2](https://github.com/lh3/minimap2.git) and [samtools](http://www.htslib.org) for alignment.\
You should sort the bam file by read name for [bam_refiner](https://github.com/yos-sk/bam_refiner.git).

```
# For ONT data
minimap2 -t 16 -ax map-ont ${reference} input.fastq | samtools view --Shb > output.unsorted
# For HiFi data
minimap2 -t 16 -ax map-hifi ${reference} input.fastq | samtools view --Shb > output.unsorted

samtools sort -@ 16 -m 2G -n output.unsorted -o output.bam
samtools index output.bam
```

# Usage: bam_refiner single reference mode
## Split bam (Optional)
Split a bam file using [split_bam](https://github.com/yos-sk/split_bam.git) if you necessary.
```
# Split bam
SIZE=`split_bam size --input-file ${INPUT_BAM}`
split_bam split \
    --input-file ${INPUT_BAM} \
    --output-dir ${OUTPUT_DIR} \
    --input-size ${SIZE} \
    --num-split 8
```
## Refine bam
You can perform array job if you split a bam file.
```
bam_refiner single \
    --input-bam input.bam \
    --output-bam output.refined.bam \
    --ref-tabix kmerposition.bed.gz \
    --kmer-size 21 \
    1> bam_refiner.tsv 2> bam_refiner.err
```
## Merge bam files (Optional)
You should merge bam files when you split a bam file.
```
samtools merge \
    -@ 8 \
    -o Merged.bam \
    *.refined.bam

samtools sort -@ 8 -o Sorted.bam Merged.bam
samtools index Sorted.bam
rm Merged.bam
```