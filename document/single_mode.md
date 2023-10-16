# How to use bam_refiner single mode
Alignments of reads to reference genomes should be done before running this tool.

# Step 1: Extract sequences of target regions

```
# Please specify the ${reference}.
samtools faidx ${reference}
samtools faidx ${reference} \
    -o target.fa \
    -r target.txt
```
You should prepare the target.txt whose format is chr:from-to.

## Step 2: Extract unique k-mer

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

## Step 2: Refine alignments with k-mers
### Split bam
Split a bam file if you necessary.
```
# Split bam
SIZE=`singularity exec ~/bin/bam_refiner/bam_refiner_single.sif 
split_bam size --input-file ${INPUT_BAM}`
split_bam split \
    --input-file ${INPUT_BAM} \
    --output-dir ${OUTPUT_DIR} \
    --input-size ${SIZE} \
    --num-split 8
```
### Refine bam
Please use "-t 1:8" qsub option and add "TASK_ID=$(( ${SGE_TASK_ID} - 1 ))" to qsub script if you split a bam file.
```
bam_refiner single \
    --input-bam ${OUTPUT_DIR}/${TASK_ID}.bam \
    --output-bam ${OUTPUT_DIR}/${TASK_ID}.refined.bam \
    --ref-tabix ${OUTPUT_DIR}/kmerposition.bed.gz \
    --kmer-size 21 \
    1> ${OUTPUT_DIR}/${TASK_ID}.marker_filter.tsv 2> ${OUTPUT_DIR}/${TASK_ID}.err.tsv
```
### Merge bam files
```
samtools merge \
    -@ 8 \
    -o Merged.bam \
    ${OUTPUT_DIR}/0.refined.bam ${OUTPUT_DIR}/1.refined.bam \
    ${OUTPUT_DIR}/2.refined.bam ${OUTPUT_DIR}/3.refined.bam \
    ${OUTPUT_DIR}/4.refined.bam ${OUTPUT_DIR}/5.refined.bam \
    ${OUTPUT_DIR}/6.refined.bam ${OUTPUT_DIR}/7.refined.bam

samtools sort -@ 8 -o Sorted.bam Merged.bam
samtools index Sorted.bam
rm Merged.bam
```