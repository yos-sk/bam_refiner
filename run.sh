#!/bin/bash

SAMPLE=$1
hap1_contig=$2
hap2_contig=$3
FASTQ=$4
OPTION_SPLIT=$5
WORK_DIR=$6
OUTPUT_DIR=$7

set -xv
set -o errexit
set -o nounset
set -o pipefail


# Step1: Mapping
mkdir -p ${WORK_DIR} 
OUTPUT_BAM_PREFIX=${WORK_DIR}/${SAMPLE}
cat ${hap1_contig} ${hap2_contig} > ${WORK_DIR}/reference.fa
minimap2 -t 16 -ax asm10 ${WORK_DIR}/reference.fa ${FASTQ} | samtools view -Shb > ${OUTPUT_BAM_PREFIX}.unsorted
samtools sort -@ 16 -m 2G -n ${OUTPUT_BAM_PREFIX}.unsorted -o ${OUTPUT_BAM_PREFIX}.bam
samtools index ${OUTPUT_BAM_PREFIX}.bam
rm ${OUTPUT_BAM_PREFIX}.unsorted

# Step2: Extract haplotype-specific unique k-mer
mkdir -p ${WORK_DIR}/meryl
meryl count k=21 threads=16 ${hap1_contig} output ${WORK_DIR}/meryl/hap1.meryl
meryl count k=21 threads=16 ${hap2_contig} output ${WORK_DIR}/meryl/hap2.meryl


meryl difference ${WORK_DIR}/meryl/hap1.meryl ${WORK_DIR}/meryl/hap2.meryl output ${WORK_DIR}/meryl/hap1.uniq.meryl
meryl difference ${WORK_DIR}/meryl/hap2.meryl ${WORK_DIR}/meryl/hap1.meryl output ${WORK_DIR}/meryl/hap2.uniq.meryl

meryl print threads=16 ${WORK_DIR}/meryl/hap1.uniq.meryl > ${WORK_DIR}/meryl/hap1.uniq.tsv
meryl print threads=16 ${WORK_DIR}/meryl/hap2.uniq.meryl > ${WORK_DIR}/meryl/hap2.uniq.tsv

for hap in hap1 hap2
do
    awk '{if ($2 <= 10) print}' ${WORK_DIR}/meryl/${hap}.uniq.tsv > ${WORK_DIR}/meryl/${hap}.cnt10.uniq.tsv
    gzip ${WORK_DIR}/meryl/${hap}.cnt10.uniq.tsv
    gzip ${WORK_DIR}/meryl/${hap}.uniq.tsv
done


mkdir -p ${OUTPUT_DIR}
python3 /tools/kmer_locate/script/kmercounts2fasta.py ${WORK_DIR}/meryl/hap1.cnt10.uniq.tsv.gz > ${WORK_DIR}/meryl/unique_kmerCounts_cnt10_hap1.fa hap1
kmer_locate --kmer-path ${WORK_DIR}/meryl/unique_kmerCounts_cnt10_${hap}.fa --input-file ${hap1_contig} --kmer-size 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap1_cnt10_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap1_cnt10_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap1_cnt10_kmerposition.bed.gz

python /tools/kmer_locate/script/kmercounts2fasta.py ${WORK_DIR}/meryl/hap2.cnt10.uniq.tsv.gz > ${WORK_DIR}/meryl/unique_kmerCounts_cnt10_hap2.fa hap2
kmer_locate --kmer-path ${WORK_DIR}/meryl/unique_kmerCounts_cnt10_hap2.fa --input-file ${hap2_contig} --kmer-size 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap2_cnt10_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap2_cnt10_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap2_cnt10_kmerposition.bed.gz

# Step 3; List up contig names
grep ">" ${hap1_contig} | sed s/\>// > ${OUTPUT_DIR}/hap1_list.txt
gzip ${OUTPUT_DIR}/hap1_list.txt
grep ">" ${hap2_contig} | sed s/\>// > ${OUTPUT_DIR}/hap2_list.txt
gzip ${OUTPUT_DIR}/hap2_list.txt

# Step 4: Refine BAM file
mkdir -p ${OUTPUT_DIR}
if [ $OPTION_SPLIT = "true" ]
then
    mkdir -p ${WORK_DIR}/split
    SIZE=`split_bam size --input-file ${OUTPUT_BAM_PREFIX}.bam`
    split_bam split \
        --input-file ${OUTPUT_BAM_PREFIX}.bam \
        --output-dir ${WORK_DIR}/split \
        --input-size ${SIZE} \
        --num-split 8
    
    for i in {0..7}; do
        bam_refiner \
            --input-bam ${WORK_DIR}/split/${i}.bam \
            --output-bam ${WORK_DIR}/split/${i}.refined.bam \
            --hap1-tabix ${OUTPUT_DIR}/hap1_cnt10_kmerposition.bed.gz \
            --hap2-tabix ${OUTPUT_DIR}/hap2_cnt10_kmerposition.bed.gz \
            --hap1-list ${OUTPUT_DIR}/hap1_list.txt.gz \
            --hap2-list ${OUTPUT_DIR}/hap2_list.txt.gz \
            --kmer-size 21 \
            1>${WORK_DIR}/split/${i}.bam_refiner.tsv 2>${WORK_DIR}/split/${i}.bam_refiner.log &
    done
    wait
    
    cat ${WORK_DIR}/split/${i}.bam_refiner.tsv > ${OUTPUT_DIR}/bam_refiner_result.tsv
    cat ${WORK_DIR}/split/${i}.bam_refiner.log > ${OUTPUT_DIR}/bam_refiner.log
    samtools merge \
        -@ 8 \
        -o ${WORK_DIR}/${SAMPLE}_bam_refined.bam \
        ${WORK_DIR}/split/*.refined.bam

    samtools sort -@ 8 -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam  ${WORK_DIR}/${SAMPLE}_bam_refined.bam 
    samtools index ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam 
else
    bam_refiner \
        --input-bam ${OUTPUT_BAM_PREFIX}.bam \
        --output-bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam \
        --hap1-tabix ${OUTPUT_DIR}/hap1_cnt10_kmerposition.bed.gz \
        --hap2-tabix ${OUTPUT_DIR}/hap2_cnt10_kmerposition.bed.gz \
        --hap1-list ${OUTPUT_DIR}/hap1_list.txt.gz \
        --hap2-list ${OUTPUT_DIR}/hap2_list.txt.gz \
        --kmer-size 21 \
        1>${OUTPUT_DIR}/bam_refiner_result.tsv 2>${OUTPUT_DIR}/log/bam_refiner.log
    
    samtools sort -@ 8 -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam 
    samtools index ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam 
    rm ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam 
fi
gzip ${OUTPUT_DIR}/bam_refiner_result.tsv
gzip ${OUTPUT_DIR}/log/bam_refiner.log

rm -rf ${WORK_DIR}