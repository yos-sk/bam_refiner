#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

while getopts "b:df:l:M:m:o:pR:r:s:t:" opt; do
  case $opt in
    b) INPUT_BAM=$OPTARG ;;
    d) DEBUG="true" ;;
    f) INPUT_FASTQ=$OPTARG ;;
    l) REGIONS=$OPTARG ;;
    m) MINIMAP_OPTION=$OPTARG ;;
    o) OUTPUT_DIR=$OPTARG ;;
    p) OPTION_SPLIT="true" ;;
    M) MIN_MARKERS=$OPTARG ;;
    R) RATIO_THRESHOLD=$OPTARG ;;
    r) REFERENCE=$OPTARG ;;
    s) SAMPLE=$OPTARG ;;
    t) THREAD=$OPTARG ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done

if [ -z "${INPUT_FASTQ:-}" ] && [ -z "${INPUT_BAM:-}" ]; then
    echo "FASTQ/BAM file is not given. Please set -f {FASTQ_FILE} with mapping or -b {BAM_FILE} without mapping"; exit 1
fi

if [ -z "${DEBUG:-}" ]; then
    echo "Do not remove a workspace"
    DEBUG="false"
fi

if [ -z "${REGIONS:-}" ]; then
    echo "Target region list is not given. Pleaset -l {REGION_LIST}"; exit 1
fi

if [ -z "${MINIMAP2_OPTION:-}" ] && [ -z "${INPUT_BAM:-}" ]; then
    echo "Minimap2 option is not given. Please set -m {MINIMAP2_OPTION}"; exit 1
fi

if [ -z "${OUTPUT_DIR:-}" ]; then
    echo "Output directory is not given. Please set -o {OUTPUT_DIR}"; exit 1
fi

if [ -z "${OPTION_SPLIT:-}" ]; then
    echo "Split option is not given. Bam_refiner will be performed without splitting a BAM file"
    OPTION_SPLIT="false"
fi

if [ -z "${RATIO_THRESHOLD:-}" ]; then
    echo "Ratio threshold is not given. Ratio threshold is set to default (0.8)"
    RATIO_THRESHOLD=0.8
fi

if [ -z "${MIN_MARKERS:-}" ]; then
    echo "Minimum marker count is not given. It is set to default (3)"
    MIN_MARKERS=3
fi

if [ -z "${REFERENCE:-}" ]; then
    echo "Reference genome is not given. Please set -r {REFERENCE}"; exit 1
fi

if [ -z "${SAMPLE:-}" ]; then
    echo "Sample name is not given. Please set -s {SAMPLE}"; exit 1
fi

if [ -z "${THREAD:-}" ]; then
    echo "Thread number is not given. Thread number is set to default (8)"
    THREAD=8
fi

WORK_DIR=${OUTPUT_DIR}/workspace
mkdir -p ${WORK_DIR} 

# Step1: Mapping
OUTPUT_BAM_PREFIX=${WORK_DIR}/${SAMPLE}
if [ -z ${INPUT_BAM:-} ] && [ ! -z ${INPUT_FASTQ:-} ]; then
    minimap2 -t ${THREAD} ${MINIMAP2_OPTION} ${REFERENCE} ${INPUT_FASTQ} | samtools view -@ ${THREAD}-Shb - > ${OUTPUT_BAM_PREFIX}.unsorted
    samtools sort -@ ${THREAD} -m 2G -n ${OUTPUT_BAM_PREFIX}.unsorted -o ${OUTPUT_BAM_PREFIX}.bam
    rm ${OUTPUT_BAM_PREFIX}.unsorted
    BAM=${OUTPUT_BAM_PREFIX}.bam
elif [ ! -z ${INPUT_BAM:-} ]; then
    samtools sort -@ ${THREAD} -m 2G -n ${INPUT_BAM} -o ${OUTPUT_BAM_PREFIX}.bam
    BAM=${OUTPUT_BAM_PREFIX}.bam 
fi

# Step2: Extract reference sequences of target regions
if [ ! -f ${REFERENCE}.fai ]; then
    samtools faidx ${REFERENCE}
fi
samtools faidx ${REFERENCE} -o ${WORK_DIR}/target_regions.fa -r ${REGIONS}

# Step3: Extract region-specific k-mer
mkdir -p ${WORK_DIR}/meryl
if [ -d ${WORK_DIR}/meryl/target.meryl ]; then
    rm -r ${WORK_DIR}/meryl/target.meryl
fi
meryl count k=21 threads=${THREAD} ${WORK_DIR}/target_regions.fa output ${WORK_DIR}/meryl/target.meryl
meryl print threads=${THREAD} ${WORK_DIR}/meryl/target.meryl > ${WORK_DIR}/meryl/target_kmers.tsv

awk '{if ($2 == 1) print}' ${WORK_DIR}/meryl/target_kmers.tsv > ${WORK_DIR}/meryl/target_kmers_uniq.tsv
gzip -f ${WORK_DIR}/meryl/target_kmers_uniq.tsv

bam_refiner locate-kmers \
    -i ${WORK_DIR}/meryl/target_kmers_uniq.tsv.gz \
    -f ${REFERENCE} \
    -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/kmerposition.bed
bgzip -f ${OUTPUT_DIR}/kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/kmerposition.bed.gz

# Step 4: Refine BAM file (multi-threaded; no physical splitting required)
bam_refiner single \
    --input-bam ${BAM} \
    --output-bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam \
    --ref-tabix ${OUTPUT_DIR}/kmerposition.bed.gz \
    --kmer-size 21 \
    --threads ${THREAD} \
    --ratio-threshold ${RATIO_THRESHOLD} \
    --min-markers ${MIN_MARKERS} \
    1>${OUTPUT_DIR}/bam_refiner_result.tsv 2>${OUTPUT_DIR}/bam_refiner.log

samtools sort \
    -@ ${THREAD} \
    -o ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam \
    ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam 
samtools index -@ ${THREAD} ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam 
rm ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam

gzip -f ${OUTPUT_DIR}/bam_refiner_result.tsv
gzip -f ${OUTPUT_DIR}/bam_refiner.log

if [ ${DEBUG} = "false" ]; then
    rm -rf ${WORK_DIR}
fi
