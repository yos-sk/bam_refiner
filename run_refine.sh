#!/bin/bash

set -xv
set -o errexit
set -o nounset
set -o pipefail

while getopts "b:df:h:i:M:m:o:pR:s:t:u:" opt; do
  case $opt in
    b) INPUT_BAM=$OPTARG ;;
    d) DEBUG="true" ;;
    f) INPUT_FASTQ=$OPTARG ;;
    h) HAP1_CONTIG=$OPTARG ;;
    i) HAP2_CONTIG=$OPTARG ;;
    m) MINIMAP_OPTION=$OPTARG ;;
    o) OUTPUT_DIR=$OPTARG ;;
    p) OPTION_SPLIT="true" ;;
    M) MIN_MARKERS=$OPTARG ;;
    R) RATIO_THRESHOLD=$OPTARG ;;
    s) SAMPLE=$OPTARG ;;
    t) THREAD=$OPTARG ;;
    u) DATA=$OPTARG ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done

if [ -z "${INPUT_FASTQ:-}" ] && [ -z "${INPUT_BAM:-}" ]; then
    echo "FASTQ/BAM file is not given. Please set -f {FASTQ_FILE} or -b {BAM_FILE}"; exit 1
fi

if [ -z "${DEBUG:-}" ]; then
    echo "Do not remove a workspace"
    DEBUG="false"
fi

if [ -z "${OUTPUT_DIR:-}" ]; then
    echo "Output directory is not given. Please set -o {OUTPUT_DIR}"; exit 1
fi

if [ -z "${MINIMAP_OPTION:-}" ]; then
    MINIMAP_OPTION=""
    echo "Not given any options for minimap2"
fi

if [ -z "${OPTION_SPLIT:-}" ]; then
    echo "Split option is not given. Bam_refiner will be performed without splitting a BAM file"
    OPTION_SPLIT="false"
fi

if [ -z "${RATIO_THRESHOLD:-}" ]; then
    echo "Ratio threshold is not given. Ratio threshold is set to default (0.5)"
    RATIO_THRESHOLD=0.5
fi

if [ -z "${MIN_MARKERS:-}" ]; then
    echo "Minimum marker count is not given. It is set to default (1 = disabled)"
    MIN_MARKERS=1
fi

if [ -z "${SAMPLE:-}" ]; then
    echo "Sample name is not given. Please set -s {SAMPLE}"; exit 1
fi

WORK_DIR=${OUTPUT_DIR}/workspace
mkdir -p ${WORK_DIR}

# Step1: Mapping
mkdir -p ${WORK_DIR} 
OUTPUT_BAM_PREFIX=${WORK_DIR}/${SAMPLE}
cat ${HAP1_CONTIG} ${HAP2_CONTIG} > ${WORK_DIR}/reference.fa

if [ -z ${INPUT_BAM:-} ] && [ ! -z ${INPUT_FASTQ:-} ]; then
    if [ ${DATA} == "hifi" ]; then
        minimap2 -t ${THREAD} -ax asm5 ${MINIMAP_OPTION} ${WORK_DIR}/reference.fa ${INPUT_FASTQ} | samtools view -Shb > ${OUTPUT_BAM_PREFIX}.unsorted
    else
        minimap2 -t ${THREAD} -ax asm10 ${MINIMAP_OPTION} ${WORK_DIR}/reference.fa ${INPUT_FASTQ} | samtools view -Shb > ${OUTPUT_BAM_PREFIX}.unsorted
    fi
elif [ ! -z ${INPUT_BAM:-} ]; then
    if [ ${DATA} == "hifi" ]; then
        samtools fastq -@ ${THREAD} -TMM,ML ${INPUT_BAM} | minimap2 -t ${THREAD} -ax asm5 -y ${MINIMAP_OPTION} ${WORK_DIR}/reference.fa - | samtools view -@ ${THREAD} -Shb - > ${OUTPUT_BAM_PREFIX}.unsorted
    else
        samtools fastq -@ ${THREAD} -TMM,ML ${INPUT_BAM} | minimap2 -t ${THREAD} -ax asm10 -y ${MINIMAP_OPTION} ${WORK_DIR}/reference.fa - | samtools view -@ ${THREAD} -Shb - > ${OUTPUT_BAM_PREFIX}.unsorted
    fi
fi
 
samtools sort -@ ${THREAD} -m 2G -n ${OUTPUT_BAM_PREFIX}.unsorted -o ${OUTPUT_BAM_PREFIX}.bam
rm ${OUTPUT_BAM_PREFIX}.unsorted

# Step2: Extract haplotype-specific unique k-mer
mkdir -p ${WORK_DIR}/meryl
meryl count k=21 threads=${THREAD} ${HAP1_CONTIG} output ${WORK_DIR}/meryl/hap1.meryl
meryl count k=21 threads=${THREAD} ${HAP2_CONTIG} output ${WORK_DIR}/meryl/hap2.meryl


meryl difference ${WORK_DIR}/meryl/hap1.meryl ${WORK_DIR}/meryl/hap2.meryl output ${WORK_DIR}/meryl/hap1.uniq.meryl
meryl difference ${WORK_DIR}/meryl/hap2.meryl ${WORK_DIR}/meryl/hap1.meryl output ${WORK_DIR}/meryl/hap2.uniq.meryl

meryl print threads=${THREAD} ${WORK_DIR}/meryl/hap1.uniq.meryl > ${WORK_DIR}/meryl/hap1.uniq.tsv
meryl print threads=${THREAD} ${WORK_DIR}/meryl/hap2.uniq.meryl > ${WORK_DIR}/meryl/hap2.uniq.tsv

for hap in hap1 hap2
do
    awk '{if ($2 == 1) print}' ${WORK_DIR}/meryl/${hap}.uniq.tsv > ${WORK_DIR}/meryl/${hap}.cnt.uniq.tsv
    gzip -f ${WORK_DIR}/meryl/${hap}.cnt.uniq.tsv
    gzip -f ${WORK_DIR}/meryl/${hap}.uniq.tsv
done


mkdir -p ${OUTPUT_DIR}
bam_refiner locate-kmers \
    -i ${WORK_DIR}/meryl/hap1.cnt.uniq.tsv.gz \
    -f ${HAP1_CONTIG} \
    -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed.gz

bam_refiner locate-kmers \
    -i ${WORK_DIR}/meryl/hap2.cnt.uniq.tsv.gz \
    -f ${HAP2_CONTIG} \
    -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed
bgzip -f ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed
tabix -p bed ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed.gz

# Step 3; List up contig names
grep ">" ${HAP1_CONTIG} | sed s/\>// > ${OUTPUT_DIR}/hap1_list.txt
gzip -f ${OUTPUT_DIR}/hap1_list.txt
grep ">" ${HAP2_CONTIG} | sed s/\>// > ${OUTPUT_DIR}/hap2_list.txt
gzip -f ${OUTPUT_DIR}/hap2_list.txt

# Step 4: Refine BAM file (multi-threaded; no physical splitting required)
bam_refiner refine \
    --input-bam ${OUTPUT_BAM_PREFIX}.bam \
    --output-bam ${OUTPUT_DIR}/${SAMPLE}_bam_refined.bam \
    --hap1-tabix ${OUTPUT_DIR}/hap1_cnt_kmerposition.bed.gz \
    --hap2-tabix ${OUTPUT_DIR}/hap2_cnt_kmerposition.bed.gz \
    --hap1-list ${OUTPUT_DIR}/hap1_list.txt.gz \
    --hap2-list ${OUTPUT_DIR}/hap2_list.txt.gz \
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

# Reads spanning no haplotype-specific locus fall back to this value instead of
# claiming a perfect 1.0. ONT is noisier, so its no-evidence reads get less credit.
if [ ${DATA} == "hifi" ]; then
    PRIOR_MEAN=0.8
else
    PRIOR_MEAN=0.6
fi

bam_refiner kmer-ratio \
    ${OUTPUT_DIR}/${SAMPLE}_bam_refined.sorted.bam \
    --threads ${THREAD} \
    --prior-mean ${PRIOR_MEAN} \
> ${OUTPUT_DIR}/${SAMPLE}_kmer_ratio.txt

if [ ${DEBUG} = "false" ]; then
    rm -rf ${WORK_DIR}
fi
