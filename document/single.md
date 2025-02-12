# How to use bam_refiner for single reference genome

## Installation
Create singularity image

```
singularity pull bam_refiner_${version}.sif docker://yosakam2/bam_refiner:${version}
```

## Input
- BAM/FASTQ file
- A list of target regions whose format is chr:start-end
- Reference genome (FASTA)
- Output directory 
- Sample name
- Thread number (optional, default 8)

## Output
- refined BAM file
- bam_refiner_result.tsv.gz
- bam_refiner.log.gz
- kmerposition.bed.gz

## How to use

If you have a BAM file, then

```
singularity exec bam_refiner_${version}.sif \
run_single.sh \
    -b BAM_FILE \
    -d \
    -l REGIONS_FILE \
    -o OUTPUT_DIR \
    -p \
    -r REFERENCE_FASTA \
    -s SAMPLE_NAME \
    -t THREAD
```

Otherwise,
```
singularity exec bam_refiner_${version}.sif \
run_single.sh \
    -f FASTQ_FILE \
    -d \
    -l REGIONS_FILE \
    -m MINIMAP_OPTION \ e.g. "-ax map-ont"
    -o OUTPUT_DIR \
    -p \
    -r REFERENCE_FASTA \
    -s SAMPLE_NAME \
    -t THREAD
```
