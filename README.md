# bam_refiner

![Refine strategy overview](images/Refine_strategy_v2.png)

Refine long-read alignments on a diploid genome assembly using haplotype-specific k-mers.

When long reads are aligned to a concatenated diploid assembly (hap1 + hap2), reads originating from one haplotype can be placed on the other because the two haplotypes are nearly identical and only a small fraction of positions are truly distinguishing. `bam_refiner` re-evaluates each alignment by counting *haplotype-unique* k-mers — k-mers that occur in only one of the two haplotypes — carried by the read, and reassigns the read to the haplotype it most strongly supports. Alignments that cannot be confidently attributed to either haplotype are filtered out.

The typical workflow is: (1) build haplotype-unique k-mer sets with [meryl](https://github.com/marbl/meryl.git), (2) locate those k-mers on each haplotype with `bam_refiner locate-kmers`, (3) align reads to the concatenated assembly, and (4) run `bam_refiner refine` to produce a haplotype-aware BAM. The `bam_refiner kmer-ratio` subcommand additionally reports a per-read hap1/hap2 k-mer ratio, which downstream tools (e.g. [PRCGAP](https://github.com/yos-sk/PRCGAP)) use for haplotype-resolved somatic variant calling.

Both PacBio HiFi and Oxford Nanopore reads are supported.


## Install
```
git clone https://github.com/yos-sk/bam_refiner.git
cd bam_refiner
cargo build --release
```

## Preparation 
### Diploid genome assembly
You can use [hifiasm](https://github.com/chhylp123/hifiasm.git) or [verkko](https://github.com/marbl/verkko.git) to perform diploid genome assembly.

## Usage
### 1. Container
A Docker image of bam_refiner is published on Docker Hub. The image bundles `bam_refiner`, `split_bam`, `minimap2`, and `samtools`, and ships the `run_refine.sh` end-to-end wrapper.

```
docker pull yosakam2/bam_refiner:${VERSION}
```

Run with Docker (mount the working directory so inputs/outputs are visible inside the container):
```
docker run --rm \
    -v ${PWD}:${PWD} -w ${PWD} \
    yosakam2/bam_refiner:${VERSION} \
    /bin/bash /tools/bam_refiner/run_refine.sh \
        -d \ # For debug mode to keep intermediate files
        -f ${FASTQ} \
        -h ${HAP1_ASSEMBLY} \ # fasta file of haplotype1 contigs
        -i ${HAP2_ASSEMBLY} \ # fasta file of haplotype2 contigs
        -o ${OUTPUT_DIR} \
        -p \ # For parallel processing of bam_refiner
        -s ${SAMPLE_NAME} \ # Sample name
        -t ${THREADS} \ # Number of threads
        -u ${DATA_TYPE} # hifi or ont
```

On HPC environments without a Docker daemon, pull the same image as a Singularity image and run it with `singularity exec`:
```
singularity pull bam_refiner_${VERSION}.sif docker://yosakam2/bam_refiner:${VERSION}

singularity exec bam_refiner_${VERSION}.sif \
    /bin/bash /tools/bam_refiner/run_refine.sh \
        -d \
        -f ${FASTQ} \
        -h ${HAP1_ASSEMBLY} \
        -i ${HAP2_ASSEMBLY} \
        -o ${OUTPUT_DIR} \
        -p \
        -s ${SAMPLE_NAME} \
        -t ${THREADS} \
        -u ${DATA_TYPE}
```

### 2. Step by step
#### Step 1: Extract haplotype-specific unique k-mer
You can use [meryl](https://github.com/marbl/meryl.git) to count kmers. 

```
for hap in hap1 hap2
do
    meryl count k=21 threads=16 ${hap}_contig.fa output ${hap}.meryl
done

meryl difference hap1.meryl hap2.meryl output hap1.uniq.meryl
meryl difference hap2.meryl hap1.meryl output hap2.uniq.meryl

meryl print threads=16 hap1.uniq.meryl > hap1.uniq.tsv
meryl print threads=16 hap2.uniq.meryl > hap2.uniq.tsv

for hap in hap1 hap2
do
    awk '{if ($2 == 1) print}' ${hap}.uniq.tsv > ${hap}.cnt.uniq.tsv
    gzip ${hap}.cnt.uniq.tsv
    gzip ${hap}.uniq.tsv
done

for hap in hap1 hap2
do
    bam_refiner locate-kmers \
        -i ${WORK_DIR}/meryl/${hap}.cnt.uniq.tsv.gz \
        -f ${hap}_contig \
        -k 21 | sort -k 1,1 -k 2,2n > ${OUTPUT_DIR}/${hap}_cnt_kmerposition.bed
    bgzip -f ${hap}_cnt_kmerposition.bed
    tabix -p bed ${hap}_cnt_kmerposition.bed.gz
done
```
#### Step 2: Align sequencing reads to the diploid genome assembly constructed in Step 1
You can use [minimap2](https://github.com/lh3/minimap2.git) and [samtools](http://www.htslib.org) for alignment.\
You should sort the bam file by read name for [bam_refiner](https://github.com/yos-sk/bam_refiner.git).

```
cat hap1_contig.fa hap2_contig.fa > reference.fa

# ONT
minimap2 -t 16 -ax asm10 reference.fa input.fastq | samtools view --Shb > output.unsorted
# HiFi
minimap2 -t 16 -ax asm5 reference.fa input.fastq | samtools view --Shb > output.unsorted

samtools sort -@ 16 -m 2G -n output.unsorted -o output.bam
samtools index output.bam
```

#### Step 3 (Optional): Split BAM file for the array job of bam_refiner
Please split the BAM file using [split_bam](https://github.com/yos-sk/split_bam.git) if necessary.
```
SIZE=`${path_to_split_bam}/split_bam size --input-file ${INPUT_BAM}`
split_bam split \
    --input-file output.bam \
    --output-dir ${OUTPUT_DIR} \
    --input-size ${SIZE} \
    --num-split 8
```

#### Step 4: Refine bam file

```
./target/release/bam_refiner refine \
    --input-bam ${INPUT_BAM} \
    --output-bam ${OUTPUT_DIR}/${OUTPUT_BAM} \
    --hap1-tabix hap1_cnt_kmerposition.bed.gz \
    --hap2-tabix hap2_cnt_kmerposition.bed.gz \
    --kmer-size 21 \
    1>output.tsv 2>log
```

#### Step 5 (Optional): Merge BAM files if you split bam file in Step 3
```
samtools merge \
        -@ 8 \
        -o output_refined.bam \
        ${OUTPUT_DIR}/*.refined.bam
```

#### Step 6: Sort refined bam file　
```
samtools sort \
    -@ 8 \
    -o output_refined.sorted.bam \
    output_refined.bam 
samtools index output_refined.sorted.bam 
```

If you want to use bam_refiner for the alignment data to the exisitng reference genome (e.g. GRCh38 or CHM13), plase try single_mode branch and see [document](https://github.com/yos-sk/bam_refiner/blob/master/document/single.md).
