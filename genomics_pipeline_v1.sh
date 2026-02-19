#!/bin/bash
# Genomics Pipeline Script for GCP VM
# Processes FASTQ → Sorted BAM → Variants

set -euo pipefail  # Exit on error

# ====================
# Configuration
# ====================
SAMPLE_ID="SRR1770413"
INPUT_BUCKET="gs://genomics-input-bucket-binuv/ecoli"
OUTPUT_BUCKET="gs://genomics-output-bucket-binuv/ecoli"
REF_BUCKET="gs://genomics-ref-bucket-binuv/ecoli"
THREADS=4
SCRATCH_DIR="/mnt/scratch"

# ====================
# Setup
# ====================
echo "=== Starting Genomics Pipeline ==="
echo "Sample: $SAMPLE_ID"
echo "Start time: $(date)"

# Create scratch directory
mkdir -p $SCRATCH_DIR
cd $SCRATCH_DIR

# ====================
# Install Dependencies
# ====================
echo "=== Installing dependencies ==="
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq wget bzip2 gcc make zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev python3 python3-pip docker.io

# Install BWA-MEM2 (faster than BWA)
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2
cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/

# Install samtools
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xf samtools-1.18.tar.bz2
cd samtools-1.18 && ./configure --prefix=/usr/local && make && make install && cd ..

# ====================
# Download Data
# ====================
echo "=== Downloading reference genome ==="
gsutil -m cp "${REF_BUCKET}/ecoli_k12_mg1655.fna*" .
gsutil -m cp "${REF_BUCKET}/ecoli_k12_mg1655.dict" .

echo "=== Downloading FASTQ files ==="
gsutil -m cp "${INPUT_BUCKET}/${SAMPLE_ID}_1.fastq.gz" .
gsutil -m cp "${INPUT_BUCKET}/${SAMPLE_ID}_2.fastq.gz" .

# ====================
# Pipeline Execution
# ====================

# Step 1: Alignment with BWA-MEM2
echo "=== Step 1: Alignment ==="
START_ALIGN=$(date +%s)

bwa-mem2 mem -t $THREADS \
  -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1" \
  ecoli_k12_mg1655.fna \
  ${SAMPLE_ID}_1.fastq.gz \
  ${SAMPLE_ID}_2.fastq.gz \
  | samtools view -@ $THREADS -bS - > ${SAMPLE_ID}.bam

END_ALIGN=$(date +%s)
echo "Alignment completed in $((END_ALIGN - START_ALIGN)) seconds"

# Step 2: Sort BAM
echo "=== Step 2: Sorting BAM ==="
START_SORT=$(date +%s)

samtools sort -@ $THREADS -o ${SAMPLE_ID}.sorted.bam ${SAMPLE_ID}.bam
samtools index ${SAMPLE_ID}.sorted.bam

END_SORT=$(date +%s)
echo "Sorting completed in $((END_SORT - START_SORT)) seconds"

# Step 3: Mark Duplicates
echo "=== Step 3: Marking Duplicates ==="
START_DEDUP=$(date +%s)

samtools markdup -@ $THREADS \
  ${SAMPLE_ID}.sorted.bam \
  ${SAMPLE_ID}.sorted.dedup.bam

samtools index ${SAMPLE_ID}.sorted.dedup.bam

END_DEDUP=$(date +%s)
echo "Duplicate marking completed in $((END_DEDUP - START_DEDUP)) seconds"

# Step 4: Variant Calling with DeepVariant (GPU-accelerated)
echo "=== Step 4: Variant Calling with DeepVariant ==="
START_VAR=$(date +%s)

# Pull DeepVariant Docker image
docker pull google/deepvariant:1.6.0-gpu

# Run DeepVariant on GPU
docker run --gpus all \
  -v "$SCRATCH_DIR":"$SCRATCH_DIR" \
  google/deepvariant:1.6.0-gpu \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=$SCRATCH_DIR/ecoli_k12_mg1655.fna \
  --reads=$SCRATCH_DIR/${SAMPLE_ID}.sorted.dedup.bam \
  --output_vcf=$SCRATCH_DIR/${SAMPLE_ID}.vcf.gz \
  --output_gvcf=$SCRATCH_DIR/${SAMPLE_ID}.g.vcf.gz \
  --num_shards=$THREADS

END_VAR=$(date +%s)
echo "Variant calling completed in $((END_VAR - START_VAR)) seconds"

# ====================
# Quality Metrics
# ====================
echo "=== Generating Quality Metrics ==="

# BAM statistics
samtools flagstat ${SAMPLE_ID}.sorted.dedup.bam > ${SAMPLE_ID}.flagstat.txt
samtools stats ${SAMPLE_ID}.sorted.dedup.bam > ${SAMPLE_ID}.stats.txt

# VCF statistics
bcftools stats ${SAMPLE_ID}.vcf.gz > ${SAMPLE_ID}.vcf.stats.txt

# ====================
# Upload Results
# ====================
echo "=== Uploading results to GCS ==="
gsutil -m cp ${SAMPLE_ID}.sorted.dedup.bam "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.sorted.dedup.bam.bai "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.vcf.gz "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.g.vcf.gz "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.*.txt "${OUTPUT_BUCKET}/"

# ====================
# Cleanup & Summary
# ====================
echo "=== Pipeline Summary ==="
echo "Sample: $SAMPLE_ID"
echo "End time: $(date)"
echo "Alignment time: $((END_ALIGN - START_ALIGN))s"
echo "Sorting time: $((END_SORT - START_SORT))s"
echo "Dedup time: $((END_DEDUP - START_DEDUP))s"
echo "Variant calling time: $((END_VAR - START_VAR))s"
echo "Total runtime: $((END_VAR - START_ALIGN))s"

# Delete scratch files
rm -rf $SCRATCH_DIR/*

# ====================
# Auto-Delete VM (Cost Saving!)
# ====================
echo "=== Auto-deleting VM ==="
ZONE=$(curl -H "Metadata-Flavor: Google" \
  http://metadata.google.internal/computeMetadata/v1/instance/zone | cut -d'/' -f4)
INSTANCE_NAME=$(hostname)

gcloud compute instances delete $INSTANCE_NAME \
  --zone=$ZONE \
  --quiet \
  --delete-disks=all

echo "=== Pipeline Complete ==="
