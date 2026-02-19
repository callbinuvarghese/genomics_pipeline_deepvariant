#!/bin/bash
# Genomics Pipeline Script for GCP VM - GPU VERSION
# Requires: Deep Learning VM image with pre-installed NVIDIA drivers
# This version does NOT auto-delete and logs errors to GCS

# ====================
# Configuration
# ====================
SAMPLE_ID="SRR1770413"
INPUT_BUCKET="gs://genomics-input-bucket-binuv/ecoli"
OUTPUT_BUCKET="gs://genomics-output-bucket-binuv/ecoli"
REF_BUCKET="gs://genomics-ref-bucket-binuv/ecoli"
THREADS=$(nproc)
SCRATCH_DIR="/mnt/scratch"

# Create scratch directory first (before setting up logging)
mkdir -p $SCRATCH_DIR

LOG_FILE="$SCRATCH_DIR/pipeline_gpu.log"

# ====================
# Error Handling
# ====================
set -euo pipefail  # Exit on error

# Redirect all output to log file AND console
exec > >(tee -a "$LOG_FILE") 2>&1

# Trap errors and upload log before exiting
trap 'ERROR_CODE=$?; echo "=== ERROR at line $LINENO ==="; echo "Exit code: $ERROR_CODE"; gsutil cp $LOG_FILE gs://genomics-ref-bucket-binuv/logs/error_gpu_$(date +%s).log || true; exit $ERROR_CODE' ERR

# ====================
# Setup
# ====================
echo "=== Starting Genomics Pipeline (GPU MODE) ==="
echo "Sample: $SAMPLE_ID"
echo "Start time: $(date)"
echo "Threads: $THREADS"

# Ensure /usr/local/bin is in PATH
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/usr/local/lib"

cd $SCRATCH_DIR
export DEBIAN_FRONTEND=noninteractive

# ====================
# Verify GPU
# ====================
echo "=== Verifying GPU availability ==="
nvidia-smi || { echo "ERROR: nvidia-smi not available. Use Deep Learning VM image!"; exit 1; }
echo "GPU detected ✓"

# ====================
# Install Docker and NVIDIA Container Toolkit
# ====================
echo "=== Installing Docker ==="
apt-get update -qq
apt-get install -y -qq docker.io curl gnupg
systemctl start docker
systemctl enable docker
echo "Docker installed ✓"

# Install NVIDIA Container Toolkit
echo "=== Installing NVIDIA Container Toolkit ==="
# Download and install GPG key (non-interactive)
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | gpg --batch --yes --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg || { echo "ERROR: Failed to download NVIDIA GPG key"; exit 1; }

# Add repository
curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
  sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
  tee /etc/apt/sources.list.d/nvidia-container-toolkit.list || { echo "ERROR: Failed to add NVIDIA repository"; exit 1; }

# Install toolkit
apt-get update -qq
apt-get install -y -qq nvidia-container-toolkit || { echo "ERROR: Failed to install nvidia-container-toolkit"; exit 1; }

# Configure Docker
nvidia-ctk runtime configure --runtime=docker || { echo "ERROR: Failed to configure NVIDIA runtime"; exit 1; }
systemctl restart docker
echo "NVIDIA Container Toolkit installed ✓"

# Verify Docker GPU support
echo "=== Verifying Docker GPU support ==="
docker run --rm --gpus all nvidia/cuda:12.0.0-base-ubuntu22.04 nvidia-smi || { echo "ERROR: Docker GPU runtime not configured"; exit 1; }
echo "Docker GPU support verified ✓"

# ====================
# Verify GCS Access
# ====================
echo "=== Verifying GCS bucket access ==="
gsutil ls $REF_BUCKET/ > /dev/null || { echo "ERROR: Cannot access REF_BUCKET"; exit 1; }
gsutil ls $INPUT_BUCKET/ > /dev/null || { echo "ERROR: Cannot access INPUT_BUCKET"; exit 1; }

# For output bucket, just verify we can write to it (create a test file)
echo "test" | gsutil cp - ${OUTPUT_BUCKET}/.verify_access || { echo "ERROR: Cannot write to OUTPUT_BUCKET"; exit 1; }
gsutil rm ${OUTPUT_BUCKET}/.verify_access || true

echo "GCS buckets accessible ✓"

# ====================
# Install Dependencies
# ====================
echo "=== Installing build dependencies ==="
apt-get install -y -qq wget bzip2 gcc make zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev libncurses5-dev libncursesw5-dev \
  python3 python3-pip libdeflate-dev pigz

echo "Build dependencies installed ✓"

# Install BWA-MEM2 (faster than BWA)
echo "Installing BWA-MEM2..."
wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2
cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/
bwa-mem2 version || { echo "ERROR: BWA-MEM2 installation failed"; exit 1; }

# Install samtools
echo "Installing samtools..."
wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2
tar -xf samtools-1.19.tar.bz2
cd samtools-1.19 && ./configure --with-libdeflate --prefix=/usr/local && make -j$THREADS && make install
cd ..
samtools --version || { echo "ERROR: samtools installation failed"; exit 1; }

# Install bcftools (required for VCF statistics)
echo "Installing bcftools..."
wget -q https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar -xf bcftools-1.19.tar.bz2
cd bcftools-1.19 && ./configure --prefix=/usr/local && make -j$THREADS && make install
cd ..
bcftools --version || { echo "ERROR: bcftools installation failed"; exit 1; }

echo "All tools installed ✓"

# ====================xa
# Download Data
# ====================
echo "=== Downloading reference genome ==="
gsutil -m cp "${REF_BUCKET}/ecoli_k12_mg1655.fna*" . || { echo "ERROR: Failed to download reference"; exit 1; }
gsutil -m cp "${REF_BUCKET}/ecoli_k12_mg1655.dict" . || { echo "ERROR: Failed to download dict"; exit 1; }
ls -lh ecoli_k12_mg1655.fna*

echo "=== Checking BWA-MEM2 index files ==="
# Check if BWA-MEM2 index exists in GCS, if not create it
if gsutil ls "${REF_BUCKET}/ecoli_k12_mg1655.fna.0123" &>/dev/null; then
  echo "BWA-MEM2 index found in GCS, downloading..."
  gsutil -m cp "${REF_BUCKET}/ecoli_k12_mg1655.fna.0123" .
  gsutil -m cp "${REF_BUCKET}/ecoli_k12_mg1655.fna.bwt.2bit.64" .
else
  echo "BWA-MEM2 index not found, generating (one-time setup)..."
  bwa-mem2 index ecoli_k12_mg1655.fna
  echo "Uploading BWA-MEM2 index to GCS for future use..."
  gsutil -m cp ecoli_k12_mg1655.fna.0123 "${REF_BUCKET}/"
  gsutil -m cp ecoli_k12_mg1655.fna.bwt.2bit.64 "${REF_BUCKET}/"
fi

echo "=== Downloading FASTQ files ==="
gsutil -m cp "${INPUT_BUCKET}/${SAMPLE_ID}_1.fastq.gz" . || { echo "ERROR: Failed to download FASTQ R1"; exit 1; }
gsutil -m cp "${INPUT_BUCKET}/${SAMPLE_ID}_2.fastq.gz" . || { echo "ERROR: Failed to download FASTQ R2"; exit 1; }
ls -lh ${SAMPLE_ID}*.fastq.gz

# ====================
# Pipeline Execution
# ====================

# Step 1: Alignment with BWA-MEM2
echo "=== Step 1: Alignment ==="
START_ALIGN=$(date +%s)

/usr/local/bin/bwa-mem2 mem -t $THREADS -K 100000000 \
  -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1" \
  ecoli_k12_mg1655.fna \
  ${SAMPLE_ID}_1.fastq.gz \
  ${SAMPLE_ID}_2.fastq.gz \
  | /usr/local/bin/samtools view -@ $THREADS -m 2G -bS - > ${SAMPLE_ID}.bam

END_ALIGN=$(date +%s)
echo "Alignment completed in $((END_ALIGN - START_ALIGN)) seconds"

# Step 2: Fixmate and Sort (piped to avoid intermediate file)
echo "=== Step 2: Fixmate + Sort (piped) ==="
START_FIXMATE_SORT=$(date +%s)

# Pipe fixmate directly into sort to avoid writing intermediate file
samtools fixmate -@ $THREADS -m ${SAMPLE_ID}.bam - | \
  samtools sort -@ $THREADS -m 2G -o ${SAMPLE_ID}.sorted.bam -

# Clean up intermediate alignment BAM to save disk space
rm -f ${SAMPLE_ID}.bam

END_FIXMATE_SORT=$(date +%s)
echo "Fixmate + Sort completed in $((END_FIXMATE_SORT - START_FIXMATE_SORT)) seconds"

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
echo "=== Step 4: Variant Calling with DeepVariant (GPU) ==="
START_VAR=$(date +%s)

# Pull DeepVariant GPU Docker image
docker pull google/deepvariant:1.6.0-gpu

echo "Running DeepVariant with GPU acceleration..."
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
echo "Fixmate + Sort time: $((END_FIXMATE_SORT - START_FIXMATE_SORT))s"
echo "Dedup time: $((END_DEDUP - START_DEDUP))s"
echo "Variant calling time: $((END_VAR - START_VAR))s"
echo "Total runtime: $((END_VAR - START_ALIGN))s"

# ====================
# Upload Logs & System Info
# ====================
echo "=== Uploading logs and diagnostics to GCS ==="

LOG_BUCKET="gs://genomics-ref-bucket-binuv/logs"
gsutil cp $LOG_FILE "${LOG_BUCKET}/success_gpu_$(basename $LOG_FILE)" || echo "Warning: Failed to upload pipeline log"

# Capture and upload system journal logs
echo "Capturing system logs..."
journalctl -u google-startup-scripts.service --no-pager > $SCRATCH_DIR/startup-script.log || true
gsutil cp $SCRATCH_DIR/startup-script.log "${LOG_BUCKET}/startup-script_gpu_$(date +%s).log" || echo "Warning: Failed to upload startup script log"

echo "Logs uploaded successfully to ${LOG_BUCKET}/"

# Delete scratch files
rm -rf $SCRATCH_DIR/*

echo "=== Pipeline Complete (GPU MODE - VM NOT deleted) ==="
echo "SSH into VM: gcloud compute ssh genomics-pipeline-vm-gpu --zone=us-central1-a"
echo "To delete: gcloud compute instances delete genomics-pipeline-vm-gpu --zone=us-central1-a"
echo "View logs: gsutil ls ${LOG_BUCKET}/"
