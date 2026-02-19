#!/bin/bash
# Genomics Pipeline Script for GCP VM - DEBUG VERSION
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

LOG_FILE="$SCRATCH_DIR/pipeline.log"

# ====================
# Error Handling
# ====================
set -euo pipefail  # Exit on error

# Redirect all output to log file AND console
exec > >(tee -a "$LOG_FILE") 2>&1

# Trap errors and upload log before exiting
trap 'ERROR_CODE=$?; echo "=== ERROR at line $LINENO ==="; echo "Exit code: $ERROR_CODE"; gsutil cp $LOG_FILE gs://genomics-ref-bucket-binuv/logs/error_$(date +%s).log || true; exit $ERROR_CODE' ERR

# ====================
# Setup
# ====================
echo "=== Starting Genomics Pipeline (DEBUG MODE) ==="
echo "Sample: $SAMPLE_ID"
echo "Start time: $(date)"
echo "Threads: $THREADS"

# Ensure /usr/local/bin is in PATH
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/usr/local/lib"

cd $SCRATCH_DIR
export DEBIAN_FRONTEND=noninteractive

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
echo "=== Installing dependencies ==="
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y -qq wget bzip2 gcc make zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev libncurses5-dev libncursesw5-dev \
  python3 python3-pip docker.io curl gnupg
apt-get install -y -qq libdeflate-dev pigz

# Install NVIDIA GPU Drivers
echo "Installing NVIDIA GPU drivers..."
# Add NVIDIA driver repository
apt-get install -y -qq ubuntu-drivers-common
ubuntu-drivers devices
# Install recommended NVIDIA driver
apt-get install -y -qq nvidia-driver-535
# Load the driver module
modprobe nvidia || true
echo "NVIDIA drivers installed ✓"

# Verify NVIDIA driver is loaded
echo "Verifying NVIDIA driver..."
nvidia-smi || { echo "WARNING: nvidia-smi not working, may need reboot"; }

# Install NVIDIA Container Toolkit for GPU support in Docker
echo "Installing NVIDIA Container Toolkit..."
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
  sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
  tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
apt-get update -qq
apt-get install -y -qq nvidia-container-toolkit
nvidia-ctk runtime configure --runtime=docker
systemctl restart docker
echo "NVIDIA Container Toolkit installed ✓"

# Wait for NVIDIA GPU to be ready
echo "Waiting for NVIDIA GPU drivers..."
MAX_WAIT=300  # 5 minutes
WAIT_TIME=0
while ! nvidia-smi &>/dev/null; do
  if [ $WAIT_TIME -ge $MAX_WAIT ]; then
    echo "ERROR: NVIDIA drivers not available after ${MAX_WAIT} seconds"
    echo "Check driver installation: journalctl -u google-startup-scripts.service"
    exit 1
  fi
  echo "Waiting for nvidia-smi... ($WAIT_TIME/${MAX_WAIT}s)"
  sleep 10
  WAIT_TIME=$((WAIT_TIME + 10))
done
echo "NVIDIA drivers loaded ✓"
nvidia-smi

# Verify GPU is accessible in Docker
echo "Verifying GPU access in Docker..."
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi || { echo "ERROR: GPU not accessible in Docker"; exit 1; }
echo "GPU accessible in Docker ✓"

# Install BWA-MEM2 (faster than BWA)
echo "Installing BWA-MEM2..."
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2
cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/
bwa-mem2 version || { echo "ERROR: BWA-MEM2 installation failed"; exit 1; }

# Install samtools
echo "Installing samtools..."
wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2
tar -xf  samtools-1.19.tar.bz2
cd samtools-1.19 && ./configure --with-libdeflate --prefix=/usr/local && make -j$THREADS && make install
cd ..
samtools --version || { echo "ERROR: samtools installation failed"; exit 1; }

# Install bcftools (required for VCF statistics)
echo "Installing bcftools..."
wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar -xf bcftools-1.19.tar.bz2
cd bcftools-1.19 && ./configure --prefix=/usr/local && make -j$THREADS && make install
cd ..
bcftools --version || { echo "ERROR: bcftools installation failed"; exit 1; }

echo "All dependencies installed ✓"

# Verify tools are in PATH
echo "=== Verifying tool paths ==="
which bwa-mem2 || { echo "ERROR: bwa-mem2 not in PATH"; exit 1; }
which samtools || { echo "ERROR: samtools not in PATH"; exit 1; }
which bcftools || { echo "ERROR: bcftools not in PATH"; exit 1; }
which docker || { echo "ERROR: docker not in PATH"; exit 1; }
echo "All tools accessible ✓"

# ====================
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
  # BWA-MEM2 reuses some BWA index files that are already in GCS
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

# Use full path to ensure samtools is found in the pipe
/usr/local/bin/bwa-mem2 mem -t $THREADS -K 100000000 \
  -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1" \
  ecoli_k12_mg1655.fna \
  ${SAMPLE_ID}_1.fastq.gz \
  ${SAMPLE_ID}_2.fastq.gz \
  | /usr/local/bin/samtools view -@ $THREADS -m 2G -bS - > ${SAMPLE_ID}.bam

END_ALIGN=$(date +%s)
echo "Alignment completed in $((END_ALIGN - START_ALIGN)) seconds"

# Step 2: Fixmate (required for markdup)
echo "=== Step 2: Fixmate ==="
START_FIXMATE=$(date +%s)

samtools fixmate -@ $THREADS -m ${SAMPLE_ID}.bam ${SAMPLE_ID}.fixmate.bam

END_FIXMATE=$(date +%s)
echo "Fixmate completed in $((END_FIXMATE - START_FIXMATE)) seconds"

# Step 3: Sort BAM
echo "=== Step 3: Sorting BAM ==="
START_SORT=$(date +%s)

samtools sort -@ $THREADS -o ${SAMPLE_ID}.sorted.bam ${SAMPLE_ID}.fixmate.bam

END_SORT=$(date +%s)
echo "Sorting completed in $((END_SORT - START_SORT)) seconds"

# Step 4: Mark Duplicates
echo "=== Step 4: Marking Duplicates ==="
START_DEDUP=$(date +%s)

samtools markdup -@ $THREADS \
  ${SAMPLE_ID}.sorted.bam \
  ${SAMPLE_ID}.sorted.dedup.bam

samtools index ${SAMPLE_ID}.sorted.dedup.bam

END_DEDUP=$(date +%s)
echo "Duplicate marking completed in $((END_DEDUP - START_DEDUP)) seconds"

# Step 5: Variant Calling with DeepVariant (GPU-accelerated)
echo "=== Step 5: Variant Calling with DeepVariant ==="
START_VAR=$(date +%s)

# Pull DeepVariant Docker image
docker pull google/deepvariant:1.6.0-gpu

echo "Docker run started "
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
echo "Fixmate time: $((END_FIXMATE - START_FIXMATE))s"
echo "Sorting time: $((END_SORT - START_SORT))s"
echo "Dedup time: $((END_DEDUP - START_DEDUP))s"
echo "Variant calling time: $((END_VAR - START_VAR))s"
echo "Total runtime: $((END_VAR - START_ALIGN))s"

# ====================
# Upload Logs & System Info
# ====================
echo "=== Uploading logs and diagnostics to GCS ==="

# Upload pipeline log
LOG_BUCKET="gs://genomics-ref-bucket-binuv/logs"
gsutil cp $LOG_FILE "${LOG_BUCKET}/success_$(basename $LOG_FILE)" || echo "Warning: Failed to upload pipeline log"

# Capture and upload system journal logs
echo "Capturing system logs..."
journalctl -u google-startup-scripts.service --no-pager > $SCRATCH_DIR/startup-script.log || true
gsutil cp $SCRATCH_DIR/startup-script.log "${LOG_BUCKET}/startup-script_$(date +%s).log" || echo "Warning: Failed to upload startup script log"

# Upload full system journal
journalctl --no-pager > $SCRATCH_DIR/full-system.log || true
gsutil cp $SCRATCH_DIR/full-system.log "${LOG_BUCKET}/system_$(date +%s).log" || echo "Warning: Failed to upload system log"

echo "Logs uploaded successfully to ${LOG_BUCKET}/"

# Delete scratch files
rm -rf $SCRATCH_DIR/*

echo "=== Pipeline Complete (DEBUG MODE - VM NOT deleted) ==="
echo "SSH into VM to investigate: gcloud compute ssh genomics-pipeline-vm-debug --zone=us-central1-a"
echo "To manually delete: gcloud compute instances delete genomics-pipeline-vm-debug --zone=us-central1-a"
echo "View logs: gsutil ls ${LOG_BUCKET}/"

# ====================
# VM remains running for debugging
# ====================
