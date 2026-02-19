#!/bin/bash
# Genomics Pipeline Script for Human Genome - GPU VERSION
# Designed to run on genomics-pipeline-human-ready image with pre-loaded reference
# Automatically deletes VM after completion (production mode)

# ====================
# Configuration - CUSTOMIZE THESE
# ====================
# Sample identifier (will be used for all output files)
SAMPLE_ID="${SAMPLE_ID:-NA12878}"  # Default to Genome in a Bottle sample

# GCS bucket paths
INPUT_BUCKET="gs://genomics-input-bucket-binuv/human"
OUTPUT_BUCKET="gs://genomics-output-bucket-binuv/human"

# Reference genome location (pre-loaded in image)
REF_DIR="/opt/genomics/references/human"
REF_FASTA="Homo_sapiens_assembly38.fasta"

# Input FASTQ files
FASTQ_R1="${SAMPLE_ID}_R1.fastq.gz"
FASTQ_R2="${SAMPLE_ID}_R2.fastq.gz"

# Compute resources
THREADS=$(nproc)

# Calculate safe memory allocation for samtools sort
# Rule: samtools uses -m per thread, so total = THREADS × MEM_PER_THREAD
# For n1-standard-16: 60 GB RAM / 16 threads = 3.75 GB per thread
# Use 2 GB to leave headroom for OS and other processes
TOTAL_RAM_GB=$(free -g | awk '/^Mem:/{print $2}')
SAFE_MEM_PER_THREAD=$((TOTAL_RAM_GB / THREADS - 1))
# Clamp between 1G and 4G
if [ "$SAFE_MEM_PER_THREAD" -lt 1 ]; then
  SAFE_MEM_PER_THREAD=1
elif [ "$SAFE_MEM_PER_THREAD" -gt 4 ]; then
  SAFE_MEM_PER_THREAD=4
fi
SAMTOOLS_MEM="${SAFE_MEM_PER_THREAD}G"

echo "Detected $THREADS threads and ${TOTAL_RAM_GB}GB RAM"
echo "Using ${SAMTOOLS_MEM} per thread for samtools sort (total: $((THREADS * SAFE_MEM_PER_THREAD))GB)"

# Scratch directory - will auto-detect Local SSD or use boot disk
SCRATCH_DIR="/mnt/disks/scratch"  # Default to Local SSD mount point
# Fallback to boot disk if Local SSD not available
[ -d "/mnt/disks/scratch" ] || SCRATCH_DIR="/mnt/scratch"

# Disk space requirements (in GB)
MIN_SCRATCH_SPACE_GB=100  # Minimum free space needed

# ====================
# Setup Logging
# ====================
mkdir -p $SCRATCH_DIR
LOG_FILE="$SCRATCH_DIR/pipeline_human_gpu.log"

# Redirect all output to log file AND console
exec > >(tee -a "$LOG_FILE") 2>&1

# ====================
# Helper Functions
# ====================
# Function to check and report disk space
check_disk_space() {
  local mount_point=$1
  local min_gb=$2
  local available_gb=$(df -BG "$mount_point" | awk 'NR==2 {print int($4)}')
  local total_gb=$(df -BG "$mount_point" | awk 'NR==2 {print int($2)}')
  local used_percent=$(df "$mount_point" | awk 'NR==2 {print $5}')

  echo "Disk space on $mount_point: ${available_gb}GB available / ${total_gb}GB total (${used_percent} used)"

  if [ "$available_gb" -lt "$min_gb" ]; then
    echo "ERROR: Insufficient disk space. Need ${min_gb}GB, have ${available_gb}GB"
    return 1
  fi
  return 0
}

# ====================
# Error Handling
# ====================
set -euo pipefail  # Exit on error

# Trap errors: upload log to GCS, then delete VM
trap 'ERROR_CODE=$?; echo "=== ERROR at line $LINENO ==="; echo "Exit code: $ERROR_CODE"; gsutil cp $LOG_FILE gs://genomics-ref-bucket-binuv/logs/error_human_gpu_$(date +%s).log || true; gcloud compute instances delete $(hostname) --zone=$(curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone | cut -d/ -f4) --quiet || true; exit $ERROR_CODE' ERR

# ====================
# Pipeline Start
# ====================
echo "=== Starting Human Genomics Pipeline (GPU MODE) ==="
echo "Sample: $SAMPLE_ID"
echo "Start time: $(date)"
echo "Threads: $THREADS"
echo "Reference: $REF_FASTA (pre-loaded)"

# Ensure paths are set
export PATH="/usr/local/bin:$PATH"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/usr/local/lib"

cd $SCRATCH_DIR

# ====================
# Setup GCP Ops Agent (Pre-installed in custom image)
# ====================
echo "=== Setting up GCP Ops Agent ==="
if systemctl is-active --quiet google-cloud-ops-agent 2>/dev/null; then
  echo "Ops Agent detected (pre-installed in image)"
  echo "Restarting to generate fresh VM identity..."
  sudo systemctl restart google-cloud-ops-agent
  echo "Ops Agent restarted with new VM identity ✓"
  echo "View metrics at: https://console.cloud.google.com/monitoring/dashboards"
else
  echo "WARNING: Ops Agent not running (expected to be pre-installed in image)"
  echo "Attempting to start..."
  sudo systemctl start google-cloud-ops-agent || echo "Failed to start Ops Agent"
fi

# ====================
# Setup Local SSD (if available)
# ====================
echo "=== Setting up scratch disk ==="
if [ -b /dev/nvme0n1 ]; then
  echo "Local NVMe SSD detected - formatting and mounting..."
  sudo mkfs.ext4 -F /dev/nvme0n1
  sudo mkdir -p /mnt/disks/scratch
  sudo mount /dev/nvme0n1 /mnt/disks/scratch
  sudo chmod 777 /mnt/disks/scratch
  SCRATCH_DIR="/mnt/disks/scratch"
  mkdir -p $SCRATCH_DIR
  LOG_FILE="$SCRATCH_DIR/pipeline_human_gpu.log"
  # Redirect logs to new location
  exec > >(tee -a "$LOG_FILE") 2>&1
  echo "Local SSD mounted at $SCRATCH_DIR ✓"
  echo "Disk performance: NVMe (170K-660K IOPS)"
elif [ -d /mnt/disks/scratch ]; then
  echo "Local SSD already mounted at /mnt/disks/scratch ✓"
  SCRATCH_DIR="/mnt/disks/scratch"
else
  echo "WARNING: No Local SSD detected. Using boot disk at /mnt/scratch"
  echo "WARNING: Performance will be significantly slower (10-20x)"
  echo "WARNING: Recommend recreating VM with --local-ssd interface=NVME"
  SCRATCH_DIR="/mnt/scratch"
  mkdir -p $SCRATCH_DIR
fi

# Verify disk space
echo ""
echo "=== Checking disk space ==="
check_disk_space "$SCRATCH_DIR" "$MIN_SCRATCH_SPACE_GB" || exit 1
echo "Disk space check passed ✓"
echo ""

# ====================
# Verify Pre-loaded Reference
# ====================
echo "=== Verifying pre-loaded reference genome ==="
if [ ! -f "${REF_DIR}/${REF_FASTA}" ]; then
  echo "ERROR: Reference genome not found at ${REF_DIR}/${REF_FASTA}"
  echo "This script requires the genomics-pipeline-human-ready image with pre-loaded reference"
  exit 1
fi

# Verify BWA-MEM2 index files exist
if [ ! -f "${REF_DIR}/${REF_FASTA}.0123" ]; then
  echo "ERROR: BWA-MEM2 index not found at ${REF_DIR}/${REF_FASTA}.0123"
  echo "This script requires the genomics-pipeline-human-ready image with pre-generated BWA-MEM2 index"
  exit 1
fi

echo "Reference genome found: ${REF_DIR}/${REF_FASTA}"
echo "BWA-MEM2 index found: ${REF_DIR}/${REF_FASTA}.0123"
ls -lh ${REF_DIR}/${REF_FASTA}*

# Create symlinks in scratch dir for easier access
ln -s ${REF_DIR}/${REF_FASTA} .
ln -s ${REF_DIR}/${REF_FASTA}.fai .
ln -s ${REF_DIR}/${REF_FASTA}.dict .
ln -s ${REF_DIR}/${REF_FASTA}.0123 .
ln -s ${REF_DIR}/${REF_FASTA}.amb .
ln -s ${REF_DIR}/${REF_FASTA}.ann .
ln -s ${REF_DIR}/${REF_FASTA}.bwt.2bit.64 .
ln -s ${REF_DIR}/${REF_FASTA}.pac .

echo "Reference verification complete ✓"

# ====================
# Verify GPU and Docker
# ====================
echo "=== Verifying GPU and Docker ==="
nvidia-smi || { echo "ERROR: GPU not available"; exit 1; }
docker run --rm --gpus all nvidia/cuda:12.0.0-base-ubuntu22.04 nvidia-smi || { echo "ERROR: Docker GPU runtime not configured"; exit 1; }
echo "GPU and Docker verified ✓"

# ====================
# Verify GCS Access
# ====================
echo "=== Verifying GCS bucket access ==="
gsutil ls $INPUT_BUCKET/ > /dev/null || { echo "ERROR: Cannot access INPUT_BUCKET"; exit 1; }
echo "test" | gsutil cp - ${OUTPUT_BUCKET}/.verify_access || { echo "ERROR: Cannot write to OUTPUT_BUCKET"; exit 1; }
gsutil rm ${OUTPUT_BUCKET}/.verify_access || true
echo "GCS buckets accessible ✓"

# ====================
# Download FASTQ Files
# ====================
echo "=== Downloading FASTQ files ==="
START_FASTQ=$(date +%s)

gsutil -m cp "${INPUT_BUCKET}/${FASTQ_R1}" . || { echo "ERROR: Failed to download FASTQ R1"; exit 1; }
gsutil -m cp "${INPUT_BUCKET}/${FASTQ_R2}" . || { echo "ERROR: Failed to download FASTQ R2"; exit 1; }

END_FASTQ=$(date +%s)
echo "FASTQ download completed in $((END_FASTQ - START_FASTQ)) seconds"
ls -lh ${SAMPLE_ID}*.fastq.gz

# ====================
# Pipeline Execution
# ====================
PIPELINE_START=$(date +%s)
LOG_BUCKET="gs://genomics-ref-bucket-binuv/logs"

# Steps 1-2 Combined: Alignment → Name Sort → Fixmate → Coordinate Sort (FULLY PIPED)
echo "=== Steps 1-2: Alignment → Sort → Fixmate → Sort (piped) ==="
check_disk_space "$SCRATCH_DIR" 50 || { echo "ERROR: Insufficient space for alignment"; exit 1; }
START_ALIGN=$(date +%s)

echo "Running fully-piped alignment and sorting pipeline..."
echo "Pipeline: BWA-MEM2 → samtools sort -n → samtools fixmate → samtools sort"
echo "Benefits: No intermediate files, reduced I/O, faster processing"

/usr/local/bin/bwa-mem2 mem -t $THREADS -K 100000000 \
  -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1" \
  $REF_FASTA \
  $FASTQ_R1 \
  $FASTQ_R2 \
  | samtools sort -n -@ $THREADS -m $SAMTOOLS_MEM - \
  | samtools fixmate -@ $THREADS -m - - \
  | samtools sort -@ $THREADS -m $SAMTOOLS_MEM -o ${SAMPLE_ID}.sorted.bam -

END_ALIGN=$(date +%s)
echo "Fully-piped alignment and sorting completed in $((END_ALIGN - START_ALIGN)) seconds"
echo "Output: ${SAMPLE_ID}.sorted.bam (coordinate-sorted, fixmate applied)"
check_disk_space "$SCRATCH_DIR" 30

# Upload progress log
gsutil cp $LOG_FILE "${LOG_BUCKET}/progress_step1-2_${SAMPLE_ID}_$(date +%s).log" || echo "Warning: Failed to upload Step 1-2 progress log"

# Clean up FASTQ files to save disk space
rm -f ${FASTQ_R1} ${FASTQ_R2}
echo "FASTQ files removed to save disk space"

# Step 3: Mark Duplicates
echo "=== Step 3: Marking Duplicates ==="
START_DEDUP=$(date +%s)

samtools markdup -@ $THREADS \
  ${SAMPLE_ID}.sorted.bam \
  ${SAMPLE_ID}.sorted.dedup.bam

samtools index ${SAMPLE_ID}.sorted.dedup.bam

rm -f ${SAMPLE_ID}.sorted.bam
echo "Sorted BAM (pre-dedup) removed to save disk space"

END_DEDUP=$(date +%s)
echo "Duplicate marking completed in $((END_DEDUP - START_DEDUP)) seconds"

# Upload progress log
gsutil cp $LOG_FILE "${LOG_BUCKET}/progress_step3_${SAMPLE_ID}_$(date +%s).log" || echo "Warning: Failed to upload Step 3 progress log"

# Step 4: Variant Calling with DeepVariant (GPU-accelerated)
echo "=== Step 4: Variant Calling with DeepVariant (GPU) ==="
START_VAR=$(date +%s)

echo "Running DeepVariant with GPU acceleration..."
echo "Note: make_examples (CPU-intensive) runs first, then call_variants (GPU) and postprocess_variants"
echo "GPU utilization will be low initially, then spike during variant calling phase"

# Mount both scratch dir (for BAM) and reference dir (for FASTA)
# This ensures Docker can access the reference genome even if it's symlinked
docker run --gpus all \
  -v "$SCRATCH_DIR":"$SCRATCH_DIR" \
  -v "$REF_DIR":"$REF_DIR":ro \
  google/deepvariant:1.6.0-gpu \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=$REF_DIR/$REF_FASTA \
  --reads=$SCRATCH_DIR/${SAMPLE_ID}.sorted.dedup.bam \
  --output_vcf=$SCRATCH_DIR/${SAMPLE_ID}.vcf.gz \
  --output_gvcf=$SCRATCH_DIR/${SAMPLE_ID}.g.vcf.gz \
  --num_shards=$THREADS

END_VAR=$(date +%s)
echo "Variant calling completed in $((END_VAR - START_VAR)) seconds"

# Upload progress log
gsutil cp $LOG_FILE "${LOG_BUCKET}/progress_step4_${SAMPLE_ID}_$(date +%s).log" || echo "Warning: Failed to upload Step 4 progress log"

# ====================
# Quality Metrics
# ====================
echo "=== Generating Quality Metrics ==="
START_QC=$(date +%s)

# BAM statistics
samtools flagstat ${SAMPLE_ID}.sorted.dedup.bam > ${SAMPLE_ID}.flagstat.txt
samtools stats ${SAMPLE_ID}.sorted.dedup.bam > ${SAMPLE_ID}.stats.txt

# VCF statistics
bcftools stats ${SAMPLE_ID}.vcf.gz > ${SAMPLE_ID}.vcf.stats.txt

# Extract key metrics for summary
TOTAL_READS=$(grep "^SN.*total.*:" ${SAMPLE_ID}.stats.txt | awk '{print $4}')
MAPPED_READS=$(grep "^SN.*mapped.*:" ${SAMPLE_ID}.stats.txt | awk '{print $4}')
DUPLICATES=$(grep "^SN.*duplicates.*:" ${SAMPLE_ID}.stats.txt | awk '{print $4}')
TOTAL_VARIANTS=$(bcftools view -H ${SAMPLE_ID}.vcf.gz | wc -l)

END_QC=$(date +%s)
echo "QC metrics generated in $((END_QC - START_QC)) seconds"

# ====================
# Upload Results
# ====================
echo "=== Uploading results to GCS ==="
START_UPLOAD=$(date +%s)

gsutil -m cp ${SAMPLE_ID}.sorted.dedup.bam "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.sorted.dedup.bam.bai "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.vcf.gz "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.vcf.gz.tbi "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.g.vcf.gz "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.g.vcf.gz.tbi "${OUTPUT_BUCKET}/"
gsutil -m cp ${SAMPLE_ID}.*.txt "${OUTPUT_BUCKET}/"

END_UPLOAD=$(date +%s)
echo "Upload completed in $((END_UPLOAD - START_UPLOAD)) seconds"

# ====================
# Pipeline Summary
# ====================
TOTAL_TIME=$((END_UPLOAD - PIPELINE_START))

echo ""
echo "==================================================="
echo "===        PIPELINE SUMMARY                     ==="
echo "==================================================="
echo "Sample:              $SAMPLE_ID"
echo "Reference:           $REF_FASTA (pre-loaded)"
echo "End time:            $(date)"
echo ""
echo "--- Timing (seconds) ---"
echo "FASTQ download:      $((END_FASTQ - START_FASTQ))s"
echo "Align+Sort+Fixmate:  $((END_ALIGN - START_ALIGN))s (fully piped)"
echo "Mark Duplicates:     $((END_DEDUP - START_DEDUP))s"
echo "Variant Calling:     $((END_VAR - START_VAR))s"
echo "QC Metrics:          $((END_QC - START_QC))s"
echo "Upload:              $((END_UPLOAD - START_UPLOAD))s"
echo "TOTAL RUNTIME:       ${TOTAL_TIME}s ($(($TOTAL_TIME / 60)) minutes)"
echo ""
echo "--- Quality Metrics ---"
echo "Total reads:         ${TOTAL_READS}"
echo "Mapped reads:        ${MAPPED_READS}"
echo "Duplicates:          ${DUPLICATES}"
echo "Total variants:      ${TOTAL_VARIANTS}"
echo ""
echo "Output location:     ${OUTPUT_BUCKET}/"
echo "==================================================="

# ====================
# Upload Logs
# ====================
echo ""
echo "=== Uploading logs to GCS ==="
LOG_BUCKET="gs://genomics-ref-bucket-binuv/logs"
CHECKPOINT_BUCKET="gs://genomics-ref-bucket-binuv/checkpoints"

gsutil cp $LOG_FILE "${LOG_BUCKET}/success_human_gpu_${SAMPLE_ID}_$(date +%s).log" || echo "Warning: Failed to upload pipeline log"

journalctl -u google-startup-scripts.service --no-pager > $SCRATCH_DIR/startup-script.log || true
gsutil cp $SCRATCH_DIR/startup-script.log "${LOG_BUCKET}/startup_human_gpu_${SAMPLE_ID}_$(date +%s).log" || echo "Warning: Failed to upload startup script log"

echo "Logs uploaded to ${LOG_BUCKET}/"

# ====================
# Cleanup and Auto-Delete VM
# ====================
echo ""
echo "=== Cleaning up scratch files ==="
rm -rf $SCRATCH_DIR/*

echo ""
echo "=== PIPELINE COMPLETE - Deleting VM ==="
ZONE=$(curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone | cut -d/ -f4)
VM_NAME=$(hostname)

echo "Deleting VM: $VM_NAME in zone: $ZONE"
gcloud compute instances delete $VM_NAME --zone=$ZONE --quiet || echo "WARNING: Failed to auto-delete VM"

echo "=== END ==="
