#!/bin/bash
# Genomics Pipeline Script for Human Genome - GPU VERSION WITH CHECKPOINTING
# Designed to run on genomics-pipeline-human-ready image with pre-loaded reference
# Automatically deletes VM after completion (production mode)
# CHECKPOINT VERSION: Uploads logs after each step + BAM file after Step 3
# Redirect all output to a local file AND the GCS bucket simultaneously
LOG_FILE="/tmp/pipeline.log"
REMOTE_LOG="gs://genomics-ref-bucket-binuv/logs/v5_live_progress.log"

exec > >(tee -a "$LOG_FILE") 2>&1

# Background loop to sync logs to GCS every 60 seconds
(while true; do gsutil cp "$LOG_FILE" "$REMOTE_LOG"; sleep 60; done) &

# ====================
# Configuration - CUSTOMIZE THESE
# ====================
# Sample identifier (will be used for all output files)
SAMPLE_ID="${SAMPLE_ID:-NA12878}"  # Default to Genome in a Bottle sample

# GCS bucket paths
INPUT_BUCKET="gs://genomics-input-bucket-binuv/human"
OUTPUT_BUCKET="gs://genomics-output-bucket-binuv/human"
CHECKPOINT_BUCKET="gs://genomics-ref-bucket-binuv/checkpoints"

# Reference genome location (pre-loaded in image)
REF_DIR="/opt/genomics/references/human"
REF_FASTA="Homo_sapiens_assembly38.fasta"

# Input FASTQ files
FASTQ_R1="${SAMPLE_ID}_R1.fastq.gz"
FASTQ_R2="${SAMPLE_ID}_R2.fastq.gz"

# Use the THREADS value from metadata if it exists; otherwise, use all available cores
THREADS="${THREADS:-$(nproc)}"
# Fetch THREADS from metadata server if not already set
if [ -z "${THREADS:-}" ]; then
  THREADS=$(curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/attributes/THREADS 2>/dev/null || echo "$(nproc)")
fi
echo "Using $THREADS threads for this run."

# --- NETWORK & SWAP HARDENING ---
echo "Forcing manual networking and adding OOM protection..."

# 1. Add Swap (64GB RAM might be peaking)
fallocate -l 32G /swapfile && chmod 600 /swapfile && mkswap /swapfile && swapon /swapfile

# 2. Manual IP (No DHCP)
INSTANCE_IP=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/network-interfaces/0/ip)

# --- NETWORK HARDENING ---
echo "Optimizing network for high-load genomics..."
# Create a network override for ens5
mkdir -p /etc/systemd/network/
cat <<EOF > /etc/systemd/network/10-ens5.network
[Match]
Name=ens5
[Network]
Address=${INSTANCE_IP}/32
DHCP=no
KeepConfiguration=yes
EOF

systemctl restart systemd-networkd
# -------------------------
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
echo "=== Starting Human Genomics Pipeline (GPU MODE WITH CHECKPOINTING) ==="
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
# Check for Checkpoint (Resume Logic)
# ====================
echo ""
echo "=== Checking for previous checkpoint ==="
RESUME_FROM_CHECKPOINT=false

if gsutil -q stat "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam"; then
  echo "✓ Checkpoint found! Previous run was interrupted."
  echo "  Found: ${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam"

  # Verify checkpoint files are complete
  if gsutil -q stat "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam.bai"; then
    echo "✓ Checkpoint BAM index found"
    echo ""
    echo "RESUMING FROM CHECKPOINT: Skipping Steps 1-3, jumping to Step 4 (Variant Calling)"
    echo "Downloading checkpoint files from GCS..."

    START_CHECKPOINT_DOWNLOAD=$(date +%s)
    gsutil -m cp "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam" . || { echo "ERROR: Failed to download checkpoint BAM"; exit 1; }
    gsutil -m cp "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam.bai" . || { echo "ERROR: Failed to download checkpoint BAI"; exit 1; }
    END_CHECKPOINT_DOWNLOAD=$(date +%s)

    echo "Checkpoint download completed in $((END_CHECKPOINT_DOWNLOAD - START_CHECKPOINT_DOWNLOAD)) seconds"
    ls -lh ${SAMPLE_ID}.sorted.dedup.bam*

    RESUME_FROM_CHECKPOINT=true

    # Set timing variables for steps we're skipping
    START_FASTQ=0
    END_FASTQ=0
    START_ALIGN=0
    END_ALIGN=0
    START_DEDUP=0
    END_DEDUP=0
    START_CHECKPOINT=0
    END_CHECKPOINT=0
  else
    echo "WARNING: Checkpoint BAM found but index is missing"
    echo "WARNING: Checkpoint may be incomplete, will run full pipeline"
    echo "Cleaning up incomplete checkpoint..."
    gsutil rm "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam" || true
  fi
else
  echo "No checkpoint found for sample ${SAMPLE_ID}"
  echo "Will run full pipeline from the beginning"
fi
echo ""

# ====================
# Download FASTQ Files (Skip if resuming from checkpoint)
# ====================
if [ "$RESUME_FROM_CHECKPOINT" = false ]; then
  echo "=== Downloading FASTQ files ==="
  START_FASTQ=$(date +%s)

  gsutil -m cp "${INPUT_BUCKET}/${FASTQ_R1}" . || { echo "ERROR: Failed to download FASTQ R1"; exit 1; }
  gsutil -m cp "${INPUT_BUCKET}/${FASTQ_R2}" . || { echo "ERROR: Failed to download FASTQ R2"; exit 1; }

  END_FASTQ=$(date +%s)
  echo "FASTQ download completed in $((END_FASTQ - START_FASTQ)) seconds"
  ls -lh ${SAMPLE_ID}*.fastq.gz
else
  echo "=== Skipping FASTQ download (resuming from checkpoint) ==="
  START_FASTQ=0
  END_FASTQ=0
fi

# ====================
# Pipeline Execution
# ====================
PIPELINE_START=$(date +%s)
LOG_BUCKET="gs://genomics-ref-bucket-binuv/logs"

# Steps 1-3: Alignment, Sorting, and Duplicate Marking (Skip if resuming from checkpoint)
if [ "$RESUME_FROM_CHECKPOINT" = false ]; then
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

  # CHECKPOINT: Upload deduplicated BAM (enables resume from variant calling)
  echo "=== CHECKPOINT: Uploading deduplicated BAM to GCS ==="
  START_CHECKPOINT=$(date +%s)
  gsutil -m cp ${SAMPLE_ID}.sorted.dedup.bam "${CHECKPOINT_BUCKET}/" || echo "Warning: Failed to upload checkpoint BAM"
  gsutil -m cp ${SAMPLE_ID}.sorted.dedup.bam.bai "${CHECKPOINT_BUCKET}/" || echo "Warning: Failed to upload checkpoint BAI"
  END_CHECKPOINT=$(date +%s)
  echo "Checkpoint uploaded in $((END_CHECKPOINT - START_CHECKPOINT)) seconds"
else
  echo "=== Skipping Steps 1-3 (resuming from checkpoint) ==="
  echo "Using previously computed BAM: ${SAMPLE_ID}.sorted.dedup.bam"
  # Timing variables already set to 0 during checkpoint download
fi

# Step 4: Variant Calling with DeepVariant (GPU-accelerated)
echo "=== Step 4: Variant Calling with DeepVariant (GPU) ==="
START_VAR=$(date +%s)

echo "Running DeepVariant with GPU acceleration..."
echo "Note: make_examples (CPU-intensive) runs first, then call_variants (GPU) and postprocess_variants"
echo "GPU utilization will be low initially, then spike during variant calling phase"

# Getting prempted so adding 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64

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
echo "Checkpoint Upload:   $((END_CHECKPOINT - START_CHECKPOINT))s"
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
echo "Checkpoint location: ${CHECKPOINT_BUCKET}/"
echo "==================================================="

# ====================
# Cleanup Checkpoint (Pipeline completed successfully)
# ====================
echo ""
echo "=== Cleaning up checkpoint files ==="
echo "Pipeline completed successfully - removing checkpoint to save storage costs"
gsutil rm "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam" 2>/dev/null || echo "No checkpoint BAM to remove"
gsutil rm "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam.bai" 2>/dev/null || echo "No checkpoint BAI to remove"
echo "Checkpoint cleanup complete"

# ====================
# Upload Logs
# ====================
echo ""
echo "=== Uploading logs to GCS ==="

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


# ➜  genomics gcloud compute instances create genomics-pipeline-human-vm-final-v4 \
#   --project=genomics-pipeline-poc \
#   --zone=us-central1-a \
#   --machine-type=g2-standard-16 \
#   --network-interface=nic-type=GVNIC,private-network-ip=10.128.0.32 \
#   --image-family=genomics-pipeline-human-ready \
#   --local-ssd=interface=NVME \
#   --maintenance-policy=TERMINATE \
#   --preemptible \
#   --scopes=https://www.googleapis.com/auth/cloud-platform \
#   --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
#   --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=NA12878,THREADS=12

# ➜  genomics gsutil ls -l gs://genomics-ref-bucket-binuv/checkpoints/
# 11820939142  2026-02-16T00:51:56Z  gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam
#    4504480  2026-02-16T00:51:58Z  gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam.bai
# TOTAL: 2 objects, 11825443622 bytes (11.01 GiB)


# ➜  genomics gcloud compute operations list \
#     --filter="targetLink:genomics-pipeline-human-vm-final-v4" \
#     --sort-by=~startTime
# NAME                                                     TYPE    TARGET                                                       HTTP_STATUS  STATUS  TIMESTAMP
# operation-1771203139931-64ae65dfb0f11-4b7d0e9a-672e7b4c  delete  us-central1-a/instances/genomics-pipeline-human-vm-final-v4  200          DONE    2026-02-15T16:52:20.130-08:00
# operation-1771198742581-64ae557e0d1d5-72607f8f-37ebc8b5  insert  us-central1-a/instances/genomics-pipeline-human-vm-final-v4  200          DONE    2026-02-15T15:39:04.051-08:00
# ➜  genomics
