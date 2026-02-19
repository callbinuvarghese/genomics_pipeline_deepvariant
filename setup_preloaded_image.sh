#!/bin/bash
# Startup script to create a pre-loaded genomics image
# This installs all tools, downloads reference genome, and generates BWA-MEM2 index
# RECOMMENDED: n1-standard-16 (60GB RAM + 20GB swap) or n1-highmem-16 (104GB RAM)
# MINIMUM: n1-standard-8 (30GB RAM + 20GB swap) - will be slower

set -e
export DEBIAN_FRONTEND=noninteractive

# Setup logging
LOG_FILE="/var/log/genomics-setup.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== Starting genomics image setup at $(date) ==="

# Install dependencies
echo "=== Installing dependencies ==="
apt-get update -qq
apt-get install -y -qq wget bzip2 gcc make zlib1g-dev libbz2-dev liblzma-dev \
  libcurl4-openssl-dev libssl-dev libncurses5-dev libncursesw5-dev \
  python3 python3-pip libdeflate-dev pigz docker.io curl gnupg

echo "Dependencies installed ✓"

# Install BWA-MEM2
echo "=== Installing BWA-MEM2 ==="
cd /tmp
wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2
cp bwa-mem2-2.2.1_x64-linux/bwa-mem2* /usr/local/bin/
bwa-mem2 version || { echo "ERROR: BWA-MEM2 installation failed"; exit 1; }
echo "BWA-MEM2 installed ✓"

# Install samtools
echo "=== Installing samtools ==="
wget -q https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2
tar -xf samtools-1.19.tar.bz2
cd samtools-1.19 && ./configure --with-libdeflate --prefix=/usr/local && make -j$(nproc) && make install
cd /tmp
samtools --version || { echo "ERROR: samtools installation failed"; exit 1; }
echo "samtools installed ✓"

# Install bcftools
echo "=== Installing bcftools ==="
wget -q https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar -xf bcftools-1.19.tar.bz2
cd bcftools-1.19 && ./configure --prefix=/usr/local && make -j$(nproc) && make install
cd /tmp
bcftools --version || { echo "ERROR: bcftools installation failed"; exit 1; }
echo "bcftools installed ✓"

# Setup Docker and NVIDIA toolkit
echo "=== Setting up Docker and NVIDIA Container Toolkit ==="
systemctl start docker
systemctl enable docker

# Download and install NVIDIA GPG key
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | gpg --batch --yes --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg || { echo "ERROR: Failed to download NVIDIA GPG key"; exit 1; }

# Add NVIDIA repository
curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
  sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
  tee /etc/apt/sources.list.d/nvidia-container-toolkit.list || { echo "ERROR: Failed to add NVIDIA repository"; exit 1; }

# Install NVIDIA Container Toolkit
apt-get update -qq
apt-get install -y -qq nvidia-container-toolkit || { echo "ERROR: Failed to install nvidia-container-toolkit"; exit 1; }
nvidia-ctk runtime configure --runtime=docker || { echo "ERROR: Failed to configure NVIDIA runtime"; exit 1; }
systemctl restart docker

echo "Docker and NVIDIA Container Toolkit installed ✓"

# Verify Docker GPU support (skip if no GPU present - this is an image builder VM)
echo "=== Checking Docker GPU support ==="
if nvidia-smi &>/dev/null; then
    echo "GPU detected, verifying Docker GPU support..."
    docker run --rm --gpus all nvidia/cuda:12.0.0-base-ubuntu22.04 nvidia-smi || { echo "ERROR: Docker GPU runtime not configured"; exit 1; }
    echo "Docker GPU support verified ✓"
else
    echo "No GPU detected (expected for image builder VM)"
    echo "NVIDIA Container Toolkit is installed and will work on GPU VMs ✓"
fi

# Pre-pull DeepVariant image
echo "=== Pulling DeepVariant GPU image ==="
docker pull google/deepvariant:1.6.0-gpu || { echo "ERROR: Failed to pull DeepVariant image"; exit 1; }
echo "DeepVariant image pulled ✓"

# Install GCP Ops Agent
echo "=== Installing GCP Ops Agent ==="
curl -sSO https://dl.google.com/cloudagents/add-google-cloud-ops-agent-repo.sh
bash add-google-cloud-ops-agent-repo.sh --also-install || { echo "ERROR: Failed to install Ops Agent"; exit 1; }
rm -f add-google-cloud-ops-agent-repo.sh
echo "GCP Ops Agent installed ✓"
echo "Note: Ops Agent will automatically start on VMs created from this image"

# Create reference directory
echo "=== Creating reference directory ==="
mkdir -p /opt/genomics/references/human
cd /opt/genomics/references/human

# Check available disk space
echo "=== Checking disk space ==="
df -h /opt/genomics
AVAILABLE_GB=$(df -BG /opt/genomics | tail -1 | awk '{print $4}' | sed 's/G//')
echo "Available space: ${AVAILABLE_GB} GB"

if [ "$AVAILABLE_GB" -lt 50 ]; then
    echo "ERROR: Insufficient disk space. Need at least 50 GB, have ${AVAILABLE_GB} GB"
    echo "BWA-MEM2 index generation requires significant temporary space"
    exit 1
fi

echo "Disk space check passed ✓"

# Download human reference genome
echo "=== Downloading human reference genome (GRCh38) ==="
echo "This will download ~3 GB of data (FASTA + all index files)..."
echo "Using wildcard to get all associated files (.fasta, .fai, .dict, etc.)"

# Download all reference files with wildcard (gets FASTA and all sidecar files)
gsutil -m cp "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta*" . || { echo "ERROR: Failed to download reference genome files"; exit 1; }
gsutil -m cp "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict" . || { echo "ERROR: Failed to download reference dict"; exit 1; }

echo "Reference genome downloaded ✓"
echo "Files downloaded:"
ls -lh

# Verify critical files exist and have reasonable sizes
echo ""
echo "=== Verifying reference genome integrity ==="

# Check FASTA file
if [ ! -f "Homo_sapiens_assembly38.fasta" ]; then
    echo "ERROR: Reference FASTA not found!"
    exit 1
fi

FASTA_SIZE=$(stat -c%s "Homo_sapiens_assembly38.fasta" 2>/dev/null || stat -f%z "Homo_sapiens_assembly38.fasta" 2>/dev/null)
FASTA_SIZE_GB=$((FASTA_SIZE / 1024 / 1024 / 1024))

echo "FASTA file size: ${FASTA_SIZE_GB} GB (${FASTA_SIZE} bytes)"

# Expected: ~3.1 GB (3100000000 bytes)
if [ "$FASTA_SIZE" -lt 3000000000 ]; then
    echo "ERROR: FASTA file too small (expected ~3.1 GB, got ${FASTA_SIZE_GB} GB)"
    echo "Download may be incomplete or corrupted"
    exit 1
fi

echo "FASTA size verification passed ✓"

# Check for required index files
if [ ! -f "Homo_sapiens_assembly38.fasta.fai" ]; then
    echo "ERROR: FASTA index (.fai) not found!"
    exit 1
fi

if [ ! -f "Homo_sapiens_assembly38.dict" ]; then
    echo "ERROR: Sequence dictionary (.dict) not found!"
    exit 1
fi

echo "All required reference files present ✓"

# Check disk space after download
echo ""
echo "Disk space after download:"
df -h /opt/genomics

# Setup swap file to prevent OOM during BWA-MEM2 indexing
echo "=== Setting up swap file to prevent OOM ==="
echo "Creating 20GB swap file (safety measure for BWA-MEM2 indexing)..."
fallocate -l 20G /swapfile || { echo "ERROR: Failed to create swap file"; exit 1; }
chmod 600 /swapfile
mkswap /swapfile || { echo "ERROR: Failed to format swap file"; exit 1; }
swapon /swapfile || { echo "ERROR: Failed to enable swap"; exit 1; }
echo "Swap enabled ✓"
free -h

# Generate BWA-MEM2 index
echo "========================================="
echo "=== Generating BWA-MEM2 index ==="
echo "Machine specs: $(nproc) CPUs, $(free -h | awk '/^Mem:/{print $2}') RAM + 20GB swap"
echo "This will take 30-60 minutes - please be patient!"
echo "Start time: $(date)"
echo "========================================="

# BWA-MEM2 index command
# Note: BWA-MEM2 doesn't have multi-threading for index generation
# The process is memory-intensive (requires ~48-64 GB RAM for human genome)
# Swap file prevents OOM killer if memory usage spikes above physical RAM
/usr/local/bin/bwa-mem2 index Homo_sapiens_assembly38.fasta || { echo "ERROR: BWA-MEM2 index generation failed"; exit 1; }

echo "========================================="
echo "BWA-MEM2 index generation completed!"
echo "End time: $(date)"
echo "========================================="

# Check disk space after indexing
echo "Disk space after indexing:"
df -h /opt/genomics

# Verify all index files were created
echo "=== Verifying index files ==="
ls -lh Homo_sapiens_assembly38.fasta.*

# Check for critical files
if [ ! -f "Homo_sapiens_assembly38.fasta.bwt.2bit.64" ]; then
    echo "ERROR: Critical index file .bwt.2bit.64 not found!"
    exit 1
fi

if [ ! -f "Homo_sapiens_assembly38.fasta.0123" ]; then
    echo "ERROR: Critical index file .0123 not found!"
    exit 1
fi

echo "All index files verified ✓"

# Set permissions
echo "=== Setting permissions ==="
chmod -R 755 /opt/genomics
echo "Permissions set ✓"

# Cleanup for image creation (reduces image size and improves stability)
echo "=== Cleaning up for image creation ==="

# Disable and remove swap file (not needed in the image)
echo "Removing swap file..."
swapoff /swapfile || echo "Warning: Failed to disable swap (may already be off)"
rm -f /swapfile
echo "Swap file removed ✓"

# Remove build-only dependencies to reduce image size
echo "Removing build tools (gcc, make, etc.) to reduce image size..."
apt-get remove -y -qq gcc g++ make cpp gcc-* dpkg-dev || echo "Warning: Some build tools may have already been removed"
apt-get autoremove -y -qq
echo "Build tools removed ✓"
echo "Estimated space saved: 200-400 MB"

# Remove temporary files from all build processes
echo "Removing temporary files..."
rm -rf /tmp/*
rm -rf /var/tmp/*

# Clean apt cache and package lists
echo "Cleaning apt cache and package lists..."
apt-get clean
rm -rf /var/lib/apt/lists/*
rm -rf /var/cache/apt/archives/*

# Clean Docker (but keep the DeepVariant image we pulled)
echo "Cleaning Docker system (keeping DeepVariant image)..."
docker system prune -f  # Remove dangling images, containers, networks (but not named images)

# Stop Docker cleanly to flush containerd state to disk
# This prevents stale container runtime state in the snapshot
echo "Stopping Docker service cleanly..."
systemctl stop docker.socket || echo "Warning: docker.socket already stopped"
systemctl stop docker || echo "Warning: docker service already stopped"
echo "Docker stopped ✓"

echo "Cleanup complete ✓"

# Upload log to GCS
echo "=== Uploading setup log to GCS ==="
gsutil cp $LOG_FILE gs://genomics-ref-bucket-binuv/logs/image-setup-$(date +%s).log || echo "Warning: Failed to upload log"

# Summary
echo ""
echo "========================================="
echo "===   SETUP COMPLETE!                 ==="
echo "========================================="
echo "Completion time: $(date)"
echo ""
echo "Installed tools:"
echo "  - BWA-MEM2 $(bwa-mem2 version 2>&1 | head -1)"
echo "  - samtools $(samtools --version | head -1)"
echo "  - bcftools $(bcftools --version | head -1)"
echo "  - Docker with NVIDIA GPU support"
echo "  - DeepVariant 1.6.0-gpu (pre-pulled)"
echo "  - GCP Ops Agent (for Cloud Monitoring/Logging)"
echo ""
echo "Reference genome:"
echo "  - Location: /opt/genomics/references/human/"
echo "  - GRCh38 (Homo_sapiens_assembly38.fasta)"
echo "  - BWA-MEM2 index: COMPLETE"
echo ""
echo "Cleanup performed:"
echo "  - Swap file removed"
echo "  - Build tools removed (gcc, make, etc.) - saved ~200-400 MB"
echo "  - Temporary files cleaned"
echo "  - Apt cache purged"
echo "  - Docker stopped cleanly (ready for image creation)"
echo ""
echo "Next steps:"
echo ""
echo "1. Verify your project ID (IMPORTANT!):"
echo "   gcloud config get-value project"
echo "   Expected: genomics-pipeline-poc"
echo ""
echo "2. Stop this VM:"
echo "   gcloud compute instances stop genomics-image-builder --zone=us-central1-a"
echo ""
echo "3. Create image with high-performance features:"
echo "   gcloud compute images create genomics-pipeline-human-ready-v1 \\"
echo "     --source-disk=genomics-image-builder \\"
echo "     --source-disk-zone=us-central1-a \\"
echo "     --family=genomics-pipeline-human-ready \\"
echo "     --project=genomics-pipeline-poc \\"
echo "     --guest-os-features=GVNIC,VIRTIO_SCSI_MULTIQUEUE"
echo ""
echo "   Guest OS features enabled:"
echo "   - GVNIC: Google Virtual NIC for high-bandwidth networking (up to 100 Gbps)"
echo "   - VIRTIO_SCSI_MULTIQUEUE: Multi-queue SCSI for faster disk I/O"
echo ""
echo "   IMPORTANT: Ensure --project matches your actual project ID!"
echo ""
echo "4. Delete this VM:"
echo "   gcloud compute instances delete genomics-image-builder --zone=us-central1-a --quiet"
echo "========================================="
echo ""
echo "VM is ready for image creation!"
echo "DO NOT delete this VM yet - create the image first!"
echo ""
