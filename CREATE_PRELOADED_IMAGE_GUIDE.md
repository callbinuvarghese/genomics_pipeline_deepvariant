# Creating a Pre-loaded Image with Reference Genome

## Overview
This guide shows how to create a VM image with **everything pre-loaded**: tools, Docker, NVIDIA toolkit, reference genome, and BWA-MEM2 index. This eliminates the need to download the reference genome on every run, saving significant time and bandwidth.

---

## Time Savings Comparison

| Approach | Reference Download | BWA-MEM2 Indexing | Tool Installation | Total Setup Time |
|----------|-------------------|-------------------|-------------------|------------------|
| **No Image** | ~5-10 min | ~30-60 min | ~3-5 min | **38-75 min** |
| **Tools Only** | ~5-10 min | ~30-60 min | 0 min | **35-70 min** |
| **Full Pre-load** | 0 min | 0 min | 0 min | **~10 sec** |

**Savings: 35-75 minutes per pipeline run!**

---

## Step-by-Step Instructions

### Step 1: Start Your Current VM
```bash
gcloud compute instances start genomics-pipeline-vm-gpu --zone=us-central1-a
```

### Step 2: SSH into the VM
```bash
gcloud compute ssh genomics-pipeline-vm-gpu --zone=us-central1-a
```

### Step 3: Download and Prepare Reference Genome

```bash
# Create reference directory
sudo mkdir -p /opt/genomics/references/human
cd /opt/genomics/references/human

# Download GRCh38 reference genome from Broad Institute
echo "Downloading reference genome (~3 GB)..."
sudo gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta .
sudo gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai .
sudo gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict .

echo "Download complete. Verifying files..."
ls -lh

# Verify files
if [ ! -f "Homo_sapiens_assembly38.fasta" ]; then
  echo "ERROR: Reference FASTA not downloaded"
  exit 1
fi

echo "Reference genome downloaded successfully ✓"
```

### Step 4: Generate BWA-MEM2 Index (30-60 minutes)

**IMPORTANT**: This is the most time-consuming step. Make sure the VM doesn't disconnect.

```bash
# Generate BWA-MEM2 index
echo "Generating BWA-MEM2 index - this will take 30-60 minutes..."
echo "Start time: $(date)"

sudo /usr/local/bin/bwa-mem2 index Homo_sapiens_assembly38.fasta

echo "End time: $(date)"
echo "BWA-MEM2 index generation complete ✓"

# Verify all index files were created
echo "Verifying index files..."
ls -lh Homo_sapiens_assembly38.fasta.*

# Expected files:
# Homo_sapiens_assembly38.fasta.0123
# Homo_sapiens_assembly38.fasta.amb
# Homo_sapiens_assembly38.fasta.ann
# Homo_sapiens_assembly38.fasta.bwt.2bit.64
# Homo_sapiens_assembly38.fasta.pac
```

### Step 5: Set Proper Permissions

```bash
# Set directory permissions
sudo chmod -R 755 /opt/genomics

# Verify permissions
ls -la /opt/genomics/references/human/

echo "Permissions set correctly ✓"
```

### Step 6: Exit SSH and Stop the VM

```bash
# Exit SSH session
exit

# Stop the VM (from your local terminal)
gcloud compute instances stop genomics-pipeline-vm-gpu --zone=us-central1-a
```

Wait for the VM to fully stop before proceeding.

### Step 7: Create the Enhanced Custom Image

```bash
gcloud compute images create genomics-pipeline-human-ready \
  --source-disk=genomics-pipeline-vm-gpu \
  --source-disk-zone=us-central1-a \
  --family=genomics-pipeline-gpu \
  --description="Complete genomics pipeline with GRCh38 reference and BWA-MEM2 index pre-loaded: Deep Learning VM (CUDA 12.8, NVIDIA 570) + BWA-MEM2 2.2.1 + samtools 1.19 + bcftools 1.19 + Docker + NVIDIA Container Toolkit + DeepVariant 1.6.0-gpu + GRCh38 reference genome + BWA-MEM2 index"
```

**Image creation takes 5-15 minutes**. You'll see output like:
```
Created [https://www.googleapis.com/compute/v1/projects/genomics-pipeline-poc/global/images/genomics-pipeline-human-ready].
```

### Step 8: Verify the Image

```bash
# Check image details
gcloud compute images describe genomics-pipeline-human-ready

# List all images in the family
gcloud compute images list --filter="family:genomics-pipeline-gpu"
```

### Step 9: Upload the Optimized Pipeline Script

```bash
# Upload the new pre-loaded pipeline script
gsutil cp genomics_pipeline_human_gpu_preloaded.sh gs://genomics-ref-bucket-binuv/scripts/
```

### Step 10: (Optional) Delete the Original VM

Since you now have the image, you can delete the VM to save costs:

```bash
gcloud compute instances delete genomics-pipeline-vm-gpu --zone=us-central1-a --quiet
```

You can always recreate VMs from the image.

---

## Running the Pipeline with Pre-loaded Image

### Quick Run Command

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-gpu \
  --image-project=genomics-pipeline-poc \
  --boot-disk-size=500GB \
  --boot-disk-type=pd-ssd \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded.sh,SAMPLE_ID=NA12878
```

**Important**: Replace `--image-project=genomics-pipeline-poc` with your actual GCP project ID.

### What Happens at Runtime

1. ✅ **VM boots** (~30 seconds) - All tools, Docker, reference genome already present
2. ✅ **Verifies reference** (~1 second) - Confirms files at `/opt/genomics/references/human/`
3. ✅ **Downloads FASTQ** (~5-15 minutes depending on file size)
4. ✅ **Runs alignment** (BWA-MEM2 uses pre-generated index)
5. ✅ **Variant calling** (DeepVariant with GPU)
6. ✅ **Uploads results** to GCS
7. ✅ **Auto-deletes VM**

**Total time for 30x WGS: ~4-6 hours** (vs 5-7 hours without pre-loaded image)

---

## What's Inside the Image

### 1. Pre-installed Tools
- BWA-MEM2 v2.2.1
- samtools v1.19
- bcftools v1.19
- Docker with NVIDIA Container Toolkit
- DeepVariant 1.6.0-gpu (Docker image pre-pulled)

### 2. Reference Genome (in `/opt/genomics/references/human/`)
- `Homo_sapiens_assembly38.fasta` (3.0 GB)
- `Homo_sapiens_assembly38.fasta.fai` (14 KB)
- `Homo_sapiens_assembly38.dict` (250 KB)

### 3. BWA-MEM2 Index Files (in `/opt/genomics/references/human/`)
- `Homo_sapiens_assembly38.fasta.0123` (~800 MB)
- `Homo_sapiens_assembly38.fasta.amb` (~6 KB)
- `Homo_sapiens_assembly38.fasta.ann` (~6 KB)
- `Homo_sapiens_assembly38.fasta.bwt.2bit.64` (~800 MB)
- `Homo_sapiens_assembly38.fasta.pac` (~800 MB)

**Total image size**: ~10-12 GB (vs ~5 GB for base Deep Learning VM)

---

## Cost Analysis

### Image Storage Costs
- **Storage**: $0.05/GB/month
- **Image size**: ~10-12 GB
- **Monthly cost**: ~$0.50-0.60/month

### Savings Per Run
- **Bandwidth saved**: ~3 GB reference download = $0.12-0.23
- **Compute saved**: 30-60 min indexing on n1-standard-16 = $0.40-0.80
- **Total savings per run**: ~$0.52-1.03

**Break-even point**: After 1 run, the image pays for itself!

---

## Alternative: Add Multiple Reference Genomes

You can extend this approach to include multiple reference genomes in one image:

```bash
# In the VM before creating the image
sudo mkdir -p /opt/genomics/references/mouse
cd /opt/genomics/references/mouse

# Download mouse reference (GRCm39)
sudo gsutil -m cp gs://gcp-public-data--broad-references/mm10/v0/Mus_musculus_assembly38.fasta .
# ... and generate BWA-MEM2 index

# Then create image with both human and mouse references
```

---

## Updating the Image

When tools or references need updates:

1. **Create a new VM from the current image**
   ```bash
   gcloud compute instances create update-genomics-image \
     --zone=us-central1-a \
     --machine-type=n1-standard-8 \
     --image-family=genomics-pipeline-gpu \
     --image-project=genomics-pipeline-poc \
     --boot-disk-size=100GB
   ```

2. **Make updates** (SSH in and update tools/references)

3. **Create a new image**
   ```bash
   gcloud compute images create genomics-pipeline-human-ready-v2 \
     --source-disk=update-genomics-image \
     --source-disk-zone=us-central1-a \
     --family=genomics-pipeline-gpu
   ```

4. **(Optional) Delete old image**
   ```bash
   gcloud compute images delete genomics-pipeline-human-ready
   ```

Using `--family=genomics-pipeline-gpu` ensures VMs always use the latest image in the family.

---

## Troubleshooting

### Issue: "Reference genome not found at /opt/genomics/references/human/"

**Cause**: VM was created from wrong image or image creation failed.

**Solution**:
```bash
# Check which image the VM is using
gcloud compute instances describe YOUR-VM-NAME --zone=us-central1-a | grep sourceImage

# Verify reference files in the image
gcloud compute ssh YOUR-VM-NAME --zone=us-central1-a --command="ls -la /opt/genomics/references/human/"
```

### Issue: "BWA-MEM2 index not found"

**Cause**: Index generation step was skipped or failed.

**Solution**: Verify all index files exist in the image:
```bash
gcloud compute ssh YOUR-VM-NAME --zone=us-central1-a --command="ls -lh /opt/genomics/references/human/Homo_sapiens_assembly38.fasta.*"
```

### Issue: Image creation is slow

**Normal**: Creating images from large disks (100GB+) takes 10-20 minutes. Be patient.

---

## Best Practices

1. ✅ **Version your images**: Use descriptive names with dates or versions
   - `genomics-pipeline-human-ready-2026-02-12`
   - `genomics-pipeline-human-ready-v1.0`

2. ✅ **Document image contents**: Use detailed `--description` flags

3. ✅ **Use image families**: Allows automatic updates without changing VM creation scripts

4. ✅ **Test before production**: Create a test VM from the image and verify all files

5. ✅ **Keep images up to date**: Update quarterly or when major tool versions change

---

## Summary

By creating a pre-loaded image, you've optimized your genomics pipeline to:

- ✅ **Save 35-75 minutes** per run (reference download + indexing)
- ✅ **Reduce bandwidth costs** (~$0.12-0.23 per run)
- ✅ **Simplify pipeline script** (no complex download logic)
- ✅ **Improve reliability** (no network dependency for reference files)
- ✅ **Enable offline testing** (reference always available)

**Total monthly cost**: ~$0.50-0.60 for image storage
**Savings per run**: ~$0.52-1.03
**Break-even**: After first run!

🎉 **Your optimized genomics pipeline is ready!**
