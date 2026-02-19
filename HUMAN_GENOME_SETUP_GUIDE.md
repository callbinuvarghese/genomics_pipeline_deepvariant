# Human Genome Pipeline Setup Guide

## Overview
This guide walks you through setting up a production-ready human genomics pipeline using GPU-accelerated DeepVariant on GCP.

---

## Part 1: Create Custom VM Image (One-Time Setup)

### Why Create a Custom Image?
Your current VM has all tools pre-installed (BWA-MEM2, samtools, bcftools, Docker, NVIDIA toolkit, DeepVariant image). Creating a custom image from this saves **~3-5 minutes** per pipeline run by skipping installation steps.

### Steps:

#### 1. Stop the current VM (optional but recommended)
```bash
gcloud compute instances stop genomics-pipeline-vm-gpu --zone=us-central1-a
```

#### 2. Create the custom image
```bash
gcloud compute images create genomics-pipeline-base-gpu \
  --source-disk=genomics-pipeline-vm-gpu \
  --source-disk-zone=us-central1-a \
  --family=genomics-pipeline-gpu \
  --description="Pre-configured genomics pipeline: Deep Learning VM (CUDA 12.8, NVIDIA 570) + BWA-MEM2 2.2.1 + samtools 1.19 + bcftools 1.19 + Docker + NVIDIA Container Toolkit + DeepVariant 1.6.0-gpu"
```

This creates an image in the `genomics-pipeline-gpu` family that you can reference in future VM creations.

#### 3. (Optional) Delete the original VM to save costs
```bash
gcloud compute instances delete genomics-pipeline-vm-gpu --zone=us-central1-a --quiet
```

#### 4. Verify the image was created
```bash
gcloud compute images describe genomics-pipeline-base-gpu
```

---

## Part 2: Prepare Human Reference Genome Data

### Required Files for GRCh38/hg38:

You need these files in `gs://genomics-ref-bucket-binuv/human/`:

1. **Reference FASTA**: `Homo_sapiens_assembly38.fasta` (~3 GB)
2. **FASTA Index**: `Homo_sapiens_assembly38.fasta.fai`
3. **Sequence Dictionary**: `Homo_sapiens_assembly38.dict`
4. **BWA-MEM2 Index** (recommended to pre-generate):
   - `Homo_sapiens_assembly38.fasta.0123`
   - `Homo_sapiens_assembly38.fasta.amb`
   - `Homo_sapiens_assembly38.fasta.ann`
   - `Homo_sapiens_assembly38.fasta.bwt.2bit.64`
   - `Homo_sapiens_assembly38.fasta.pac`

### Where to Download Reference Files:

#### Option 1: Broad Institute (Recommended)
```bash
# Download from Broad's public bucket
gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta gs://genomics-ref-bucket-binuv/human/
gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai gs://genomics-ref-bucket-binuv/human/
gsutil -m cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict gs://genomics-ref-bucket-binuv/human/
```

#### Option 2: NCBI (Alternative)
Download from: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/

### Pre-Generate BWA-MEM2 Index (Highly Recommended)

**IMPORTANT**: Generating BWA-MEM2 index for human genome takes **30-60 minutes**. Do this once and reuse:

```bash
# Option A: Generate locally if you have a powerful machine
bwa-mem2 index Homo_sapiens_assembly38.fasta
gsutil -m cp Homo_sapiens_assembly38.fasta.* gs://genomics-ref-bucket-binuv/human/

# Option B: Use a temporary GCP VM to generate
gcloud compute instances create bwa-index-builder \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --image-family=genomics-pipeline-gpu \
  --image-project=YOUR_PROJECT_ID \
  --boot-disk-size=200GB \
  --scopes=cloud-platform \
  --metadata=startup-script='#!/bin/bash
    cd /tmp
    gsutil cp gs://genomics-ref-bucket-binuv/human/Homo_sapiens_assembly38.fasta .
    bwa-mem2 index Homo_sapiens_assembly38.fasta
    gsutil -m cp Homo_sapiens_assembly38.fasta.* gs://genomics-ref-bucket-binuv/human/
    gcloud compute instances delete $(hostname) --zone=us-central1-a --quiet
  '
```

---

## Part 3: Prepare Input FASTQ Files

Upload your human FASTQ files to the input bucket:

```bash
# Example: Upload paired-end FASTQ files
gsutil -m cp NA12878_R1.fastq.gz gs://genomics-input-bucket-binuv/human/
gsutil -m cp NA12878_R2.fastq.gz gs://genomics-input-bucket-binuv/human/
```

**Naming convention**: The script expects `${SAMPLE_ID}_R1.fastq.gz` and `${SAMPLE_ID}_R2.fastq.gz`

---

## Part 4: Upload the Human Pipeline Script

```bash
# Upload the new human genome script
gsutil cp genomics_pipeline_human_gpu.sh gs://genomics-ref-bucket-binuv/scripts/
```

---

## Part 5: Run the Human Genome Pipeline

### Using the Custom Image (Recommended)

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-gpu \
  --image-project=YOUR_PROJECT_ID \
  --boot-disk-size=500GB \
  --boot-disk-type=pd-ssd \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu.sh,SAMPLE_ID=NA12878
```

**Key changes from E. coli pipeline:**
- `--machine-type=n1-standard-16` (more CPU cores for human genome)
- `--boot-disk-size=500GB` (human genome data is much larger)
- `--image-family=genomics-pipeline-gpu` (uses your custom image)
- `--metadata=...,SAMPLE_ID=NA12878` (sets the sample ID)

### Using the Original Deep Learning Image (No Custom Image)

If you don't create the custom image, use the original Deep Learning VM image and the E. coli script will install everything:

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=common-cu128-ubuntu-2204-nvidia-570 \
  --image-project=deeplearning-platform-release \
  --boot-disk-size=500GB \
  --boot-disk-type=pd-ssd \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu.sh,SAMPLE_ID=NA12878
```

---

## Part 6: Monitor the Pipeline

### Watch startup script logs
```bash
gcloud compute ssh genomics-pipeline-human-vm --zone=us-central1-a \
  --command="sudo journalctl -u google-startup-scripts.service -f"
```

### SSH and tail the pipeline log
```bash
gcloud compute ssh genomics-pipeline-human-vm --zone=us-central1-a
sudo tail -f /mnt/scratch/pipeline_human_gpu.log
```

---

## Part 7: Review Results

After the pipeline completes (VM auto-deletes), check the outputs:

```bash
# List all outputs
gsutil ls gs://genomics-output-bucket-binuv/human/

# Download specific results
gsutil cp gs://genomics-output-bucket-binuv/human/NA12878.vcf.gz .
gsutil cp gs://genomics-output-bucket-binuv/human/NA12878.flagstat.txt .

# View VCF stats
gsutil cat gs://genomics-output-bucket-binuv/human/NA12878.vcf.stats.txt | head -50

# Check logs
gsutil ls gs://genomics-ref-bucket-binuv/logs/
```

---

## Key Differences: E. coli vs Human Genome

| Parameter | E. coli | Human |
|-----------|---------|-------|
| **Reference size** | 4.6 MB | ~3 GB |
| **FASTQ size** | ~500 MB | ~50-100 GB |
| **Machine type** | n1-standard-4 | n1-standard-16 |
| **Boot disk** | 100 GB | 500 GB |
| **Alignment time** | ~78s | ~2-4 hours |
| **Total runtime** | ~8 minutes | ~4-8 hours |
| **Cost per run (preemptible)** | ~$0.10 | ~$2-5 |

---

## Cost Optimization Tips

### 1. Use Preemptible VMs
Already configured with `--preemptible`. Saves 60-80% vs on-demand.

### 2. Use Custom Image
Saves 3-5 minutes of installation time = cost savings.

### 3. Right-size the Machine
- **30x coverage WGS**: n1-standard-16 or n1-standard-32
- **Exome sequencing**: n1-standard-8
- **Targeted panels**: n1-standard-4

### 4. Pre-generate BWA-MEM2 Index
Saves 30-60 minutes per run.

### 5. Consider Spot VMs
GCP's newer alternative to preemptible VMs with better availability.

---

## Troubleshooting

### Issue: "Out of disk space"
**Solution**: Increase `--boot-disk-size` or modify script to delete intermediate files more aggressively.

### Issue: "VM was preempted"
**Solution**: Implement checkpointing or use non-preemptible VMs for critical runs.

### Issue: "BWA-MEM2 index generation takes too long"
**Solution**: Pre-generate the index once and upload to GCS (see Part 2).

### Issue: "DeepVariant running slowly"
**Solution**:
- Verify GPU is being used: check logs for "Using GPU" messages
- Try larger GPU: `--accelerator=type=nvidia-tesla-v100,count=1`

---

## Next Steps

1. **Create the custom image** (saves time and money)
2. **Download and upload human reference genome** to GCS
3. **Pre-generate BWA-MEM2 index** (saves 30-60 min per run)
4. **Upload sample FASTQ files** to input bucket
5. **Run the pipeline** with the command in Part 5
6. **Review results** and adjust machine type/disk size as needed

---

## Sample Test Data

For testing, use Genome in a Bottle (GIAB) sample NA12878:
- Available in SRA: `SRR1596637` (downsampled to 30x coverage)
- Reference truth set available for validation
- Widely used benchmark for variant calling

```bash
# Download GIAB sample (example)
# This is a large dataset - consider using a smaller subset for testing
# See: https://github.com/genome-in-a-bottle/giab_data_indexes
```

---

## Additional Resources

- **DeepVariant Documentation**: https://github.com/google/deepvariant
- **BWA-MEM2 Paper**: https://ieeexplore.ieee.org/document/8820962
- **GATK Best Practices**: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932
- **GCP Genomics**: https://cloud.google.com/life-sciences/docs/resources/public-datasets
