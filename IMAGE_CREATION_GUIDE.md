# Custom Image Creation Guide

## Overview

This guide covers creating the custom genomics pipeline image with high-performance features enabled.

## Prerequisites

1. **Builder VM must be stopped** (not deleted!)
2. **All setup steps completed** (check the setup log)
3. **Sufficient quota** for custom images in your project

## Image Creation Command

### Standard Command (Recommended)

```bash
gcloud compute images create genomics-pipeline-human-ready-v1 \
  --source-disk=genomics-image-builder \
  --source-disk-zone=us-central1-a \
  --family=genomics-pipeline-human-ready \
  --project=genomics-pipeline-poc \
  --guest-os-features=GVNIC,VIRTIO_SCSI_MULTIQUEUE \
  --description="Human genomics pipeline image: Ubuntu 22.04 + CUDA 12.8 + NVIDIA 570 + BWA-MEM2 2.2.1 + samtools 1.19 + bcftools 1.19 + DeepVariant 1.6.0-gpu + GCP Ops Agent + Pre-loaded GRCh38 reference genome with BWA-MEM2 index. Features: GVNIC for high-bandwidth networking, VIRTIO_SCSI_MULTIQUEUE for faster disk I/O."
```

## Guest OS Features Explained

### GVNIC (Google Virtual NIC)

**What it is:**
- Google's next-generation virtual network interface
- Enables high-bandwidth networking for supported machine types

**Benefits:**
- ✅ **Up to 100 Gbps bandwidth** (vs 32 Gbps with virtIO)
- ✅ **Lower latency** for network operations
- ✅ **Faster FASTQ downloads** from GCS buckets
- ✅ **Faster result uploads** to GCS

**When it matters:**
- Downloading large FASTQ files (30x WGS = ~75 GB)
- Uploading large BAM files (~100-200 GB)
- Using high-bandwidth machine types (n2, n2d, c2, c2d)

**Performance impact:**
- 75 GB FASTQ download: ~10 seconds (100 Gbps) vs ~30 seconds (32 Gbps)
- **3x faster network I/O** on supported machines

### VIRTIO_SCSI_MULTIQUEUE

**What it is:**
- Multi-queue block device driver
- Enables parallel disk I/O operations across multiple CPU cores

**Benefits:**
- ✅ **Higher IOPS** (Input/Output Operations Per Second)
- ✅ **Better disk throughput** on multi-core machines
- ✅ **Scales with CPU count** (more cores = more queues)
- ✅ **Reduced I/O latency**

**When it matters:**
- Writing large BAM files (~100-200 GB)
- Sorting operations (high random I/O)
- Using machines with many vCPUs (16+)

**Performance impact:**
- 16 vCPUs: ~60% more IOPS with multiqueue
- Better utilization of SSD performance

## Machine Type Compatibility

### GVNIC Support

| Machine Family | GVNIC Support | Max Bandwidth |
|----------------|---------------|---------------|
| n1 | ❌ No | 32 Gbps (virtIO) |
| n2, n2d | ✅ Yes | 50-100 Gbps |
| c2, c2d | ✅ Yes | 100 Gbps |
| m1, m2, m3 | ✅ Yes | 100 Gbps |
| a2 (GPU) | ✅ Yes | 100 Gbps |

**For n1-standard-16:**
- GVNIC feature is ignored (falls back to virtIO)
- Still safe to include in image
- Future-proofs for n2 machine upgrades

### VIRTIO_SCSI_MULTIQUEUE Support

- ✅ **All machine types** support this feature
- Automatically scales queue count with vCPU count
- No downside to enabling

## Image Creation Timeline

| Step | Duration | Notes |
|------|----------|-------|
| Stop VM | ~30 sec | Wait for clean shutdown |
| Create image | ~5-10 min | Snapshots the 200GB disk |
| Image ready | - | Can now launch VMs |

## Verification

After image creation, verify it exists:

```bash
# List images in the family
gcloud compute images list --filter="family=genomics-pipeline-human-ready" --project=genomics-pipeline-poc

# Describe the image
gcloud compute images describe genomics-pipeline-human-ready-v1 --project=genomics-pipeline-poc

# Check guest OS features
gcloud compute images describe genomics-pipeline-human-ready-v1 \
  --project=genomics-pipeline-poc \
  --format="value(guestOsFeatures)"
```

**Expected output:**
```
name: GVNIC
name: VIRTIO_SCSI_MULTIQUEUE
```

## Using the Image

### VM Creation with High-Performance Networking

#### Option 1: n1-standard-16 (Current)

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --project=genomics-pipeline-poc \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-human-ready \
  --image-project=genomics-pipeline-poc \
  --boot-disk-size=100GB \
  --boot-disk-type=pd-ssd \
  --local-ssd=interface=NVME \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=NA12878
```

**Note:** n1 machines will use virtIO networking (32 Gbps), but VIRTIO_SCSI_MULTIQUEUE is active.

#### Option 2: n2-standard-16 (High-Bandwidth Networking)

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --project=genomics-pipeline-poc \
  --zone=us-central1-a \
  --machine-type=n2-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-human-ready \
  --image-project=genomics-pipeline-poc \
  --boot-disk-size=100GB \
  --boot-disk-type=pd-ssd \
  --local-ssd=interface=NVME \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --network-tier=PREMIUM \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=NA12878
```

**Benefits:**
- ✅ GVNIC active (up to 100 Gbps)
- ✅ VIRTIO_SCSI_MULTIQUEUE active
- ✅ Newer generation CPU (Cascade Lake or Ice Lake)

**Cost:** ~10% more expensive than n1

## Performance Comparison

### FASTQ Download (75 GB)

| Setup | Bandwidth | Download Time |
|-------|-----------|---------------|
| n1 + virtIO | 32 Gbps | ~30 seconds |
| n2 + GVNIC | 100 Gbps | ~10 seconds |
| **Savings** | - | **20 seconds** |

### Disk I/O (100 GB BAM write)

| Setup | IOPS | Write Time |
|-------|------|------------|
| Single queue | ~30K IOPS | ~200 seconds |
| Multiqueue (16 cores) | ~50K IOPS | ~120 seconds |
| **Savings** | - | **80 seconds** |

### Total Pipeline Impact

For 30x WGS with 75 GB FASTQ:
- **Network savings:** ~20 seconds (download) + ~15 seconds (upload) = **35 seconds**
- **Disk I/O savings:** ~80 seconds (sorting/writing BAM)
- **Total savings:** ~2 minutes per genome

**At scale (100 genomes):**
- Time saved: ~200 minutes (3.3 hours)
- Cost saved at $0.76/hr: ~$2.50

## Cost Considerations

### Image Storage

- **Size:** ~29 GB
- **Cost:** ~$1.45/month ($0.05/GB/month)
- **One-time creation cost:** ~$0.76 (n1-standard-16 × 60 min)

### Network Tier

To use GVNIC effectively:
- Use `--network-tier=PREMIUM` (default)
- Avoid `--network-tier=STANDARD` (limits bandwidth)

## GCP Ops Agent Identity Management

### The Issue

When you create an image with Ops Agent pre-installed, the agent has a unique ID tied to the builder VM. If not handled correctly, new VMs might report metrics using the builder VM's identity for the first few minutes.

### The Solution

**Already implemented in pipeline scripts:**

Both pipeline scripts automatically restart the Ops Agent on VM boot to generate a fresh identity:

```bash
# From genomics_pipeline_human_gpu_preloaded_checkpointed.sh
if systemctl is-active --quiet google-cloud-ops-agent; then
  echo "Ops Agent detected (pre-installed in image)"
  echo "Restarting to generate fresh VM identity..."
  sudo systemctl restart google-cloud-ops-agent
  echo "Ops Agent restarted with new VM identity ✓"
fi
```

**What this does:**
- ✅ Detects pre-installed Ops Agent
- ✅ Restarts the service on first boot
- ✅ Generates unique VM-specific identity
- ✅ Ensures metrics are correctly attributed to the new VM

**Verification:**

After VM starts, check that it has its own identity:

```bash
# Check Ops Agent status
sudo systemctl status google-cloud-ops-agent

# View metrics in Cloud Console - should show metrics under the new VM name
# NOT under "genomics-image-builder"
```

## Troubleshooting

### Issue: Guest OS features not appearing

**Symptom:** `gcloud compute images describe` doesn't show GVNIC/VIRTIO_SCSI_MULTIQUEUE

**Solution:**
```bash
# Delete and recreate image with features
gcloud compute images delete genomics-pipeline-human-ready-v1 --project=genomics-pipeline-poc
gcloud compute images create genomics-pipeline-human-ready-v1 \
  --source-disk=genomics-image-builder \
  --source-disk-zone=us-central1-a \
  --family=genomics-pipeline-human-ready \
  --project=genomics-pipeline-poc \
  --guest-os-features=GVNIC,VIRTIO_SCSI_MULTIQUEUE
```

### Issue: VM doesn't use GVNIC

**Symptom:** Network bandwidth limited to 32 Gbps on n2 machine

**Possible causes:**
1. Image doesn't have GVNIC feature
2. Machine type doesn't support GVNIC (e.g., n1)
3. Using STANDARD network tier

**Solution:**
1. Verify image has GVNIC: `gcloud compute images describe ...`
2. Use supported machine type: n2, c2, a2
3. Use PREMIUM network tier (default)

## Best Practices

1. ✅ **Always include both features** - No downside, future-proofs image
2. ✅ **Use descriptive image names** - Include version number
3. ✅ **Add detailed description** - Document what's in the image
4. ✅ **Use image families** - Easier to reference latest version
5. ✅ **Test the image** - Launch a test VM before production use

## Next Steps

After creating the image:

1. **Test VM launch:**
   ```bash
   gcloud compute instances create test-vm \
     --image-family=genomics-pipeline-human-ready \
     --image-project=genomics-pipeline-poc \
     --zone=us-central1-a
   ```

2. **Verify features:**
   ```bash
   gcloud compute instances describe test-vm --zone=us-central1-a \
     --format="value(networkInterfaces[0].nicType)"
   ```
   Expected: `GVNIC` (on n2) or `VIRTIO_NET` (on n1)

3. **Delete test VM:**
   ```bash
   gcloud compute instances delete test-vm --zone=us-central1-a --quiet
   ```

4. **Run production pipeline!**

---

## Related Documentation

- [setup_preloaded_image.sh](setup_preloaded_image.sh) - Image setup script
- [HUMAN_GENOME_SETUP_GUIDE.md](HUMAN_GENOME_SETUP_GUIDE.md) - Complete setup guide
- [IMPROVEMENTS_SUMMARY.md](IMPROVEMENTS_SUMMARY.md) - All pipeline improvements
