# Local SSD Setup for Genomics Pipeline

## Why Use Local SSD?

### Performance Comparison

| Metric | Standard PD | SSD PD | **Local NVMe SSD** |
|--------|-------------|--------|-------------------|
| **Read IOPS** | 0.75-7.5K | 15K-60K | **170K-660K** |
| **Write IOPS** | 1.5K | 15K-30K | **120K-240K** |
| **Throughput** | 240 MB/s | 1,200 MB/s | **2,400 MB/s** |
| **Latency** | ~10ms | ~1ms | **<1ms** |
| **Cost/GB-month** | $0.04 | $0.17 | **$0.08** |

### Pipeline Performance Impact

For human WGS (30x coverage):

| Step | Standard PD | Local SSD | **Speedup** |
|------|-------------|-----------|-------------|
| **Alignment + Name Sort** | ~3.5 hours | **~2 hours** | **1.75x** |
| **Fixmate + Sort** | ~45 min | **~20 min** | **2.25x** |
| **Mark Duplicates** | ~30 min | **~15 min** | **2x** |
| **DeepVariant** | ~2 hours | **~1.5 hours** | **1.33x** |
| **TOTAL** | **~6.5 hours** | **~4 hours** | **1.6x faster** |

**Cost savings from faster completion often offset the slightly higher Local SSD cost!**

## Local SSD Sizing

### Size Options

Local SSDs come in **375 GB increments**:
- 1x Local SSD = 375 GB
- 2x Local SSD = 750 GB
- 3x Local SSD = 1,125 GB
- etc.

### Recommended Sizes by Coverage

| Use Case | FASTQ Size | Peak Disk Usage | **Recommended** |
|----------|------------|-----------------|-----------------|
| **Test data (4x)** | ~10 GB | ~40 GB | **1x (375 GB)** |
| **Exome (30x)** | ~15 GB | ~60 GB | **1x (375 GB)** |
| **WGS (30x)** | ~75 GB | ~300 GB | **1x (375 GB)** |
| **WGS (60x)** | ~150 GB | ~600 GB | **2x (750 GB)** |
| **WGS (100x)** | ~250 GB | ~1,000 GB | **3x (1,125 GB)** |

**Rule of thumb**: Peak disk usage = FASTQ size × 4

## VM Creation with Local SSD

### Option 1: Single Local SSD (375 GB) - Recommended for Most Workflows

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-human-ready \
  --image-project=YOUR_PROJECT_ID \
  --boot-disk-size=100GB \
  --boot-disk-type=pd-ssd \
  --local-ssd interface=NVME \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=NA12878
```

**Key changes from previous command:**
- ✅ `--local-ssd interface=NVME` - Adds 375 GB Local NVMe SSD
- ✅ `--boot-disk-size=100GB` - Reduced (we now use Local SSD for scratch)
- ✅ Script auto-detects and mounts Local SSD

### Option 2: Multiple Local SSDs (750 GB+) - For High Coverage

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-human-ready \
  --image-project=YOUR_PROJECT_ID \
  --boot-disk-size=100GB \
  --boot-disk-type=pd-ssd \
  --local-ssd interface=NVME \
  --local-ssd interface=NVME \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=NA12878
```

**Note**: Add `--local-ssd interface=NVME` multiple times for additional 375 GB increments.

If using multiple Local SSDs, you'll need to set up RAID for better performance:

```bash
# In your startup script, before running the pipeline
sudo mdadm --create /dev/md0 --level=0 --raid-devices=2 /dev/nvme0n1 /dev/nvme0n2
sudo mkfs.ext4 -F /dev/md0
sudo mkdir -p /mnt/disks/scratch
sudo mount /dev/md0 /mnt/disks/scratch
sudo chmod 777 /mnt/disks/scratch
```

## Auto-Detection in Pipeline Script

The updated pipeline script automatically:

1. **Detects Local SSD**: Checks for `/dev/nvme0n1`
2. **Formats and mounts**: Creates ext4 filesystem at `/mnt/disks/scratch`
3. **Falls back to boot disk**: Uses `/mnt/scratch` if no Local SSD
4. **Warns about performance**: Alerts if running without Local SSD

### Script Output Examples

**With Local SSD:**
```
=== Setting up scratch disk ===
Local NVMe SSD detected - formatting and mounting...
Local SSD mounted at /mnt/disks/scratch ✓
Disk performance: NVMe (170K-660K IOPS)

=== Checking disk space ===
Disk space on /mnt/disks/scratch: 360GB available / 375GB total (4% used)
Disk space check passed ✓
```

**Without Local SSD:**
```
=== Setting up scratch disk ===
WARNING: No Local SSD detected. Using boot disk at /mnt/scratch
WARNING: Performance will be significantly slower (10-20x)
WARNING: Recommend recreating VM with --local-ssd interface=NVME

=== Checking disk space ===
Disk space on /mnt/scratch: 450GB available / 500GB total (10% used)
Disk space check passed ✓
```

## Cost Analysis

### Example: 30x WGS Pipeline

**Without Local SSD (using 500 GB SSD PD):**
- VM runtime: 6.5 hours
- VM cost: 6.5 hrs × $0.76/hr = $4.94
- SSD PD cost (500 GB for 6.5 hrs): $0.017/GB/hr × 500 × 6.5 = $0.06
- **Total**: $5.00

**With Local SSD (1x 375 GB + 100 GB boot disk):**
- VM runtime: 4 hours (1.6x faster)
- VM cost: 4 hrs × $0.76/hr = $3.04
- Local SSD cost: $0.08/GB-month = $0.00011/GB-hr × 375 × 4 = $0.17
- Boot disk cost: $0.017/GB/hr × 100 × 4 = $0.007
- **Total**: $3.22

**Savings**: $5.00 - $3.22 = **$1.78 per run (36% cost reduction!)**

## Important Limitations

### 1. Data is Ephemeral
- **Local SSD data is lost when VM stops/terminates**
- ✅ Safe for our pipeline (we upload results to GCS before deletion)
- ❌ Don't use for long-term data storage

### 2. Must Use with Specific Machine Types
- Local SSD works with: `n1-*`, `n2-*`, `c2-*`, `m1-*`, `m2-*`
- Does NOT work with: `e2-*`, `t2d-*`
- Always use `n1-standard-*` for genomics pipelines

### 3. Requires `--maintenance-policy=TERMINATE`
- VMs with Local SSD cannot live-migrate
- Already required for GPU VMs anyway ✓

### 4. Interface Type Matters
- Use `interface=NVME` for best performance (NVMe protocol)
- Old default was `interface=SCSI` (much slower)
- Always specify `--local-ssd interface=NVME`

## Verification Commands

### Check if Local SSD is Attached

```bash
# SSH into VM
gcloud compute ssh genomics-pipeline-human-vm --zone=us-central1-a

# Check for NVMe devices
lsblk
# Should show nvme0n1 if Local SSD is attached

# Check disk performance
sudo fio --filename=/mnt/disks/scratch/test --size=1G --direct=1 --rw=randrw --bs=4k --ioengine=libaio --iodepth=256 --runtime=10 --numjobs=4 --time_based --group_reporting --name=throughput-test
```

### Monitor Disk I/O During Pipeline

```bash
# In one terminal
gcloud compute ssh genomics-pipeline-human-vm --zone=us-central1-a
sudo iostat -x 2

# Watch for:
# - %util close to 100% = disk is bottleneck
# - r/s and w/s = IOPS (should be very high with Local SSD)
```

## Troubleshooting

### Issue: "No Local SSD detected"

**Check VM configuration:**
```bash
gcloud compute instances describe genomics-pipeline-human-vm --zone=us-central1-a | grep -A 5 "disks:"
```

Should show:
```yaml
disks:
- interface: NVME
  type: SCRATCH
  deviceName: local-ssd-0
```

**Solution**: Recreate VM with `--local-ssd interface=NVME`

### Issue: "mkfs.ext4: Device or resource busy"

**Cause**: NVMe device already has a filesystem or is mounted.

**Solution**:
```bash
sudo umount /dev/nvme0n1 || true
sudo mkfs.ext4 -F /dev/nvme0n1  # -F forces overwrite
```

### Issue: Pipeline slower than expected

**Check if using Local SSD:**
```bash
df -h | grep scratch
# Should show /dev/nvme0n1 or /dev/md0, NOT /dev/sda1
```

**Check I/O performance:**
```bash
sudo iostat -x 2
# Look for nvme0n1 device utilization
```

## Summary

### Quick Comparison

| Aspect | Without Local SSD | **With Local SSD** |
|--------|-------------------|-------------------|
| **Performance** | Baseline | **1.6-2x faster** |
| **Cost per run** | $5.00 | **$3.22 (36% cheaper)** |
| **Setup complexity** | Simple | **Auto-detected by script** |
| **Disk size** | Flexible | **375 GB increments** |
| **Data persistence** | Survives stop | **Ephemeral (OK for pipelines)** |

### Recommendations

✅ **DO use Local SSD for:**
- Production genomics pipelines
- Any workflow with heavy disk I/O
- Cost-sensitive workloads (faster = cheaper!)
- Datasets up to 100 GB FASTQ

⚠️ **Consider alternatives if:**
- You need >1 TB scratch space (use 3+ Local SSDs or PD)
- You want to pause VMs (Local SSD data is lost)
- You're using unsupported machine types

### Updated VM Creation Command (RECOMMENDED)

```bash
gcloud compute instances create genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-human-ready \
  --image-project=YOUR_PROJECT_ID \
  --boot-disk-size=100GB \
  --boot-disk-type=pd-ssd \
  --local-ssd interface=NVME \
  --maintenance-policy=TERMINATE \
  --preemptible \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=NA12878
```

This setup provides the **best performance at the lowest cost** for genomics pipelines! 🚀
