# Checkpoint and Resume Guide

## Overview

The checkpointed version of the pipeline (`genomics_pipeline_human_gpu_preloaded_checkpointed.sh`) includes automatic checkpoint/resume functionality to handle VM preemptions gracefully.

## How It Works

### Checkpoint Creation

After completing Step 3 (Mark Duplicates), the pipeline:
1. ✅ Uploads `${SAMPLE_ID}.sorted.dedup.bam` to `gs://genomics-ref-bucket-binuv/checkpoints/`
2. ✅ Uploads `${SAMPLE_ID}.sorted.dedup.bam.bai` (BAM index)
3. ✅ Continues to Step 4 (Variant Calling)

### Resume Detection

On pipeline start, the script:
1. Checks if checkpoint exists in GCS
2. Verifies both BAM and index files are present
3. If valid checkpoint found:
   - Downloads checkpoint files
   - **Skips Steps 1-3** (saves ~2 hours for 30x WGS!)
   - Jumps directly to Step 4 (Variant Calling)

### Automatic Cleanup

On successful completion:
- Checkpoint files are automatically deleted from GCS
- Saves storage costs (~$0.002/GB-month)

## When Resume Happens

### Scenario 1: VM Preempted During Step 3

```
Run 1 (Preempted):
  ✓ Step 1: Alignment (complete)
  ✓ Step 2: Fixmate+Sort (complete)
  ✓ Step 3: Mark Duplicates (complete)
  ✓ Checkpoint uploaded
  ✗ Step 4: Variant Calling (preempted)

Run 2 (Resumed):
  ✓ Checkpoint detected!
  ⏭ Skipping Steps 1-3
  ✓ Step 4: Variant Calling (complete)
  ✓ Results uploaded
  ✓ Checkpoint cleaned up
```

**Time saved**: ~2 hours (for 30x WGS)

### Scenario 2: VM Preempted During Step 1 or 2

```
Run 1 (Preempted):
  ✗ Step 1: Alignment (interrupted)

Run 2 (Full Run):
  ✗ No checkpoint detected
  ✓ Full pipeline from beginning
```

**Time saved**: None (must restart from beginning)

### Scenario 3: Successful First Run

```
Run 1 (Success):
  ✓ All steps complete
  ✓ Checkpoint cleaned up automatically
```

**No resume needed**

## Usage

### Standard Usage (Preemptible VM)

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

**If preempted**: Simply run the same command again with the same `SAMPLE_ID`!

### Manual Resume

If you need to manually trigger a resume:

```bash
# 1. Check if checkpoint exists
gsutil ls gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam*

# 2. If exists, just run the same VM creation command
# The script will automatically detect and resume
```

### Force Fresh Start (Ignore Checkpoint)

To start from scratch even if checkpoint exists:

```bash
# Delete the checkpoint first
gsutil rm gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam*

# Then create VM as usual
```

Or use a different `SAMPLE_ID`:
```bash
--metadata=...,SAMPLE_ID=NA12878_run2
```

## Log Output Examples

### First Run (Creating Checkpoint)

```
=== Checking for previous checkpoint ===
No checkpoint found for sample NA12878
Will run full pipeline from the beginning

=== Downloading FASTQ files ===
FASTQ download completed in 126 seconds

=== Steps 1-2: Alignment → Sort → Fixmate → Sort (piped) ===
Fully-piped alignment and sorting completed in 3748 seconds

=== Step 3: Marking Duplicates ===
Duplicate marking completed in 285 seconds

=== CHECKPOINT: Uploading deduplicated BAM to GCS ===
Checkpoint uploaded in 45 seconds

=== Step 4: Variant Calling with DeepVariant (GPU) ===
[PREEMPTED HERE - VM dies]
```

### Second Run (Resuming from Checkpoint)

```
=== Checking for previous checkpoint ===
✓ Checkpoint found! Previous run was interrupted.
  Found: gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam
✓ Checkpoint BAM index found

RESUMING FROM CHECKPOINT: Skipping Steps 1-3, jumping to Step 4 (Variant Calling)
Downloading checkpoint files from GCS...
Checkpoint download completed in 38 seconds

=== Skipping FASTQ download (resuming from checkpoint) ===

=== Skipping Steps 1-3 (resuming from checkpoint) ===
Using previously computed BAM: NA12878.sorted.dedup.bam

=== Step 4: Variant Calling with DeepVariant (GPU) ===
Variant calling completed in 1856 seconds

=== Pipeline SUMMARY ===
FASTQ download:      0s (skipped)
Align+Sort+Fixmate:  0s (skipped)
Mark Duplicates:     0s (skipped)
Checkpoint Upload:   0s (resumed from existing)
Variant Calling:     1856s
TOTAL RUNTIME:       1900s (32 minutes)

=== Cleaning up checkpoint files ===
Pipeline completed successfully - removing checkpoint to save storage costs
```

## Benefits

### Time Savings

| Coverage | Step 1-3 Time | Step 4 Time | **Time Saved if Preempted After Step 3** |
|----------|---------------|-------------|------------------------------------------|
| 4x (test) | ~25 min | ~20 min | **~25 min (55% saved)** |
| 30x WGS | ~2.5 hours | ~1.5 hours | **~2.5 hours (63% saved)** |
| 60x WGS | ~5 hours | ~3 hours | **~5 hours (63% saved)** |

### Cost Savings

**Example: 30x WGS on preemptible VM**

**Without checkpoint** (preempted after Step 3):
- Run 1: 2.5 hrs × $0.76/hr = $1.90 (wasted)
- Run 2: 4 hrs × $0.76/hr = $3.04 (full restart)
- **Total**: $4.94

**With checkpoint** (preempted after Step 3):
- Run 1: 2.5 hrs × $0.76/hr = $1.90 (creates checkpoint)
- Run 2: 1.5 hrs × $0.76/hr = $1.14 (resumes from checkpoint)
- Checkpoint storage: $0.002/GB-month × 12 GB × 0.1 month = $0.002
- **Total**: $3.04

**Savings**: $4.94 - $3.04 = **$1.90 (38% cost reduction if preempted!)**

## Technical Details

### Checkpoint Location

```
gs://genomics-ref-bucket-binuv/checkpoints/${SAMPLE_ID}.sorted.dedup.bam
gs://genomics-ref-bucket-binuv/checkpoints/${SAMPLE_ID}.sorted.dedup.bam.bai
```

### Checkpoint Size

| Coverage | FASTQ Size | Checkpoint BAM Size | Upload Time (Local SSD) |
|----------|------------|---------------------|------------------------|
| 4x | 9.2 GB | ~10 GB | ~30 sec |
| 30x | 75 GB | ~100 GB | ~4 min |
| 60x | 150 GB | ~200 GB | ~8 min |

### Resume Detection Logic

```bash
if gsutil -q stat "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam"; then
  if gsutil -q stat "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam.bai"; then
    # Both files exist - safe to resume
    RESUME_FROM_CHECKPOINT=true
  else
    # Incomplete checkpoint - clean up and start fresh
    gsutil rm "${CHECKPOINT_BUCKET}/${SAMPLE_ID}.sorted.dedup.bam"
  fi
fi
```

**Safety check**: Both BAM and index must exist to resume

### What Gets Skipped on Resume

✅ **Skipped** (saves time):
- FASTQ download
- Alignment (Step 1)
- Fixmate + Sort (Step 2)
- Mark Duplicates (Step 3)
- Checkpoint upload

✅ **Still Runs** (required):
- Reference verification
- GPU/Docker verification
- GCS access verification
- Checkpoint download
- Variant Calling (Step 4)
- QC metrics
- Result upload
- Checkpoint cleanup

## Monitoring Checkpoints

### List All Checkpoints

```bash
gsutil ls -lh gs://genomics-ref-bucket-binuv/checkpoints/
```

### Check Specific Sample

```bash
gsutil ls -lh gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam*
```

### Verify Checkpoint Integrity

```bash
# Download and validate
gsutil cp gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam .
samtools quickcheck NA12878.sorted.dedup.bam
# No output = valid file
```

### Manual Cleanup (if needed)

```bash
# Remove all checkpoints
gsutil rm gs://genomics-ref-bucket-binuv/checkpoints/*

# Remove specific sample
gsutil rm gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam*
```

## Troubleshooting

### Issue: Checkpoint exists but resume fails

**Symptoms**:
```
✓ Checkpoint found!
ERROR: Failed to download checkpoint BAM
```

**Cause**: Checkpoint file is corrupted or incomplete

**Solution**:
```bash
# Delete corrupted checkpoint
gsutil rm gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam*

# Run pipeline again (will start fresh)
```

### Issue: Pipeline always resumes even though I want fresh start

**Cause**: Old checkpoint exists for this `SAMPLE_ID`

**Solution**:
```bash
# Option 1: Delete checkpoint
gsutil rm gs://genomics-ref-bucket-binuv/checkpoints/NA12878.sorted.dedup.bam*

# Option 2: Use different SAMPLE_ID
--metadata=...,SAMPLE_ID=NA12878_run2
```

### Issue: Resumed run has different results than first run would have

**Cause**: Resume is working correctly! The first run was interrupted.

**Note**: Resume produces identical results to a full run because:
- Same BAM file (from checkpoint)
- Same DeepVariant model and parameters
- Deterministic variant calling

**Verification**:
```bash
# Run twice with different SAMPLE_IDs, compare VCFs
diff sample1.vcf sample2.vcf
# Should be identical (ignoring timestamps)
```

## Best Practices

### 1. Use Checkpointing for Large Datasets

**Recommended for**:
- ✅ 30x+ WGS (alignment takes >2 hours)
- ✅ Preemptible VMs
- ✅ Production pipelines

**Not needed for**:
- ❌ Small test data (<20 GB FASTQ)
- ❌ Non-preemptible VMs
- ❌ One-off analyses

### 2. Monitor Checkpoint Storage Costs

```bash
# Check checkpoint bucket size
gsutil du -sh gs://genomics-ref-bucket-binuv/checkpoints/

# Set lifecycle policy to auto-delete old checkpoints
cat > lifecycle.json <<EOF
{
  "lifecycle": {
    "rule": [
      {
        "action": {"type": "Delete"},
        "condition": {"age": 7}
      }
    ]
  }
}
EOF
gsutil lifecycle set lifecycle.json gs://genomics-ref-bucket-binuv/checkpoints/
```

**Automatically deletes checkpoints older than 7 days**

### 3. Name Samples Consistently

Use unique, descriptive `SAMPLE_ID`s:
- ✅ Good: `NA12878_run1`, `patient123_wgs_30x`
- ❌ Bad: `test`, `sample1`, `NA12878` (reused)

### 4. Test Resume Logic

```bash
# Simulate preemption
gcloud compute instances create test-checkpoint-vm \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-human-ready \
  --image-project=YOUR_PROJECT_ID \
  --boot-disk-size=100GB \
  --boot-disk-type=pd-ssd \
  --local-ssd interface=NVME \
  --maintenance-policy=TERMINATE \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=test_checkpoint

# Wait for Step 3 to complete, then manually delete VM
gcloud compute instances delete test-checkpoint-vm --zone=us-central1-a --quiet

# Verify checkpoint created
gsutil ls gs://genomics-ref-bucket-binuv/checkpoints/test_checkpoint*

# Run again - should resume
gcloud compute instances create test-checkpoint-vm-resume \
  --zone=us-central1-a \
  --machine-type=n1-standard-16 \
  --accelerator=type=nvidia-tesla-t4,count=1 \
  --image-family=genomics-pipeline-human-ready \
  --image-project=YOUR_PROJECT_ID \
  --boot-disk-size=100GB \
  --boot-disk-type=pd-ssd \
  --local-ssd interface=NVME \
  --maintenance-policy=TERMINATE \
  --scopes=cloud-platform \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded_checkpointed.sh,SAMPLE_ID=test_checkpoint

# Check logs - should see "RESUMING FROM CHECKPOINT"
```

## Comparison: With vs Without Checkpointing

| Feature | Standard Script | **Checkpointed Script** |
|---------|----------------|------------------------|
| **Resume on preemption** | ❌ No | ✅ Yes |
| **Time saved if preempted** | 0% | **Up to 63%** |
| **Cost saved if preempted** | 0% | **Up to 38%** |
| **Storage overhead** | $0 | **~$0.002/sample** |
| **Script complexity** | Simple | Moderate |
| **Best for** | Small datasets, test runs | **Large datasets, production** |

## Summary

✅ **Use checkpointed version when:**
- Running large datasets (30x+ WGS)
- Using preemptible VMs
- Pipeline takes >2 hours
- Cost optimization is important

❌ **Use standard version when:**
- Small test datasets (<20 GB)
- Non-preemptible VMs
- Quick one-off analyses
- Simplicity preferred

**The checkpoint/resume feature provides automatic fault tolerance at minimal cost, making it ideal for production genomics pipelines on preemptible VMs!** 🎯
