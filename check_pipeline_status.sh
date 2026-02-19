#!/bin/bash
# Check pipeline status without SSH

echo "=== VM Status ==="
gcloud compute instances describe genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --format="value(status)"

echo ""
echo "=== Latest Pipeline Output (last 50 lines) ==="
gcloud compute instances get-serial-port-output genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --port=1 2>/dev/null | tail -50

echo ""
echo "=== GCS Output Files ==="
gsutil ls -lh gs://genomics-output-bucket-binuv/human/ 2>/dev/null || echo "No output files yet"

echo ""
echo "=== Recent Log Files ==="
gsutil ls -lh gs://genomics-ref-bucket-binuv/logs/ 2>/dev/null | tail -5 || echo "No log files yet"

