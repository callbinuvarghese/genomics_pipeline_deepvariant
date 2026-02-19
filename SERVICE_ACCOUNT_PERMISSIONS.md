# Service Account Permissions for Genomics Pipeline

## Overview

The genomics pipeline scripts include VM self-deletion functionality to automatically clean up resources after completion. This requires specific IAM permissions for the service account attached to the VM.

## Current Service Account

```
genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com
```

## Required Permissions for VM Self-Deletion

### Option 1: Use Predefined Role (Recommended)

Grant the **Compute Instance Admin (v1)** role:

```bash
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/compute.instanceAdmin.v1"
```

**Permissions included:**
- `compute.instances.delete` (required for self-deletion)
- `compute.instances.get`
- `compute.instances.list`
- `compute.zones.get`
- Plus other instance management permissions

### Option 2: Create Custom Role with Minimal Permissions (More Secure)

For least-privilege access, create a custom role with only the required permissions:

```bash
# Create custom role
gcloud iam roles create genomicsPipelineSelfDelete \
  --project=genomics-pipeline-poc \
  --title="Genomics Pipeline Self-Delete Role" \
  --description="Allows VM to delete itself after pipeline completion" \
  --permissions=compute.instances.delete,compute.instances.get,compute.zones.get \
  --stage=GA

# Bind custom role to service account
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="projects/genomics-pipeline-poc/roles/genomicsPipelineSelfDelete"
```

**Minimum required permissions:**
- `compute.instances.delete` - Delete the VM instance
- `compute.instances.get` - Get VM metadata (zone, name)
- `compute.zones.get` - Access zone information

## Existing Permissions Check

Check what permissions the service account currently has:

```bash
# View IAM policy for the project
gcloud projects get-iam-policy genomics-pipeline-poc \
  --flatten="bindings[].members" \
  --filter="bindings.members:serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --format="table(bindings.role)"
```

## Test VM Self-Deletion

To test if the service account has proper permissions:

```bash
# Create a test VM
gcloud compute instances create test-self-delete \
  --zone=us-central1-a \
  --machine-type=e2-micro \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --scopes=cloud-platform \
  --metadata=startup-script='#!/bin/bash
    sleep 30
    ZONE=$(curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone | cut -d/ -f4)
    VM_NAME=$(hostname)
    gcloud compute instances delete $VM_NAME --zone=$ZONE --quiet
  '

# Watch for the VM to self-delete after ~30 seconds
watch -n 2 'gcloud compute instances list --filter="name=test-self-delete"'
```

If the VM successfully deletes itself, the permissions are configured correctly.

## Troubleshooting

### Error: "Required 'compute.instances.delete' permission"

**Cause:** Service account lacks the `compute.instances.delete` permission.

**Solution:** Add the role using Option 1 or Option 2 above.

### Error: "Failed to auto-delete VM" (in pipeline logs)

**Cause:** Either:
1. Missing permissions (see above)
2. Network/API connectivity issues
3. VM is already deleted (rare race condition)

**Solution:**
1. Check service account permissions
2. Check Cloud Logging for detailed error messages:
   ```bash
   gcloud logging read "resource.type=gce_instance AND severity>=ERROR" --limit=50 --format=json
   ```

### VM Not Deleting But No Error

**Cause:** The trap command might not be executing properly.

**Debug:**
```bash
# SSH into VM while pipeline is running
gcloud compute ssh genomics-pipeline-human-vm --zone=us-central1-a

# Check if trap is set
trap -p ERR

# Check if gcloud is accessible
which gcloud
gcloud --version
```

## Security Considerations

### Principle of Least Privilege

The custom role (Option 2) follows the principle of least privilege by granting only the minimum required permissions. This is recommended for production environments.

### Scope Limitations

The service account should also have appropriate scopes. Ensure VMs are created with:

```bash
--scopes=cloud-platform
```

Or more specifically:
```bash
--scopes=https://www.googleapis.com/auth/compute
```

### Alternative: Manual Cleanup

If security policies prohibit VM self-deletion, you can:

1. **Remove the auto-delete code** from the script
2. **Set up Cloud Scheduler** to delete VMs older than X hours
3. **Use Cloud Functions** triggered by log entries to clean up completed VMs

Example Cloud Scheduler approach:

```bash
# Create a cleanup script
cat > cleanup_vms.sh <<'EOF'
#!/bin/bash
# Delete VMs older than 4 hours with specific labels
gcloud compute instances list \
  --filter="labels.pipeline=genomics AND creationTimestamp<$(date -u -d '4 hours ago' '+%Y-%m-%dT%H:%M:%S')" \
  --format="value(name,zone)" | \
while read name zone; do
  echo "Deleting $name in $zone"
  gcloud compute instances delete $name --zone=$zone --quiet
done
EOF

# Schedule to run every hour
gcloud scheduler jobs create http cleanup-genomics-vms \
  --schedule="0 * * * *" \
  --uri="https://YOUR_CLOUD_FUNCTION_URL" \
  --http-method=POST
```

## Additional Required Permissions

Beyond self-deletion, the service account also needs:

### Cloud Storage (already configured)
- `storage.buckets.get`
- `storage.objects.create`
- `storage.objects.get`
- `storage.objects.list`

Usually granted via: `roles/storage.objectAdmin`

### Cloud Logging (for Ops Agent)
- `logging.logEntries.create`
- `monitoring.timeSeries.create`

Usually granted via: `roles/logging.logWriter` and `roles/monitoring.metricWriter`

### Check All Current Roles

```bash
gcloud projects get-iam-policy genomics-pipeline-poc \
  --flatten="bindings[].members" \
  --filter="bindings.members:serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --format="yaml"
```

## Recommended Full Setup

```bash
# Service account should have these roles:
# 1. Compute Instance Admin (for self-deletion)
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/compute.instanceAdmin.v1"

# 2. Storage Object Admin (for reading/writing files)
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/storage.objectAdmin"

# 3. Logging Writer (for Ops Agent)
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/logging.logWriter"

# 4. Monitoring Metric Writer (for Ops Agent)
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/monitoring.metricWriter"
```

## Summary

For the VM self-deletion feature to work properly:

✅ **Grant `roles/compute.instanceAdmin.v1`** (or custom role with `compute.instances.delete` permission)
✅ **Ensure VM is created with `--scopes=cloud-platform`**
✅ **Test with a dummy VM** before running full pipeline
✅ **Monitor Cloud Logging** for any permission errors

The pipeline will gracefully handle permission failures by logging a warning instead of failing the pipeline, so you'll still get your results even if auto-deletion fails.
