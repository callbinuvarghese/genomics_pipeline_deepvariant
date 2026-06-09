# VM Scopes and Permissions Guide

## Overview

GCP VMs require both **Service Account IAM roles** AND **OAuth scopes** to access Google Cloud APIs. Both must be configured correctly for the pipeline to work.

## Current Configuration ✅

Our VM creation command uses:
```bash
--scopes=cloud-platform
```

This is **shorthand** for:
```bash
--scopes=https://www.googleapis.com/auth/cloud-platform
```

Both are equivalent and provide full access to all Google Cloud APIs (limited by IAM roles).

## Understanding Scopes vs IAM Roles

### Scopes (OAuth 2.0)
- **What**: OAuth 2.0 scopes that define the maximum permissions a VM can request
- **Where**: Configured at VM creation time (`--scopes` flag)
- **Cannot be changed**: After VM is created, scopes are immutable
- **Think of it as**: "What the VM *can ask for*"

### IAM Roles
- **What**: Actual permissions granted to the service account
- **Where**: Configured in IAM policy (`gcloud projects add-iam-policy-binding`)
- **Can be changed**: Anytime, even while VM is running
- **Think of it as**: "What the VM *is allowed to do*"

### The Rule
**Effective permissions = MIN(scopes, IAM roles)**

Both must permit an action for it to work!

## Scope Options

### Option 1: `cloud-platform` (Current - RECOMMENDED)

```bash
--scopes=cloud-platform
# OR
--scopes=https://www.googleapis.com/auth/cloud-platform
```

**Pros:**
- ✅ Works with all Google Cloud services
- ✅ Simple - one scope covers everything
- ✅ No need to update if pipeline needs change
- ✅ Security still controlled by IAM roles

**Cons:**
- ⚠️ Grants broad scope (but still limited by IAM roles)

**Use when**: You want flexibility and will control access via IAM roles (our case)

### Option 2: Specific Scopes (Most Restrictive)

```bash
--scopes=storage-rw,compute-rw,logging-write,monitoring-write
# OR full URLs:
--scopes=https://www.googleapis.com/auth/devstorage.read_write,https://www.googleapis.com/auth/compute,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write
```

**Pros:**
- ✅ Most restrictive (defense in depth)
- ✅ Explicit about what APIs are accessible

**Cons:**
- ❌ Must update if pipeline needs change
- ❌ Harder to debug scope vs IAM issues
- ❌ More verbose

**Use when**: Maximum security is required, pipeline is stable

### Option 3: Default Scopes (NOT RECOMMENDED)

```bash
# No --scopes flag (uses default)
# Default = storage-ro, logging-write, monitoring-write, service-control, service-management
```

**Pros:**
- ⚠️ None for our use case

**Cons:**
- ❌ `storage-ro` is read-only (pipeline needs write!)
- ❌ Missing compute scope (can't self-delete)
- ❌ Will cause pipeline failures

**Use when**: Never for our pipeline!

## Required Scopes for Our Pipeline

### Minimum Required Scopes

| Service | Scope | Why Needed |
|---------|-------|------------|
| **Cloud Storage** | `storage-rw` or `devstorage.read_write` | Download FASTQ, upload results |
| **Compute Engine** | `compute-rw` or `compute` | VM self-deletion |
| **Cloud Logging** | `logging-write` | Ops Agent logs |
| **Cloud Monitoring** | `monitoring-write` | Ops Agent metrics |

### Recommended Configuration

**Production (Recommended):**
```bash
--scopes=cloud-platform
```
✅ Simple, flexible, secure (when combined with proper IAM roles)

**High-Security Environment:**
```bash
--scopes=storage-rw,compute-rw,logging-write,monitoring-write
```
✅ Most restrictive, explicit

## Required IAM Roles

The service account also needs these IAM roles:

```bash
# Storage access
roles/storage.objectAdmin  # OR roles/storage.admin (superset)

# VM self-deletion
roles/compute.instanceAdmin.v1

# Ops Agent
roles/logging.logWriter
roles/monitoring.metricWriter
```

**Note:** `roles/storage.admin` includes all permissions from `roles/storage.objectAdmin` plus bucket-level operations. Either role works for the pipeline.

See [SERVICE_ACCOUNT_PERMISSIONS.md](SERVICE_ACCOUNT_PERMISSIONS.md) for details.

## Verification Commands

### Check VM Scopes

```bash
# List scopes for running VM
gcloud compute instances describe genomics-pipeline-human-vm \
  --zone=us-central1-a \
  --format="value(serviceAccounts[0].scopes)"
```

**Expected output** (with `--scopes=cloud-platform`):
```
https://www.googleapis.com/auth/cloud-platform
```

### Check Service Account IAM Roles

```bash
# List roles for service account
gcloud projects get-iam-policy genomics-pipeline-poc \
  --flatten="bindings[].members" \
  --filter="bindings.members:serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --format="table(bindings.role)"
```

**Expected output:**
```
ROLE
roles/compute.instanceAdmin.v1
roles/logging.logWriter
roles/monitoring.metricWriter
roles/storage.objectAdmin  # OR roles/storage.admin (superset with all objectAdmin permissions)
```

### Test Access from VM

```bash
# SSH into VM
gcloud compute ssh genomics-pipeline-human-vm --zone=us-central1-a

# Test Cloud Storage access
gsutil ls gs://genomics-input-bucket-binuv/
# Should succeed

# Test Compute Engine access
gcloud compute instances list --limit=1
# Should succeed

# Test token scopes
gcloud auth list
gcloud auth application-default print-access-token
curl -H "Metadata-Flavor: Google" \
  http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/scopes
```

## Common Issues and Solutions

### Issue 1: "403 Forbidden" when accessing GCS

**Possible causes:**
1. Missing storage scope
2. Service account lacks IAM role
3. Wrong bucket permissions

**Debug:**
```bash
# Check scopes
gcloud compute instances describe VM_NAME --zone=ZONE \
  --format="value(serviceAccounts[0].scopes)"

# Check IAM roles
gcloud projects get-iam-policy PROJECT_ID \
  --flatten="bindings[].members" \
  --filter="bindings.members:serviceAccount:SERVICE_ACCOUNT"

# Check bucket IAM
gsutil iam get gs://BUCKET_NAME
```

**Solution:**
```bash
# If scope is missing: Recreate VM with correct scopes
# If IAM role is missing:
gcloud projects add-iam-policy-binding PROJECT_ID \
  --member="serviceAccount:SERVICE_ACCOUNT" \
  --role="roles/storage.objectAdmin"
```

### Issue 2: "403 Forbidden" when deleting VM

**Possible cause:** Missing compute scope or IAM role

**Solution:**
```bash
# Grant IAM role
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/compute.instanceAdmin.v1"

# If scope is missing: Recreate VM with --scopes=cloud-platform
```

### Issue 3: Ops Agent can't send metrics

**Possible cause:** Missing logging/monitoring scopes

**Solution:**
```bash
# Grant IAM roles
gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/logging.logWriter"

gcloud projects add-iam-policy-binding genomics-pipeline-poc \
  --member="serviceAccount:genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com" \
  --role="roles/monitoring.metricWriter"

# Recreate VM with correct scopes if needed
```

## Scope Comparison Table

| Scope Alias | Full Scope URL | Purpose |
|-------------|----------------|---------|
| `cloud-platform` | `https://www.googleapis.com/auth/cloud-platform` | Full access to all GCP APIs |
| `storage-rw` | `https://www.googleapis.com/auth/devstorage.read_write` | Read/write Cloud Storage |
| `storage-ro` | `https://www.googleapis.com/auth/devstorage.read_only` | Read-only Cloud Storage |
| `storage-full` | `https://www.googleapis.com/auth/devstorage.full_control` | Full control Cloud Storage |
| `compute-rw` | `https://www.googleapis.com/auth/compute` | Compute Engine read/write |
| `compute-ro` | `https://www.googleapis.com/auth/compute.readonly` | Compute Engine read-only |
| `logging-write` | `https://www.googleapis.com/auth/logging.write` | Write logs |
| `monitoring` | `https://www.googleapis.com/auth/monitoring` | Monitoring (read/write) |
| `monitoring-write` | `https://www.googleapis.com/auth/monitoring.write` | Write metrics only |

## Best Practices

### 1. Use `cloud-platform` for Development/Testing
```bash
--scopes=cloud-platform
```
- ✅ Easy to debug
- ✅ No scope-related issues
- ✅ Security still enforced by IAM

### 2. Consider Specific Scopes for Production
```bash
--scopes=storage-rw,compute-rw,logging-write,monitoring-write
```
- ✅ Principle of least privilege
- ✅ Defense in depth
- ⚠️ Requires maintenance if pipeline changes

### 3. Never Use Default Scopes
- ❌ Missing required scopes
- ❌ Will cause failures

### 4. Document Scope Requirements
- 📝 List required scopes in README
- 📝 Include scope verification in setup guide
- 📝 Add scope check to pre-flight tests

### 5. Audit Regularly
```bash
# Monthly audit: Check all VMs have correct scopes
gcloud compute instances list --format="table(name,zone,serviceAccounts[0].scopes)"
```

## Updated VM Creation Commands

### Option 1: Simple (Recommended)

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
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded.sh,SAMPLE_ID=NA12878
```

### Option 2: Explicit Scopes (High Security)

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
  --scopes=storage-rw,compute-rw,logging-write,monitoring-write \
  --service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com \
  --metadata=startup-script-url=gs://genomics-ref-bucket-binuv/scripts/genomics_pipeline_human_gpu_preloaded.sh,SAMPLE_ID=NA12878
```

Both work! Option 1 is simpler, Option 2 is more restrictive.

## Security Considerations

### Scopes are NOT Security Boundaries

**Common misconception:**
> "If I use specific scopes, the VM is more secure"

**Reality:**
- Scopes only limit what APIs can be called
- **IAM roles determine actual permissions**
- A VM with `cloud-platform` scope but no IAM roles can do nothing
- A VM with specific scopes but broad IAM roles can do a lot

**Best practice:**
1. Use `cloud-platform` scope for simplicity
2. Grant minimal IAM roles (principle of least privilege)
3. Use separate service accounts for different workloads
4. Regularly audit IAM policies

### Defense in Depth

For maximum security, layer controls:
1. ✅ Specific scopes (`storage-rw,compute-rw,logging-write,monitoring-write`)
2. ✅ Minimal IAM roles (only what's needed)
3. ✅ VPC Service Controls (if available)
4. ✅ Organization policies
5. ✅ Audit logging enabled

## Summary

### Current Configuration ✅

```bash
--scopes=cloud-platform
--service-account=genomics-pipeline@genomics-pipeline-poc.iam.gserviceaccount.com
```

**Is this secure?** YES!
- Scopes grant broad API access
- IAM roles control actual permissions
- Service account has only required roles
- Pipeline runs with least privilege

### Quick Checklist

- [x] VM created with `--scopes=cloud-platform` ✅
- [x] Service account has required IAM roles ✅
- [x] Scopes verified with `gcloud compute instances describe` ✅
- [x] IAM roles verified with `gcloud projects get-iam-policy` ✅
- [x] Pipeline tested and working ✅

**The current scope configuration is correct and secure!** 🔒

---

**Related Documentation:**
- [SERVICE_ACCOUNT_PERMISSIONS.md](SERVICE_ACCOUNT_PERMISSIONS.md) - IAM roles setup
- [IMPROVEMENTS_SUMMARY.md](IMPROVEMENTS_SUMMARY.md) - Overall pipeline improvements
- [LOCAL_SSD_SETUP.md](LOCAL_SSD_SETUP.md) - Local SSD configuration
