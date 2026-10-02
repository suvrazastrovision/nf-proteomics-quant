# AWS DIA Proteomics Tutorial

[Back to project overview](../README.md)

Complete the sections in order. Use placeholders for bucket names, profiles, SSH keys, and IP addresses. Never commit credentials or private scientific data.

## Project contribution

This project documents the AWS execution and reproducibility layer's validated design. The contribution covers secure S3-to-EC2 data movement, least-privilege IAM access, repeatable EC2 setup, bounded Nextflow resources, automated DIA quantification, execution reporting, provenance, checksums, and public validation data. The nf-diann workflow and DIA-NN engine are upstream open-source projects. 

## Validated environment

| Component | Validated value |
|---|---|
| AWS Region | `eu-central-1` |
| Operating system | Ubuntu 24.04 LTS, x86_64 |
| EC2 example | `m7i.4xlarge` — 16 vCPU, 64 GiB RAM |
| Root storage | 100 GiB encrypted gp3 EBS |
| Nextflow | `26.04.6` |
| nf-diann | `0.4`, commit `0446401d9d62b1d2acda186f697c6ef14a982327` |
| DIA-NN | `2.6.1` in `ghcr.io/lehtiolab/nf-diann:0.3` |

The EC2 size is a validated example for this test. Actual memory, storage, runtime, and cost depend on the dataset.

## Published demonstration

The validated test input, complete result archive, and SHA-256 checksums are published on Zenodo:

**DOI:** [10.5281/zenodo.23065890](https://doi.org/10.5281/zenodo.23065890)

- `testsample.raw` - DIA-MS test input
- `run-001-results.zip` - DIA-NN outputs, reports, logs, and metadata
- `SHA256SUMS.txt` - integrity checksums

Large scientific files remain outside Git.

## Repository layout

```text
.
├── configuration/
│   ├── base.config
│   └── ec2.config
├── docs/
│   └── tutorial.md
├── scripts/
│   └── ec2-setup.sh
├── assets/
├── lib/
├── workflows/
├── bin/
├── main.nf
├── modules.nf
├── nextflow.config
├── nextflow_schema.json
├── LICENSE
└── README.md
```

Runtime data, `.nextflow*`, `work/`, results, RAW files, FASTA files, private keys, and AWS credentials must not be committed.

## Prerequisites

- AWS account with MFA enabled
- private S3 bucket
- EC2 SSH key pair
- EC2 IAM role allowed to list the bucket, read `input/` and `reference/`, and write `results/`
- Thermo `.raw` file and matching protein FASTA
- AWS CLI on the local computer

Never commit AWS credentials, SSH keys, account IDs, private data, or presigned URLs.

## 1. Upload data to S3

Set local WSL variables:

```bash
export AWS_REGION="eu-central-1"
export S3_BUCKET="your-unique-bucket-name"
export RUN_ID="run-001"
```

Use this bucket layout:

```text
s3://your-unique-bucket-name/
|-- input/
|-- reference/
`-- results/
```

Upload and verify:

```bash
aws s3 cp testsample.raw "s3://$S3_BUCKET/input/testsample.raw"
aws s3 cp mouse_reference.fasta "s3://$S3_BUCKET/reference/mouse_reference.fasta"
aws s3 ls "s3://$S3_BUCKET/input/"
aws s3 ls "s3://$S3_BUCKET/reference/"
```

Keep S3 Block Public Access enabled. EC2 accesses the objects through its IAM role.

## 2. Launch EC2

Launch Ubuntu 24.04 x86_64 in the bucket's Region. The validated test used:

- `m7i.4xlarge` with 16 vCPU and 64 GiB RAM;
- 100 GiB encrypted gp3 root storage;
- SSH restricted to the operator's current IP;
- the least-privilege S3 IAM role as its instance profile;
- IMDSv2 required;
- Delete on termination enabled for temporary root storage.

Connect from WSL:

```bash
chmod 600 ~/.ssh/your-key.pem
ssh -i ~/.ssh/your-key.pem ubuntu@EC2_PUBLIC_IP
```

## 3. Install the environment

After copying or cloning this repository to `~/dia`:

```bash
cd ~/dia
bash scripts/ec2-setup.sh
```

Install AWS CLI v2 on the x86_64 instance:

```bash
cd /tmp
curl -fsSLo awscliv2.zip https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
unzip -q awscliv2.zip
sudo ./aws/install
rm -rf aws awscliv2.zip
```

Exit SSH and reconnect so Docker group membership takes effect, then verify:

```bash
export PATH="$HOME/.local/bin:$PATH"
java -version
nextflow -version
docker run --rm hello-world
aws sts get-caller-identity
```

The AWS identity should be an assumed EC2 role. EC2 needs no static access key.

## 4. Download inputs to EC2

```bash
mkdir -p ~/proteomics/{input,reference,results,work}

aws s3 cp "s3://$S3_BUCKET/input/testsample.raw" \
  ~/proteomics/input/testsample.raw

aws s3 cp "s3://$S3_BUCKET/reference/mouse_reference.fasta" \
  ~/proteomics/reference/mouse_reference.fasta

test -s ~/proteomics/input/testsample.raw && echo "RAW OK"
grep -m 1 '^>' ~/proteomics/reference/mouse_reference.fasta
```

## 5. Create the nf-diann input definition

Create `~/proteomics/input/input_definition.tsv` as a tab-separated file:

```tsv
file_path	sample	create_lib	train_quantums
/home/ubuntu/proteomics/input/testsample.raw	Mouse_sample	1	1
```

`file_path` identifies the EC2-local RAW file. `create_lib=1` uses it for empirical-library creation, and `train_quantums=1` uses it for quantUMS training.

Verify the tabs:

```bash
cat -T ~/proteomics/input/input_definition.tsv
```

Each separator should display as `^I`.

## 6. Run automated quantification

```bash
export RUN_ID="run-001"
mkdir -p \
  ~/proteomics/results/$RUN_ID/{diann,nextflow,metadata} \
  ~/proteomics/work/$RUN_ID

tmux new -s nfdiann
```

Run from the repository root:

```bash
cd ~/dia
export PATH="$HOME/.local/bin:$PATH"

nextflow \
  -log "$HOME/proteomics/results/$RUN_ID/metadata/.nextflow.log" \
  run . \
  -profile ec2 \
  -resume \
  --name "$RUN_ID" \
  --input "$HOME/proteomics/input/input_definition.tsv" \
  --tdb "$HOME/proteomics/reference/mouse_reference.fasta" \
  --ms1acc 10 \
  --ms2acc 10 \
  --ec2_cpus 16 \
  --ec2_memory '60 GB' \
  --output_pred_lib \
  --output_emp_lib \
  --outputquant \
  --outputreport \
  -work-dir "$HOME/proteomics/work/$RUN_ID" \
  -output-dir "$HOME/proteomics/results/$RUN_ID/diann" \
  -with-report "$HOME/proteomics/results/$RUN_ID/nextflow/execution-report.html" \
  -with-timeline "$HOME/proteomics/results/$RUN_ID/nextflow/timeline.html" \
  -with-dag "$HOME/proteomics/results/$RUN_ID/nextflow/dag.html"
```

The `10 ppm` settings match this Orbitrap test and must be reviewed for other instruments or acquisition methods.

Detach from tmux with `Ctrl+B`, then `D`. Reattach with:

```bash
tmux attach -t nfdiann
```

Monitor from another SSH session:

```bash
tail -f ~/proteomics/results/$RUN_ID/metadata/.nextflow.log
```

## 7. Inspect results

The requested run can produce:

- `report.html` and `report_small.html`;
- `report.parquet` and quantitative TSV matrices;
- predicted and empirical spectral libraries;
- per-file `.quant` files;
- DIA-NN logs;
- Nextflow report, timeline, DAG, trace, and `.nextflow.log`.

```bash
find ~/proteomics/results/$RUN_ID -type f -printf '%P\n' | sort
```

Open `diann/report.html` for the scientific report and `nextflow/execution-report.html` for workflow performance.

## 8. Preserve provenance

Record the run date and ID, nf-diann commit, Nextflow version, DIA-NN and container versions, EC2 type, parameters, input filenames, reference FASTA, and checksums under `results/$RUN_ID/metadata/`.

```bash
sha256sum \
  ~/proteomics/input/testsample.raw \
  ~/proteomics/reference/mouse_reference.fasta \
  > ~/proteomics/results/$RUN_ID/metadata/input-SHA256SUMS.txt
```

## 9. Archive results in S3

```bash
aws s3 sync \
  "$HOME/proteomics/results/$RUN_ID/" \
  "s3://$S3_BUCKET/results/$RUN_ID/"

aws s3 ls "s3://$S3_BUCKET/results/$RUN_ID/" \
  --recursive --summarize
```

Verify this archive before terminating EC2.

## 10. Download results locally

From local WSL:

```bash
export AWS_PROFILE="your-local-login-profile"
export AWS_REGION="eu-central-1"
export S3_BUCKET="your-unique-bucket-name"
export RUN_ID="run-001"

LOCAL_RESULTS="/mnt/c/Users/YOUR_WINDOWS_USER/Desktop/proteomics_DIA/$RUN_ID"
mkdir -p "$LOCAL_RESULTS"

aws s3 sync \
  "s3://$S3_BUCKET/results/$RUN_ID/" \
  "$LOCAL_RESULTS/"

test -s "$LOCAL_RESULTS/diann/report.html" && echo "RESULTS OK"
```

## 11. Finish the run

After confirming both the S3 archive and local download:

1. terminate compute that is no longer needed;
2. confirm that no unattached EBS volume remains;
3. retain required S3 inputs and results;
4. review old S3 objects and versions periodically.

Stopping EC2 does not eliminate every charge. EBS volumes, snapshots, Elastic IP addresses, and S3 objects may continue to incur costs.

## Reusing `.quant` files

nf-diann can reuse `.quant` files through the input definition's `quantfile` column or `--quantdir`. RAW and `.quant` basenames must correspond. Reuse only files produced with compatible DIA-NN, library, FASTA, and analysis settings, and record their provenance.

## Troubleshooting

| Problem | Check |
|---|---|
| Docker permission denied | Reconnect after setup, then run `docker run --rm hello-world` |
| `nextflow: command not found` | Run `export PATH="$HOME/.local/bin:$PATH"` |
| S3 access denied | Run `aws sts get-caller-identity` and inspect the EC2 role policy |
| Input-definition error | Run `cat -T input_definition.tsv` and verify EC2-local paths |
| Interrupted run | Repeat the identical command with `-resume` and the same work directory |
| Out of memory or disk | Inspect `free -h`, `df -h`, and `docker system df` |

## Security decisions

- EC2 uses an IAM role instead of stored access keys.
- S3 Block Public Access remains enabled for operational data.
- SSH is limited to the operator's current IP.
- IMDSv2 and EBS encryption are enabled.
- Credentials, private keys, RAW data, FASTA files, work directories, and results must stay outside Git.
- Public sharing uses Zenodo and a persistent DOI instead of a public operational bucket.

## Expected completion state

The run is complete when the scientific report opens locally, the complete result directory exists in S3, the provenance files and checksums are retained, and unnecessary EC2/EBS resources are terminated.

## Attribution and citation

This AWS workflow was developed and validated by **Suvra Nath** as a portfolio project in cloud bioinformatics and automated proteomics quantification.

Scientific execution uses [lehtiolab/nf-diann](https://github.com/lehtiolab/nf-diann), [DIA-NN](https://github.com/vdemichev/DiaNN), and [Nextflow](https://www.nextflow.io/). Their respective authorship and citation guidance apply.

For the demonstration dataset and validated results, cite:

> Suvra Nath. *Automated Proteomics(DIA) Quant: Cloud Platform Engineering*. Zenodo. https://doi.org/10.5281/zenodo.23065890

The upstream copyright notice and MIT terms remain in [LICENSE](../LICENSE).
