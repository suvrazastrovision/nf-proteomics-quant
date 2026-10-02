# Automated Proteomics(DIA) Quant

[![Documentation](https://img.shields.io/badge/Documentation-Tutorial-2ea44f?logo=readthedocs&logoColor=white)](docs/tutorial.md)
[![AWS](https://img.shields.io/badge/AWS-EC2%20%7C%20S3-FF9900?logo=amazonwebservices&logoColor=white)](https://aws.amazon.com/)
[![Ubuntu](https://img.shields.io/badge/Ubuntu-24.04-E95420?logo=ubuntu&logoColor=white)](https://ubuntu.com/)
[![DIA-NN](https://img.shields.io/badge/DIA--NN-GitHub-181717?logo=github&logoColor=white)](https://github.com/vdemichev/DiaNN)
[![Nextflow](https://img.shields.io/badge/Nextflow-Workflow-23BFC2?logo=nextflow&logoColor=white)](https://www.nextflow.io/)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.23065890.svg)](https://doi.org/10.5281/zenodo.23065890)

This repo presents a reproducible cloud-based workflow for automated DIA proteomics quantification using AWS, Nextflow and DIA-NN. By automating data processing, workflow execution, and result generation, the framework reduces manual intervention, minimizes processing errors, improves reproducibility, and enables consistent analysis across large proteomics datasets. The approach provides a scalable foundation for standardized and efficient high-throughput proteomics analysis.

```mermaid
flowchart LR
    A[Local computer] -->|RAW and FASTA| B[(Amazon S3)]
    B -->|IAM role| C[Amazon EC2]
    C --> D[Nextflow]
    D --> E[nf-diann]
    E --> F[DIA-NN quantification]
    F --> G[Reports and logs]
    G -->|Archive| B
    B -->|Download| A
```

