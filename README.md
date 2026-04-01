# NUBE
## Introduction
NUBE (NUBE is an Ultrafast Bioinformatic analyses calculation Engine) is the ECDC's calculation engine for genomic analyses. It relies on containerization and version control to ensure reproducibility, and on a workflow management system to guarantee portability, orchestration, parallelization, and robustness.
The system can leverage cloud infrastructure for scalability, although this is not a strict requirement. NUBE is also modular by design, allowing it to be adapted to different analytical needs.

The current implementation relies on Docker for containerization, git for version control, NextFlow as workflow management system and Microsoft Azure as cloud provider.

### A modular system 
NUBE was designed following a modular approach. Currently there is a main module (main.nf), which is designed to handle genome assemblies and cgMLST allele calling operations. Other specialised modules were created for other type of analyses including QC, AMR or advanced typing analyses. Such modules can be found on the `modules/` folder. The idea is that the user should be able to strip part of the workflow to use it for individual task or as a base for ad-hoc analysis.

### Signals
Analyses are triggered through signals, which are .json files which follow a defined specification. NUBE constantly checks for signals appearing on `data/signals/` and triggers the appropriate analyses based on type of input data and the expected outputs, which are defined the "experiment_list" portion of the signal. Here is an example of a signal:

```
{
    "C834F03B-3894-5360-8C00-EABC8F3064D8": {
        "sequencing_technology": [
            "ILLUMINA"
        ],
        "project": "LEGIISO",
        "organism": "LEGIISO",
        "experiment_list": [
            "assembly",
            "allele_call",
            "qc"
        ],
        "schemas": [
            "IncHI1B_pNDM-MAR"
        ],
        "reads": [
            [
                "S3|ecdc-epc-prod/LEGIISO/FR/025171512401_S11_R1.fastq.gz",
                "S3|ecdc-epc-prod/LEGIISO/FR/025171512401_S11_R2.fastq.gz"
            ]
        ]
    }
}
```

The `schemas` input is optional and may include the name of one or more schemas to use in allele calling instead of the default schema for the organism.
