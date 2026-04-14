# NUBE
## Introduction
NUBE (NUBE is an Ultrafast Bioinformatic analyses calculation Engine) is the ECDC's calculation engine for genomic analyses. It relies on containerization and version control to ensure reproducibility, and on a workflow management system to guarantee portability, orchestration, parallelization, and robustness.
The system can leverage cloud infrastructure for scalability, although this is not a strict requirement. NUBE is also modular by design, allowing it to be adapted to different analytical needs.

The current implementation relies on Docker for containerization, git for version control, NextFlow as workflow management system and Microsoft Azure as cloud provider.

### A modular system 
NUBE was designed following a modular approach. Currently there is a main module (main.nf), which is designed to handle genome assemblies and cgMLST allele calling operations. Other specialised modules were created for other type of analyses including QC, AMR or advanced typing analyses. Such modules can be found on the `modules/` folder. The idea is that the user should be able to strip part of the workflow to use it for individual task or as a base for ad-hoc analysis.

### Signals
Analyses are triggered through signals in json format. NUBE constantly checks for signals appearing on `data/signals/` and triggers the appropriate analyses based on type of input data and the expected outputs, which are defined in the "experiment_list" portion of the signal. 

Each json should contain the information needed to run one sample, with the sample ID specified as the key as seen below in the signal example.

Mandatory fields in the json input:
- Defined input in the format of accessions, reads, sequences, or assembly:
  - `accessions`*: SRA or ENA accession to start analysis from publically available reads.
  - `reads`*: Path to fastq file(s), available locally or in cloud.
  - `assembly`: Path to assembly, available locally or in cloud.
  - `sequences`: Path to sequences, available locally or in cloud. Currently used for the experiment `add_to_alignment`.
- `sequencing_technology`*: Sequencing technology used. [ILLUMINA, IONTORRENT, ONT or PACBIO]
- `organism`: Which organism to tailor the analysis for. Available options can be found in the [settings](https://github.com/EU-ECDC/ecdc_nube/blob/main/settings.json).
- `experiment_list`*: Specify one or several post-assembly analyses to perform. The following analyses are available:
  - `qc`: Run quality check on the assembly.
  - `allele_call`: Run allele calling with cgMLST or custom schemas.
  - `add_to_alignment`: Add sequences to preexisting phylogenetic tree (HAVISO and POLIISO only).
  - `amrfinder`: Run antimicrobial resistance gene profiling (CAMPISO, ECOLIISO, SALMISO only).
  - `mlst`: Run multilocus sequence typing (CAMPISO only).
- `project`: The project name that will be used for naming of the subdirectory for analysis output.

Optional fields in the json input:
- `schemas`*: The name of one (or more) schemas to use in allele calling instead of the default schema for the organism.
- `flags`*: Not yet supported.

\* Should be specified as a list.

An example of a signal:
```
{
    "C834F03B-3894-5360-8C00-EABC8F3064D8": {
        "sequencing_technology": [
            "ILLUMINA"
        ],
        "project": "LEGIISO",
        "organism": "LEGIISO",
        "experiment_list": [
            "allele_call",
            "qc"
        ],
        "schemas": [
            "IncHI1B_pNDM-MAR"
        ],
        "reads": [
            [
                "S3|path/to/LEGIISO/025171512401_S11_R1.fastq.gz",
                "S3|path/to/LEGIISO/025171512401_S11_R2.fastq.gz"
            ]
        ]
    }
}
```

Depending on the type of data, the input can be formatted as any of the below options:
```
        "reads": [
            [
                "S3|path/to/LEGIISO/025171512401_S11_R1.fastq.gz",
                "S3|path/to/LEGIISO/025171512401_S11_R2.fastq.gz"
            ]
        ]
```

```
        "accessions": [
            "SRR000001"
        ]
```

```
        "assembly": "az://path/to/LEGIISO/assemblies/025171512401.fasta"
```

```
        "sequences": "az://path/to/LEGIISO/sequences/025171512401.fasta"
```
