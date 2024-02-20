# DoBSeq Nextflow Pipeline

Usage:

1. Install nextflow [manually](https://www.nextflow.io/docs/latest/getstarted.html) or using conda:

```Bash
conda install nextflow
```

2. Clone repository

```Bash
git clone https://github.com/madscort/DoBSeqWF.git .
```

3. Run pipeline without input data:

```Bash
nextflow run main.nf -profile (standard/esrum/ngc) -stub
```

3. Run pipeline with tiny test data:

```Bash
nextflow run main.nf -profile (standard/esrum/ngc),test
```

4. Run pipeline with input data (see test data for file contents):

```Bash
nextflow run main.nf                                       \
  -profile (standard/esrum/ngc)                            \
  --pooltable <path to pool fastq file table>              \
  --decodetable <path to pool decode tsv>                  \
  --reference_genome <path to indexed reference genome>    \
  --bedfile <path to bedfile with target regions>          \
  --ploidy <integer>
```

Workflow repository contents:

```Bash
DoBSeqWF                                    
├── LICENSE
├── README.md
├── assets
│   ├── data
│   │   ├── reference_genomes
│   │   │   └── small
│   │   │       └── small_reference.*
│   │   └── test_data
│   │       ├── decodetable.tsv
│   │       ├── pools
│   │       │   └── *.fq.gz
│   │       ├── pooltable.tsv
│   │       ├── snvlist.tsv
│   │       └── target_calling.bed
│   └── helper_scripts
│       └── simulator.py                  # Script for simulating minimal pipeline data
├── bin                                   # Executable pipeline scripts
│   └── <script>.*
├── conf
│   └── profiles.config                   # Configuration profiles for compute environments
├── envs
│   └── <name>/
│       └── environment.yaml              # Conda environment definitions
├── main.nf                               # Main workflow
├── modules/
│   └── <module>.nf                       # Module scripts
└── nextflow.config                       # Workflow parameters
```

