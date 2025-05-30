# DoBSeq Workflow

### General usage:

### 1. Install nextflow [manually](https://www.nextflow.io/docs/latest/getstarted.html) or using [conda](https://www.anaconda.com/docs/getting-started/miniconda/install):

```Bash
conda install nextflow
```

### 2. In a clean folder - clone this repository:

```Bash
git clone https://github.com/RasmussenLab/DoBSeqWF.git .
```

### 3. Run pipeline without any data (dry-run) using [docker](https://docs.docker.com/desktop), [apptainer](https://apptainer.org/docs/admin/main/installation.html#installation-on-linux) or conda environments:

```Bash
# Docker
nextflow run main.nf -profile test -stub -with-docker
# Apptainer
nextflow run main.nf -profile test -stub -with-apptainer
# Conda
nextflow run main.nf -profile test -stub -with-conda
```

### 4. Run pipeline with minimal test data:

```Bash
# Docker
nextflow run main.nf -profile test -with-docker
# Apptainer
nextflow run main.nf -profile test -with-apptainer
# Conda
nextflow run main.nf -profile test -with-conda
```

### 5. Run pipeline with input data:
The pipeline has multiple optional configurations found in ```nextflow.config```.
Configurations can be supplied a configuration file, see ```config.json```, and run with ```nextflow run main.nf -params-file config.json```, or parameters can be added directly from the commandline:
```Bash
nextflow run main.nf                                       \
  (-with-docker/-with-apptainer/-with-conda)               \
  --pooltable <path to pool fastq file table>              \
  --decodetable <path to pool decode tsv>                  \
  --reference_genome <path to indexed reference genome>    \
  --bedfile <path to bedfile with target regions>          \
  --ploidy <integer>
```
The ```pooltable.tsv``` should connect (user assigned) pool id's to input FASTQ files; one entry for each pool.
```Bash
pool_row_1  path/to/sample1_R1.fq.gz  path/to/sample1_R2.fq.gz
pool_column_1  path/to/sample2_R1.fq.gz  path/to/sample2_R2.fq.gz
```
The ```decodetable.tsv``` should map (user assigned) individual id's in the matrix to the corresponding row and column id's of each pool; one entry for each element in the matrix.
```Bash
individual1  pool_row_1  pool_column_1
```
### 6. Pipeline output
The workflow will output a results folder containing multiple config dependent output files:
```Bash
results
├── pinpointables.vcf           # Merged VCF file containing all assigned variants
├── cram/                       # CRAM files for each pool
├── logs/                       # Log files for each process
├── variants/                   # VCF files for each pool
├── variant_tables/             # TSV files converted from pool VCFs
└── pinpoint_variants/
    ├── all_pins/               # All pinpointables for each sample in individual vcfs (*note)
    ├── unique_pins/            # All unique pinpointables for each sample in individual vcfs (*note)
    ├── *_merged.vcf.gz         # All pinpointables for all samples in a single vcf without sample information
    ├── summary.tsv             # Variant counts for each sample
    └── lookup.tsv              # Variant to sample lookup table
```
A central files is the ```pinpointables.vcf```. This file contains all individually assigned variants. Since each variant contains information from two pools, these a presented as the sample columns: ROW and COLUMN.

# Workflow repository contents:

```Bash
DoBSeqWF                                    
├── LICENSE
├── VERSION
├── README.md
├── assets
│   ├── data
│   │   ├── reference_genomes
│   │   │   └── small
│   │   │       └── small_reference.*
│   │   └── test_data                     # Test data
│   │       ├── coordtable.tsv
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
│   ├── container.config                  # Container registry addresses for computational tools.
│   ├── ngc.config                        # Configuration profiles the NGC-HPC compute enviornment.
│   └── profiles.config                   # Configuration profile for default environment.
├── envs
│   └── <name>/
│       └── environment.yaml              # Conda environment definitions
├── main.nf                               # Main workflow
├── modules/
│   └── <module>.nf                       # Module scripts
├── subworkflows/
│   └── <subworkflow>.nf                  # Module scripts
├── next.pbs                              # Helper script for running on NGC-HPC
└── nextflow.config                       # Workflow parameters
```


## Usage on NGC-HPC

### 0. Create qsub helper script

Create a wrapper script for qsub, so you don't have to keep track of working directory, group etc. again.
First do `mkdir ~/bin`. Then save the following script as a file named `~/bin/myqsub` and make it executable by `chmod +x ~/bin/myqsub`.

```
#!/bin/bash

qsub -W group_list=icope_staging_r -A icope_staging_r -d $(pwd) "$@"
```

Add ~/bin to your path. You can have this done on log-in by appending the following line to your `~/.bashrc`:

```
export PATH="$PATH:$HOME/bin"
```

### 1. In a clean folder - clone this repository from the local NGC-HPC folder:

```Bash
git clone /ngc/projects/icope_staging_r/git/predisposed/.git .
```

### 2. Pipeline preview with test data:

```Bash
bash next.pbs -params-file test_config.json -stub
```

### 3. Run pipeline with test data locally (all heavy jobs are still submitted with qsub):

```Bash
bash next.pbs -params-file test_config.json
```

### 4. Submit pipeline with test data:

```Bash
myqsub next.pbs -F "-params-file test_config.json"
```

## Running the pipeline with _real_ data

### 1. Organization
While the pipeline is still under development, it make sense to create new clones for each pipeline run, to keep track of possible changes done while running it. I propose this folder structure:

```Bash
predisposed
├── git/DoBSeqWF                            # Temporary local workflow repository
├── resources                               # Reference genome and target files.
├── data/                                   # Raw data for each batch
│   ├── <batch_id_I>/
│   │   └── *.fq.gz
│   ├── <batch_id_II>/
│   │   └── *.fq.gz
│   └── <batch_id_III>/
│       └── *.fq.gz
│   └── ...
└── processed_data/                         # Processed data for each batch
    ├── <batch_id_I>/
    │   ├── DoBSeqWF/                       # Clone repository here
    │   │   ├── config.json                 # Configuration file
    │   │   ├── pooltable.tsv               # Pool table (create with helper script)
    │   │   └── decodetable.tsv             # Decode table (we need a convention for this)
    │   └── results
    │       ├── cram/                       # CRAM files for each pool
    │       ├── logs/                       # Log files for each process
    │       ├── variants/                   # VCF files for each pool
    │       ├── variant_tables/             # TSV files converted from pool VCFs
    │       └── pinpoint_variants/
    │           ├── all_pins/               # All pinpointables for each sample in individual vcfs (*note)
    │           ├── unique_pins/            # All unique pinpointables for each sample in individual vcfs (*note)
    │           ├── *_merged.vcf.gz         # All pinpointables for all samples in a single vcf without sample information
    │           ├── summary.tsv             # Variant counts for each sample
    │           └── lookup.tsv              # Variant to sample lookup table
    ├── <batch_id_II>/
    │   ├── DoBSeqWF/
    │   └── results/
    ├── <batch_id_III>/
    │   ├── DoBSeqWF/
        └── results/
    └── ...
```
(*note) Each pinpointable variant can be represented by the horizontal or the vertical pools. In order not to loose any information, there are, _for now_, 6 vcf files for each sample. Four with representations from either dimension named {sample}\_{pool}\_{unique/all}\_pins.vcf.gz and 2 with all pins merged named {sample}\_{unique/all}.vcf.gz.

### 2. Prepare folders

```Bash
cd /ngc/projects2/dp_00005/data/predisposed/
mkdir -p data/<batch_id> processed_data/<batch_id>
mv /ssi/fastq/data /ngc/projects2/dp_00005/data/predisposed/data/<batch_id>/
```

### 2. Clone pipeline repository

```Bash
cd processed_data/<batch_id>
git clone /ngc/projects2/dp_00005/data/predisposed/git/DoBSeqWF
```

### 3. Create pooltable.tsv:

```Bash
cd DoBSeqWF
bash assets/helper_scripts/create_pooltable.sh ../../../data/<batch_id>/
```

### 4. Adjust configuration file:

Fill out config.json with the correct paths and parameters. Decodetable is not needed for mapping only. Look into nextflow.config for possible parameters to set in the conifg.json.

### 5. Run pipeline

```Bash
myqsub next.pbs -F "-params-file config.json"
```

### 6. Monitor progress

```Bash
tail nextflow.log
```

If the pipeline fails - it is likely due to resource constraints. Adjust as needed in the conf/profiles.config file under NGC, and rerun the PBS script. Be aware that any direct edits of the workflow scripts, ie. modules and subworkflows, can lead to complete re-run of the pipeline.

