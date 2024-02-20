# DoBSeqWF

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
│       └── simulator.py
├── bin
│   └── <script>.*
├── conf
│   └── profiles.config
├── envs
│   └── <name>/
│       └── environment.yaml
├── main.nf
├── modules/
│   └── <module>.nf
└── nextflow.config