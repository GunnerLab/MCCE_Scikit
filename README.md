# MCCE_Scikit (alpha)
## [MCCE](https://github.com/GunnerLab/Stable-MCCE)'s scientific tookit to perform analyses on its output files.

# Current project structure:
```
.
├── LICENSE
├── README.md
├── msa_minimal_envir.yaml
├── notebooks
│   └── ms_sampling.ipynb
├── pyproject.toml
├── src
│   ├── VERSION
│   └── mcceSK
│       ├── __init__.py
│       ├── base.py
│       ├── constants.py
│       ├── mcce_io.py
│       └── ms_sampling.py
└── tests
    └── data
        ├── head3.lst
        ├── ms_out
        │   ├── pH5eH0ms
        │   │   #├── MCi files not uploaded
        │   │   ├── header
        │   │   └── pdbs_from_ms
        │   │    # Example: i=MONTERUN index, j=selected ms
        │   │       ├── mc{i}_ms{j}.pdb
        │   │       ├── [..]
        │   │       └── mc0_ms91192.pdb
        │   ├── pH5eH0ms.txt
        │   └── pH6eH0ms.txt
        ├── run.prm.record
        └── step2_out.pdb
```
