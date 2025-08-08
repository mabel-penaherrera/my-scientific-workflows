# Scientific Workflows Portfolio

This repository showcases two research projects involving automated scientific workflows, simulation/analysis pipelines, and data provenance.

## Projects

### 1. High-Mass X-ray Binary Candidate Identification and Cross-matching
**Research artifact.** Automated cross-matching and probabilistic scoring of multiwavelength catalogs to identify high-mass X-ray binary (HMXB) candidates near the Galactic Center. This work was developed as part of a published astrophysics study; the notebook reflects project-specific data, assumptions, and catalogs.
[→ See Project](./high_mass_xray_binary_pipeline/)

### 2. Biomolecular Force Field Benchmarking  
Reusable, GPU-accelerated molecular dynamics simulation and analysis of short peptide conformational behavior across molecular mechanics force fields (e.g., Amber14) using OpenMM. Produces minimized vs dynamic structures, RMSD (full-atom and backbone), torsion statistics, Ramachandran plots, and hydrogen-bond summaries.  
[→ See Project](./biomolecular_sim_benchmark/)

## Repository Layout

```
.
├── README.md # top-level README (this file)
├── high_mass_xray_binary_pipeline/
│   ├── hmxb_crossmatching.ipynb     # annotated Jupyter notebook (primary analysis)
│   ├── hmxb_crossmatching.pdf       # static export of the notebook for easy viewing
│   ├── hmxb_crossmatching.html       # HTML-rendered notebook (optional for previews)
│   ├── README.md                   # HMXB-specific README
│   ├── requirements.txt           # dependencies for HMXB analysis
│   └── example_data/              # minimal / synthetic data for testing
│   └── outputs/              # generated figures from the notebook
│   └── LICENSE
│   └── GET_DATA.md # instructions to fetch real input catalogs
└── biomolecular_sim_benchmark/
    ├── run_simulations.py          # simulation driver
    ├── analysis.py                # trajectory analysis
    ├── example.sh                # canonical example wrapper
    ├── README.md                 # benchmark-specific README
    ├── launch_parallel.py       # batch launcher for multiple configs
    ├── batch_config.json        # example batch config file
    ├── requirements.txt          # dependencies for simulation+analysis
    ├── test_run/                 # example outputs   
    └── LICENSE 
  
```

## Getting Started

Each subproject is self-contained. Navigate into the relevant directory and follow its README for detailed setup and execution instructions.

Example:
```bash
# Biomolecular pipeline
cd biomolecular_sim_benchmark
./example.sh

# HMXB pipeline
cd ../high_mass_xray_binary_pipeline
# open the notebook or follow instructions in that README
```

## Notes

- The HMXB project is tied to a specific published study; its inputs and assumptions are research-specific. Refer to its README for citation, data requirements, and context.  
- The biomolecular pipeline is designed for broader reuse and benchmarking across sequences and force fields, with built-in reproducibility (seeded dynamics).

## Contact

Open an issue on GitHub or reach out via the contact method provided in the individual project READMEs.
