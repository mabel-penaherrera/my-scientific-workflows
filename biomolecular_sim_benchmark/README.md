# Peptide Force Field Benchmarking Pipeline

This project benchmarks classical biomolecular force fields by simulating short peptides in water using [OpenMM](http://openmm.org/). It builds capped peptide structures, runs explicit-solvent Langevin dynamics, and quantifies backbone conformational behavior through RMSD, torsion statistics (ϕ/ψ), Ramachandran sampling, and hydrogen bonding.

---

## Features

- Constructs capped peptide from a 1-letter FASTA sequence.  
- Supports multiple force fields (`amber14`, `charmm36`, or custom OpenMM XML).  
- Solvates and neutralizes, with optional ionic strength.  
- Performs energy minimization and Langevin dynamics with user-configurable seed, temperature, friction, timestep, and length.  
- Batch mode via JSON + parallel launcher.
- Analysis includes:
  - Full-atom vs backbone RMSD with comparative plots.  
  - Backbone torsion angles (ϕ/ψ) in degrees and circular means.  
  - Ramachandran sampling visualization.  
  - Hydrogen bond detection (Baker–Hubbard) with adjustable sensitivity.  
- One-command example (`example.sh`) for end-to-end reproduction.

---

## Quickstart

### 1. Install dependencies (example via conda)
```bash
conda create -n biosim python=3.9 -y
conda activate biosim
conda install -c conda-forge openmm=8.0 mdtraj=1.9.7 numpy matplotlib biopython pdbfixer peptidebuilder
```

### 2. Run the canonical example (peptide `ACDE` with Amber14)
```bash
# Simulation
mkdir -p test_run
python run_simulations.py \
  --sequence ACDE \
  --forcefield amber14 \
  --output test_run/ACDE_amber14.pdb \
  --steps 5000 \
  --temperature 300 \
  --friction 1 \
  --timestep 2 \
  --seed 42
  
---

## Parallel Batch Mode

You can run multiple simulations in parallel using the provided JSON config and launcher script.

### 1. Edit `batch_config.json`
```json
[
  {
    "sequence": "ACDE",
    "forcefield": "amber14",
    "output": "outputs/ACDE_amber14.pdb"
  },
  {
    "sequence": "ACDE",
    "forcefield": "charmm36",
    "output": "outputs/ACDE_charmm36.pdb"
  },
  {
    "sequence": "GGR",
    "forcefield": "amber14",
    "output": "outputs/GGR_amber14.pdb"
  }
]

### 1. Launch all jobs in parallel
python launch_parallel.py

---

## Analysis
python analysis.py test_run/ACDE_amber14.pdb \
  --output_prefix test_run/ACDE_amber14 \
  --hb-freq 0.1
```

Or simply:
```bash
./example.sh
```
---

## Example Outputs (from the above run)

- `ACDE_amber14_minimized.pdb` — minimized starting geometry.  
- `ACDE_amber14.pdb` — production trajectory (multiple MODEL frames).  
- `ACDE_amber14_rmsd_comparison.png` — full-atom vs backbone RMSD.  
- `ACDE_amber14_rmsd_full.png` / `ACDE_amber14_rmsd_backbone.png` — individual RMSD traces.  
- `ACDE_amber14_phi_deg.csv` & `ACDE_amber14_psi_deg.csv` — per-frame backbone torsions in degrees.  
- `ACDE_amber14_torsion_summary.txt` — circular means of φ/ψ.  
- `ACDE_amber14_ramachandran.png` — conformational sampling.  
- `ACDE_amber14_hbonds.txt` — hydrogen bond summary.

---

## Interpretation Notes

- **Full-atom vs backbone RMSD:** Backbone RMSD isolates peptide backbone motion (~1 Å), while full-atom drift reflects broader structural relaxation.  
- **Torsion circular means:** Provide compact fingerprints of average backbone geometry.  
- **Ramachandran plot:** Shows allowed sampling regions (example run samples extended-like conformations).  
- **Hydrogen bonding:** Short flexible peptides may lack persistent H-bonds; longer simulations or different sequences can reveal them.

---

## Extending / Benchmarking

- Change sequence: `--sequence GGGG` or any other 1-letter peptide.  
- Swap force field: `--forcefield charmm36` or supply custom XML.  
- Modify solvent conditions: `--ionic-strength`, `--no-neutralize`.  
- Ensure reproducibility: use `--seed` for deterministic Langevin dynamics.  
- Compare runs: vary seeds or force fields and aggregate torsion / RMSD statistics.
- Launch many jobs: populate batch_config.json and run launch_parallel.py.

---

## Installation / Dependencies

Suggested `requirements.txt`:
```
openmm==8.0
mdtraj==1.9.7
numpy
matplotlib
biopython
pdbfixer
PeptideBuilder
```

---

## Reproducibility

Recommended to capture run metadata per execution (e.g., in `run_metadata.json`):
```json
{
  "sequence": "ACDE",
  "forcefield": "amber14",
  "steps": 5000,
  "temperature": 300,
  "friction": 1.0,
  "timestep_fs": 2.0,
  "seed": 42,
  "water_model": "tip3p"
}
```
---

## Troubleshooting

- **Only one frame**: Increase `--steps` or adjust reporter frequency.  
- **Sparse Ramachandran**: Need multiple frames and sufficient residue count.  
- **No hydrogen bonds**: Short peptides may not form stable H-bonds; extend sampling or tune `--hb-freq`.  
- **Force field errors**: Ensure XMLs are correctly available or paths are valid.

---

## Example alternative usage

```bash
python run_simulations.py --sequence GGGG --forcefield charmm36 \
  --output test_run/GGGG_charmm36.pdb --steps 20000 \
  --temperature 300 --friction 1 --timestep 2 --seed 123

python analysis.py test_run/GGGG_charmm36.pdb --output_prefix test_run/GGGG_charmm36
```

---

##  References

- [OpenMM Documentation](http://docs.openmm.org/)
- [PeptideBuilder GitHub](https://github.com/mtien/PeptideBuilder)
- [PDBFixer](https://github.com/openmm/pdbfixer)

---

## License

This project is licensed under the **MIT License**. See `LICENSE` for details.

---

## Contact

Open an issue on GitHub or reach out at andrea.mabel13 [at] gmail [dot] com