#!/usr/bin/env bash
set -euo pipefail

# Example benchmark run for peptide ACDE with Amber14
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

python analysis.py test_run/ACDE_amber14.pdb \
  --output_prefix test_run/ACDE_amber14 \
  --hb-freq 0.1
