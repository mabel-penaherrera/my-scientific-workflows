"""
Simulation driver for short peptide molecular dynamics benchmarking.

Given a peptide sequence and a force field, this builds a capped peptide, solvates it,
optionally neutralizes it with ions, and runs a short MD trajectory with energy
minimization and dynamics. Outputs:
  - Minimized structure (PDB)
  - Dynamic trajectory (PDB)
  - Console status/logging
  
"""

import argparse
import os
import sys
import tempfile
import warnings

from openmm import unit, Platform
from openmm.app import (
    ForceField,
    Modeller,
    Simulation,
    PDBReporter,
    StateDataReporter,
    PDBFile,
)
from openmm import LangevinIntegrator
import openmm.app as app  # for PME, HBonds

from pdbfixer import PDBFixer
import PeptideBuilder
from Bio.PDB import PDBIO


def build_peptide_structure(sequence: str, phi=None, psi_im1=None):
    """
    Build a capped peptide structure from a FASTA 1-letter sequence using PeptideBuilder.

    Parameters
    sequence : str
        Peptide sequence (e.g., "ACDE").
    phi : list or None
        Sequence of phi angles in degrees for residues (len = len(sequence)-1).
    psi_im1 : list or None
        Sequence of psi_{i-1} angles in degrees. Same length as phi.

    Returns
    structure : Bio.PDB.Structure.Structure
        The peptide structure with ACE and NME caps.
    """
    if len(sequence) < 1:
        raise ValueError("Sequence must have at least one residue.")

    # Default backbone angles: extended/alpha-like (tunable)
    if phi is None:
        phi = [-120] * (len(sequence) - 1)
    if psi_im1 is None:
        psi_im1 = [140] * (len(sequence) - 1)

    structure = PeptideBuilder.make_structure(sequence, phi, psi_im1)
    return structure


def fasta_to_fixed_modeller(
    sequence: str,
    forcefield_obj: ForceField,
    padding: unit.Quantity = 1.0 * unit.nanometer,
    neutralize: bool = True,
    ionic_strength: unit.Quantity = 0.0 * unit.molar,
    water_model: str = "tip3p",
    pH: float = 7.0
) -> Modeller:
    """
    Prepares a solvated peptide system for molecular dynamics simulations.
    
    Standard workflow:
    1. Build peptide structure from sequence
    2. Fix missing atoms/protonation states at target pH
    3. Solvate with appropriate water model
    4. Neutralize charge with ions
    
    Parameters:
        sequence: Peptide sequence (1-letter FASTA code)
        forcefield_obj: Initialized OpenMM ForceField object
        padding: Solvent padding distance (default: 1.0 nm)
        neutralize: Add counterions to neutralize system charge (default: True)
        ionic_strength: Additional ion concentration (default: 0.0 M)
        water_model: Water model for solvation (default: "tip3p")
        pH: Target pH for protonation state (default: 7.0)
    
    Returns:
        Prepared OpenMM Modeller object with solvated system
    
    Raises:
        ValueError: For invalid input parameters
        TypeError: For incorrect unit types
        RuntimeError: For failures in building, fixing, or solvation steps
    """
    # Validate input parameters
    # Sequence validation
    if not sequence or not isinstance(sequence, str):
        raise ValueError("Peptide sequence must be a non-empty string")
    
    # Unit type validation
    if not isinstance(padding, unit.Quantity):
        raise TypeError("Padding must be a unit.Quantity")
    if not isinstance(ionic_strength, unit.Quantity):
        raise TypeError("Ionic strength must be a unit.Quantity")
    
    # Scientific parameter validation
    if padding.value_in_unit(unit.nanometer) <= 0:
        raise ValueError("Padding must be positive")
    if ionic_strength.value_in_unit(unit.molar) < 0:
        raise ValueError("Ionic strength cannot be negative")
    if pH < 0 or pH > 14:
        raise ValueError("pH must be between 0-14")
    
    # Water model compatibility check
    valid_water_models = ["tip3p", "tip4p", "spce", "tip4pew", "tip5p"]
    if water_model.lower() not in valid_water_models:
        warnings.warn(f"Unrecognized water model: {water_model}. "
                      f"Supported models: {', '.join(valid_water_models)}")

    tmp_name = None
    try:
        # Build peptide structure from sequence
        structure = build_peptide_structure(sequence)
        
        # Create temporary file for PDBFixer processing
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
            tmp_name = tmp.name
            io = PDBIO()
            io.set_structure(structure)
            io.save(tmp_name)
        
        # Fix structural issues and add hydrogens at target pH
        fixer = PDBFixer(filename=tmp_name)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=pH)
        
    except Exception as e:
        raise RuntimeError(f"Peptide preparation failed: {e}") from e
    finally:
        # Guaranteed cleanup of temporary file
        if tmp_name and os.path.exists(tmp_name):
            try:
                os.remove(tmp_name)
            except OSError:
                pass  # Best-effort cleanup

    try:
        # Create solvated system with ions
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.addSolvent(
            forcefield_obj,
            model=water_model,
            padding=padding,
            neutralize=neutralize,
            ionicStrength=ionic_strength,
        )
        
        # Ensure solvation actually created atoms
        if modeller.topology.getNumAtoms() == 0:
            raise RuntimeError("Solvation failed: resulted in empty topology")
        
        return modeller
        
    except Exception as e:
        raise RuntimeError(f"Solvation failed: {e}") from e

def create_system_and_simulation(
    modeller,
    forcefield_obj,
    temperature=300 * unit.kelvin,
    friction=1 / unit.picosecond,
    timestep=0.002 * unit.picoseconds,
    platform_name=None,
    output_pdb="output.pdb",
    steps=1000,
    report_interval=100,
    seed=None,
):
    """
    Given a prepared Modeller and force field, create the system, integrate, minimize, and run.

    Outputs both a minimized structure and a dynamic trajectory.
    """
    # System creation with standard benchmarking options
    system = forcefield_obj.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
    )

    # Integrator
    integrator = LangevinIntegrator(temperature, friction, timestep)
    if seed is not None:
        integrator.setRandomNumberSeed(seed)

    # Platform selection with fallback
    simulation = None
    if platform_name:
        try:
            platform = Platform.getPlatformByName(platform_name)
            simulation = Simulation(modeller.topology, system, integrator, platform)
        except Exception as e:
            warnings.warn(f"Failed to load platform '{platform_name}': {e}. Falling back to default.")
    if simulation is None:
        simulation = Simulation(modeller.topology, system, integrator)

    simulation.context.setPositions(modeller.positions)

    # Energy minimization
    print("[INFO] Minimizing energy...")
    simulation.minimizeEnergy()

    # Ensure output directory exists
    out_dir = os.path.dirname(output_pdb)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Write minimized structure separately for transparency
    min_path = output_pdb.replace(".pdb", "_minimized.pdb")
    print(f"[INFO] Writing minimized structure to {min_path}...")
    state = simulation.context.getState(getPositions=True)
    with open(min_path, "w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    # Reporters: trajectory + progress info
    simulation.reporters.append(PDBReporter(output_pdb, report_interval))
    simulation.reporters.append(
        StateDataReporter(
            sys.stdout,
            report_interval,
            step=True,
            temperature=True,
            potentialEnergy=True,
            density=True,
            progress=True,
            remainingTime=True,
            totalSteps=steps,
        )
    )

    print(f"[INFO] Running dynamics for {steps} steps...")
    simulation.step(steps)
    print(f"[INFO] Trajectory written to {output_pdb}")


def resolve_forcefield(name_or_path):
    """
    Normalize a scalar force field identifier to actual XMLs.
    """
    if name_or_path.lower() == "amber14":
        return ForceField("amber14-all.xml", "amber14/tip3p.xml")
    elif name_or_path.lower() == "charmm36":
        return ForceField("charmm36.xml", "charmm36/water.xml")
    else:
        return ForceField(name_or_path)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run short peptide MD with specified force field for benchmarking."
    )
    parser.add_argument(
        "--sequence",
        type=str,
        required=True,
        help="Peptide sequence (1-letter FASTA, e.g., ACDE).",
    )
    parser.add_argument(
        "--forcefield",
        type=str,
        required=True,
        help="Force field key: 'amber14', 'charmm36', or path to custom OpenMM XML.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output.pdb",
        help="Output trajectory (PDB) filename.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=1000,
        help="Number of integration steps to take after minimization.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Temperature in Kelvin.",
    )
    parser.add_argument(
        "--friction",
        type=float,
        default=1.0,
        help="Friction coefficient in 1/ps.",
    )
    parser.add_argument(
        "--timestep",
        type=float,
        default=2.0,
        help="Timestep in femtoseconds.",
    )
    parser.add_argument(
        "--platform",
        type=str,
        default=None,
        help="OpenMM platform to use (e.g., CUDA, CPU). Defaults to automatic.",
    )
    parser.add_argument(
        "--padding",
        type=float,
        default=1.0,
        help="Solvent padding in nanometers.",
    )
    parser.add_argument(
        "--no-neutralize",
        action="store_true",
        help="Do not add counterions to neutralize system charge.",
    )
    parser.add_argument(
        "--ionic-strength",
        type=float,
        default=0.0,
        help="Additional salt concentration in molar (e.g., 0.15 for 150 mM).",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for Langevin integrator."
    )
    parser.add_argument(
        "--water-model",
        type=str,
        default="tip3p",
        help="Water model to pass to addSolvent (e.g., 'tip3p').",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    ff = resolve_forcefield(args.forcefield)
    padding = args.padding * unit.nanometer
    neutralize = not args.no_neutralize
    ionic_strength = args.ionic_strength * unit.molar

    modeller = fasta_to_fixed_modeller(
        args.sequence,
        ff,
        padding=padding,
        neutralize=neutralize,
        ionic_strength=ionic_strength,
        water_model=args.water_model,
    )

    if modeller is None:
        raise RuntimeError("Modeller creation failed; got None. Check previous error messages.")

    create_system_and_simulation(
        modeller,
        ff,
        temperature=args.temperature * unit.kelvin,
        friction=(args.friction / unit.picosecond),
        timestep=args.timestep * unit.femtoseconds,
        platform_name=args.platform,
        output_pdb=args.output,
        steps=args.steps,
        report_interval=max(1, args.steps // 10),
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
