"""
Trajectory analysis for peptide simulations.

Computes:
  - Full-atom and backbone RMSD (reference = first frame)
  - Backbone torsion angles (phi/psi) in degrees
  - Ramachandran plot (if enough data)
  - Hydrogen bond summary (best-effort, adjustable frequency)

Outputs PNGs and CSVs with a given prefix. Warns if data are too sparse.
"""

import argparse
import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import sys


def ensure_dir(path):
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def circular_mean_deg(angles_deg):
    """Compute circular mean of angles given in degrees, returns degrees."""
    angles_rad = np.deg2rad(angles_deg)
    sin_sum = np.sum(np.sin(angles_rad), axis=0)
    cos_sum = np.sum(np.cos(angles_rad), axis=0)
    mean_angle = np.arctan2(sin_sum, cos_sum)
    return np.rad2deg(mean_angle)  # in degrees


def analyze_trajectory(pdb_file, traj_file, output_prefix, hb_freq):
    print(f"[INFO] Starting analysis on {traj_file or pdb_file}")
    if traj_file is None:
        traj = md.load(pdb_file)
    else:
        traj = md.load(traj_file, top=pdb_file)

    print(f"[INFO] Loaded trajectory with {traj.n_frames} frame(s).")
    if traj.n_frames < 1:
        raise RuntimeError("Trajectory has no frames.")

    # Ensure output directory exists
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        ensure_dir(out_dir)

    # Full-atom RMSD
    print("[INFO] Computing full-atom RMSD...")
    rmsd_full = md.rmsd(traj, traj, 0)  # nanometers
    rmsd_full_A = rmsd_full * 10.0  # angstroms

    np.savetxt(
        f"{output_prefix}_rmsd_full_nm.csv",
        rmsd_full,
        delimiter=",",
        header="RMSD_full (nm)",
        comments="",
    )
    np.savetxt(
        f"{output_prefix}_rmsd_full_A.csv",
        rmsd_full_A,
        delimiter=",",
        header="RMSD_full (Å)",
        comments="",
    )

    # Backbone RMSD
    print("[INFO] Computing backbone RMSD...")
    rmsd_bb_A = None
    try:
        backbone_sel = traj.topology.select("backbone")  # N, CA, C
        if len(backbone_sel) == 0:
            raise ValueError("No backbone atoms found for selection.")
        traj_bb = traj.atom_slice(backbone_sel)
        rmsd_bb = md.rmsd(traj_bb, traj_bb, 0)
        rmsd_bb_A = rmsd_bb * 10.0
        np.savetxt(
            f"{output_prefix}_rmsd_bb_nm.csv",
            rmsd_bb,
            delimiter=",",
            header="RMSD_backbone (nm)",
            comments="",
        )
        np.savetxt(
            f"{output_prefix}_rmsd_bb_A.csv",
            rmsd_bb_A,
            delimiter=",",
            header="RMSD_backbone (Å)",
            comments="",
        )
    except Exception as e:
        print(f"[WARNING] Backbone RMSD failed: {e}")

    # Plot RMSD comparison
    plt.figure()
    frames = np.arange(traj.n_frames)
    plt.plot(frames, rmsd_full_A, marker="o", linestyle="-", label="full atom")
    if rmsd_bb_A is not None:
        plt.plot(frames, rmsd_bb_A, marker="s", linestyle="--", label="backbone")
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD vs Frame (ref = frame 0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_rmsd_comparison.png")
    plt.close()

    # Individual RMSD plots
    plt.figure()
    plt.plot(frames, rmsd_full_A, marker="o", linestyle="-")
    plt.xlabel("Frame")
    plt.ylabel("RMSD (Å)")
    plt.title("Full-atom RMSD vs Frame (ref = frame 0)")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_rmsd_full.png")
    plt.close()

    if rmsd_bb_A is not None:
        plt.figure()
        plt.plot(frames, rmsd_bb_A, marker="s", linestyle="--")
        plt.xlabel("Frame")
        plt.ylabel("RMSD (Å)")
        plt.title("Backbone RMSD vs Frame (ref = frame 0)")
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_rmsd_backbone.png")
        plt.close()

    # If only one frame, skip torsion and Ramachandran
    if traj.n_frames < 2:
        print("[WARNING] Only one frame present; skipping torsion and Ramachandran analysis.")
        return

    # Phi / Psi torsions
    print("[INFO] Computing phi/psi torsions...")
    phi_indices, phi_angles = md.compute_phi(traj)
    psi_indices, psi_angles = md.compute_psi(traj)

    phi_deg = np.degrees(phi_angles)
    psi_deg = np.degrees(psi_angles)

    phi_headers = [f"phi_res{i+1}" for i in range(phi_deg.shape[1])]
    psi_headers = [f"psi_res{i+1}" for i in range(psi_deg.shape[1])]

    np.savetxt(
        f"{output_prefix}_phi_deg.csv",
        phi_deg,
        delimiter=",",
        header=",".join(phi_headers),
        comments="",
    )
    np.savetxt(
        f"{output_prefix}_psi_deg.csv",
        psi_deg,
        delimiter=",",
        header=",".join(psi_headers),
        comments="",
    )

    # Summary statistics of torsions (circular mean)
    phi_mean = circular_mean_deg(phi_deg)
    psi_mean = circular_mean_deg(psi_deg)
    with open(f"{output_prefix}_torsion_summary.txt", "w") as f:
        f.write("Circular mean of backbone torsions (deg):\n")
        for i in range(len(phi_mean)):
            f.write(f"Residue {i+1}: phi={phi_mean[i]:.2f}, psi={psi_mean[i]:.2f}\n")

    # Ramachandran
    if phi_deg.size > 0 and psi_deg.size > 0:
        print("[INFO] Creating Ramachandran plot...")
        plt.figure()
        plt.scatter(
            phi_deg.flatten(),
            psi_deg.flatten(),
            s=20,
            alpha=0.6,
            edgecolors="none",
        )
        plt.xlabel("Phi (deg)")
        plt.ylabel("Psi (deg)")
        plt.title("Ramachandran Plot")
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.grid(True, linestyle=":", linewidth=0.5)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_ramachandran.png")
        plt.close()
    else:
        print("[WARNING] Insufficient torsion data for Ramachandran plot.")

    # Hydrogen bonds
    print("[INFO] Computing hydrogen bonds (Baker-Hubbard)...")
    try:
        hbonds = md.baker_hubbard(traj, freq=hb_freq)
        with open(f"{output_prefix}_hbonds.txt", "w") as f:
            if len(hbonds) > 0:
                f.write("Detected hydrogen bonds (donor - hydrogen - acceptor):\n")
                for h in hbonds:
                    donor, hydrogen, acceptor = h
                    f.write(
                        f"{traj.topology.atom(donor)} - {traj.topology.atom(hydrogen)} - {traj.topology.atom(acceptor)}\n"
                    )
                print(f"[INFO] Wrote {len(hbonds)} hydrogen bonds.")
            else:
                f.write("No hydrogen bonds detected.\n")
                print("[INFO] No hydrogen bonds detected.")
    except Exception as e:
        print(f"[WARNING] Hydrogen bond calculation failed: {e}")
        with open(f"{output_prefix}_hbonds.txt", "w") as f:
            f.write(f"Error during hydrogen bond analysis: {e}\n")

    print("[INFO] Analysis complete.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze peptide trajectory: RMSD, phi/psi, Ramachandran, H-bonds."
    )
    parser.add_argument(
        "pdb",
        type=str,
        help="Reference PDB/topology file (often the output from simulation).",
    )
    parser.add_argument(
        "--trajectory",
        type=str,
        default=None,
        help="Optional separate trajectory file (e.g., DCD, XTC). If omitted, uses PDB.",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Prefix for output files (e.g., test_run/ACDE_amber14).",
    )
    parser.add_argument(
        "--hb-freq",
        type=float,
        default=0.1,
        help="Frequency parameter for Baker-Hubbard hydrogen bond detection.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    print("DEBUG: analysis2.py started")
    analyze_trajectory(args.pdb, args.trajectory, args.output_prefix, args.hb_freq)


if __name__ == "__main__":
    main()