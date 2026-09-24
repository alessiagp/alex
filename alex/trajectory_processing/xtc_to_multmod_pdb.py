#!/usr/bin/env python3
import argparse
import sys
import os
import MDAnalysis as mda

def main():
    parser = argparse.ArgumentParser(description="Convert MD trajectories to multi-model PDBs for RING analysis.")
    parser.add_argument("-p", "--top", required=True, help="Input topology (.pdb, .gro, .tpr)")
    parser.add_argument("-x", "--traj", required=True, help="Input trajectory (.xtc, .dcd)")
    parser.add_argument("-o", "--out", required=True, help="Output multi-model PDB file")
    parser.add_argument("-s", "--stride", type=int, default=10, help="Frame stride (default: 10)")
    parser.add_argument("--sel", type=str, default="protein", help="Atom selection (default: 'protein')")
    
    args = parser.parse_args()

    if not os.path.isfile(args.top) or not os.path.isfile(args.traj):
        sys.exit(f"Error: Could not find {args.top} or {args.traj}")

    print(f"Loading universe...")
    u = mda.Universe(args.top, args.traj)
    
    selection = u.select_atoms(args.sel)
    if len(selection) == 0:
        sys.exit(f"Error: The selection '{args.sel}' contains 0 atoms.")

    total_frames = len(u.trajectory)
    processed_frames = len(range(0, total_frames, args.stride))
    print(f"Extracting {processed_frames} frames (stride {args.stride}) out of {total_frames} total frames...")

    # Write the multiframe PDB
    with mda.Writer(args.out, multiframe=True) as W:
        for ts in u.trajectory[::args.stride]:
            W.write(selection)
            
    print(f"Success! Saved to {args.out}")

if __name__ == "__main__":
    main()