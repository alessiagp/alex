#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import argparse
import sys


def traj_pruning(sliced_part_traj, desired_N_frames):
    """
    Reduce the number of frames in an MDtraj trajectory to exactly desired_N_frames.
    """
    tot_frames = len(sliced_part_traj)
    print(f"Total frames in the current trajectory: {tot_frames}")

    if tot_frames > desired_N_frames:
        frame_ndx = np.round(
            np.linspace(0, tot_frames - 1, desired_N_frames)
        ).astype(int)

        # Remove duplicates
        frame_ndx = np.unique(frame_ndx)

        # Ensure exact number of frames
        while len(frame_ndx) < desired_N_frames:
            frame_ndx = np.append(frame_ndx, tot_frames - 1)
            frame_ndx = np.unique(frame_ndx)

        reduced_traj = sliced_part_traj[frame_ndx]
        print(f"Reduced trajectory length: {len(reduced_traj)}")
        return reduced_traj
    else:
        print(f"Trajectory already has {tot_frames} frames (≤ {desired_N_frames})")
        return sliced_part_traj


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prune an MD trajectory to a fixed number of frames."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input trajectory file (e.g. .xtc, .dcd)"
    )

    parser.add_argument(
        "-t", "--topology",
        required=True,
        help="Topology file (e.g. .pdb, .gro)"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output trajectory file"
    )

    parser.add_argument(
        "-n", "--nframes",
        type=int,
        required=True,
        help="Desired number of frames"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    print("Loading trajectory...")
    try:
        traj = md.load(args.input, top=args.topology)
    except Exception as e:
        print(f"Error loading trajectory: {e}")
        sys.exit(1)

    pruned_traj = traj_pruning(traj, args.nframes)

    print("Saving pruned trajectory...")
    try:
        pruned_traj.save(args.output)
    except Exception as e:
        print(f"Error saving trajectory: {e}")
        sys.exit(1)

    print("Done.")


if __name__ == "__main__":
    main()