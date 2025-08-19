import numpy as np
import logging
from pathlib import Path
from collections import Counter
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import os
import sys
import MDAnalysis as mda
import random

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class KLProbabilities:
    def __init__(self, dist_mat: str, gro_file: str, traj_file: str, save_dir: str,  struct_type: str, struct: str, L: int):
        self.dist_mat = Path(dist_mat)
        self.gro_file = Path(gro_file)
        self.traj_file = Path(traj_file)
        self.struct = struct
        self.L = L
        
        self.save_dir = Path(save_dir)
        self.struct_type = struct_type
        self.save_dir.mkdir(parents=True, exist_ok=True)

        self.t = mda.Universe(str(self.gro_file), str(self.traj_file))
        self.dist_matrix = np.load(self.dist_mat)

    def clustering(self, mat1: np.ndarray) -> np.ndarray:
        """
        Perform clustering on the RMSD matrix.
        """
        if mat1.shape[0] != mat1.shape[1]:
            raise ValueError("Distance matrix must be square for clustering.")

        logging.info("Starting clustering.")
        Z1 = linkage(squareform(mat1), method='average')
        cl_labels = fcluster(Z1, t=self.L, criterion='maxclust')
        logging.info("Clustering completed.")
        return cl_labels

    def labeling(self, cl_labels: np.ndarray) -> np.ndarray:
        """
        Compute microstate probabilities based on clustering labels.
        Select representative frames and write them to a new trajectory.
        """
        label_counts = Counter(cl_labels)
        total = sum(label_counts.values())

        # Select one random frame per cluster
        unique_classes = list(label_counts.keys())
        cl_selection = np.array([random.choice(np.where(cl_labels == class_label)[0]) for class_label in unique_classes])

        logging.info(f"Selected frames: {cl_selection}")
        logging.info(f"Number of selected frames: {len(cl_selection)}")

        # Compute label probabilities AFTER selecting frames
        label_probs = np.array([label_counts[class_label] / total for class_label in unique_classes])

        logging.info(f"Length of label_probs: {len(label_probs)}")
        logging.info("Label probabilities computed.")

        # Output file path
        output_traj = self.save_dir / f"{self.struct}-{self.struct_type}-KL_frames.xtc"

        # Ensure there are atoms before writing
        if len(self.t.atoms) == 0:
            raise ValueError("No atoms selected, cannot write trajectory.")

        # Write selected frames to a new trajectory
        with mda.Writer(str(output_traj), n_atoms=self.t.atoms.n_atoms) as writer:
            for frame in cl_selection:
                self.t.trajectory[frame]
                writer.write(self.t)

        logging.info(f"Selected frames saved to: {output_traj}")

        return label_probs

    def savefile(self, microstates_probs: np.ndarray):
        """
        Save microstate probabilities to a text file, ensuring they sum to 1.
        """
        microstates_probs /= np.sum(microstates_probs)  # Normalize probabilities to sum to 1
        logging.info(f"Sum of microstate probabilities: {np.sum(microstates_probs)}")
        output_file = self.save_dir / f"{self.struct}-{self.struct_type}_microst_p.txt"
        np.savetxt(output_file, microstates_probs, fmt="%.6f")
        logging.info(f"Microstate probabilities saved to {output_file}")

    def processing(self):
        """
        Execute the KL probabilities calculation pipeline.
        """
        logging.info("Starting clustering.")
        cl_labels = self.clustering(self.dist_matrix)
        microstates_probs = self.labeling(cl_labels)
        self.savefile(microstates_probs)

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform KL clustering on molecular dynamics trajectories.")
    
    parser.add_argument("-m", "--dist_mat", required=True, type=str, help="Path to the distance matrix (.npy).")
    parser.add_argument("-g", "--gro", required=True, type=str, help="Path to the GRO file.")
    parser.add_argument("-x", "--xtc", required=True, type=str, help="Path to the XTC trajectory file.")
    parser.add_argument("-s", "--save_dir", required=True, type=str, help="Directory to save outputs.")
    parser.add_argument("-c", "--struct", required=True, type=str, help="Molecule name.")
    parser.add_argument("-t", "--struct_type", required=True, type=str, help="Type of the structure (e.g. apo or name of ligand for holo).")
    parser.add_argument("-l", "--clusters", required=True, type=int, help="Number of clusters (L).")
    
    args = parser.parse_args()

    dist_mat = Path(args.dist_mat)
    gro_file = Path(args.gro)
    traj_file = Path(args.xtc)
    save_dir = Path(args.save_dir)
    struct_type = args.struct_type
    struct = args.struct
    L = args.clusters

    if not dist_mat.exists():
        raise FileNotFoundError(f"The distance matrix file '{dist_mat}' does not exist.")
    if not gro_file.exists():
        raise FileNotFoundError(f"The GRO file '{gro_file}' does not exist.")
    if not traj_file.exists():
        raise FileNotFoundError(f"The XTC file '{traj_file}' does not exist.")

    save_dir.mkdir(parents=True, exist_ok=True)

    klp = KLProbabilities(
        dist_mat=str(dist_mat),
        gro_file=str(gro_file),
        traj_file=str(traj_file),
        save_dir=str(save_dir),
        struct_type=struct_type,
        struct=struct,
        L=L
    )
    
    klp.processing()

