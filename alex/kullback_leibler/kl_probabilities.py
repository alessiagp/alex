import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from pathlib import Path
from collections import Counter
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align
import os
import sys
import random

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class KLProbabilities:
    def __init__(self, gro_ref: str, traj_ref: str, selection: str, chosen_struct: str, L: int, save_dir: str, save_prefix: str):
        self.gro_file = Path(gro_ref)
        self.traj_file = Path(traj_ref)
        self.selection = selection.strip()
        self.chosen_struct = chosen_struct
        self.L = L
        
        self.save_dir = Path(save_dir)
        self.save_prefix = save_prefix
        self.save_dir.mkdir(parents=True, exist_ok=True)

        # Initialize MDAnalysis Universe
        self.t = mda.Universe(str(self.gro_file), str(self.traj_file))
        self.dist_matrix = None

        # Validate selection
        try:
            selected_atoms = self.t.select_atoms(self.selection)
            if len(selected_atoms) == 0:
                raise ValueError(f"Selection '{self.selection}' did not match any atoms.")
        except Exception as e:
            raise ValueError(f"Selection error: {e}")

    def rmsd_matrix(self) -> np.ndarray:
        """
        Align the trajectory with respect to the average structure and compute the diffusion matrix.
        """
        logging.info(f"Starting alignment for {self.chosen_struct}")
        avg = align.AverageStructure(self.t, select=self.selection).run()
        align.AlignTraj(self.t, avg.results.universe, select=self.selection, in_memory=True).run()
        logging.info("Alignment completed.")

        logging.info("Calculating diffusion matrix.")
        self.dist_matrix = diffusionmap.DistanceMatrix(self.t, select=self.selection).run().results.dist_matrix
        
        return self.dist_matrix

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
        label_probs = np.array([label_counts[el] / total for el in cl_labels])

        logging.info(f"Length of label_probs: {len(label_probs)}")
        logging.info("Label probabilities computed.")

        # Select one random frame per cluster
        unique_classes = list(label_counts.keys())
        cl_selection = np.array([random.choice(np.where(cl_labels == class_label)[0]) for class_label in unique_classes])

        logging.info(f"Selected frames: {cl_selection}")
        logging.info(f"Number of selected frames: {len(cl_selection)}")

        # Output file path
        output_traj = self.save_dir / f"{self.save_prefix}-{self.chosen_struct}-KL_frames.xtc"

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
        Save microstate probabilities to a text file.
        """
        output_file = self.save_dir / f"{self.save_prefix}-{self.chosen_struct}_microst_p.txt"
        np.savetxt(output_file, microstates_probs, fmt="%.6f")
        logging.info(f"Microstate probabilities saved to {output_file}")

    def processing(self):
        """
        Execute the KL probabilities calculation pipeline.
        """
        logging.info("Starting RMSD matrix calculation.")
        dist_matrix = self.rmsd_matrix()
        
        logging.info("Starting clustering.")
        cl_labels = self.clustering(dist_matrix)
        microstates_probs = self.labeling(cl_labels)
        self.savefile(microstates_probs)

if __name__ == '__main__':
    if len(sys.argv) < 5:
        raise ValueError("Usage: kl_probabilities_serial.py <gro_filename> <traj_filename> <chosen_struct> <L> <selection> <save_prefix>")

    chosen_struct = sys.argv[3].strip()
    try:
        L = int(sys.argv[4])
        if L <= 0:
            raise ValueError("The 'L' argument must be a positive integer.")
    except ValueError as e:
        raise ValueError("Invalid 'L' argument. It should be a positive integer.") from e
    
    selection = sys.argv[5].strip()
    save_prefix = sys.argv[6].strip()

    main_dir = Path.cwd().parent
    save_dir = main_dir / chosen_struct / "KL_probs"
    datadir = main_dir / chosen_struct / "data"
    
    gro_file = datadir / f"{sys.argv[1].strip()}"
    traj_file = datadir / f"{sys.argv[2].strip()}"
    
    if not gro_file.exists():
        raise FileNotFoundError(f"The GRO file '{gro_file}' does not exist.")
    if not traj_file.exists():
        raise FileNotFoundError(f"The XTC file '{traj_file}' does not exist.")
    
    save_dir.mkdir(exist_ok=True)
    
    klp = KLProbabilities(
        gro_ref=str(gro_file),
        traj_ref=str(traj_file),
        selection=selection,
        chosen_struct=chosen_struct,
        L=L,
        save_dir=str(save_dir),
        save_prefix=save_prefix
    )
    klp.processing()
