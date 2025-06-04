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

class LoadFiles:
    def __init__(self, save_dir: str, membrane = None, struct_type = None, analysis_mode = 'rmsd_mat'):
        self.membrane = membrane
        self.struct_type = struct_type #ligand

        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

        self.analysis_mode = analysis_mode



class RMSDMatrixEquilibration(LoadFiles):
    def __init__(self, gro_file: str, traj_file: str, selection: str): 
        self.analysis_mode = 'rmsd_mat'
        self.gro_file = Path(gro_file)
        self.traj_file = Path(traj_file)
        self.selection = selection.strip()

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

    def calc_rmsd_matrix(self) -> np.ndarray:
        """
        Align the trajectory with respect to the average structure and compute the diffusion matrix.
        """
        logging.info(f"Starting alignment for {self.chosen_struct}")
        avg = align.AverageStructure(self.t, select=self.selection).run()
        align.AlignTraj(self.t, avg.results.universe, select=self.selection, in_memory=True).run()
        logging.info("Alignment completed.")

        logging.info("Calculating diffusion matrix.")
        self.dist_matrix = diffusionmap.DistanceMatrix(self.t, select=self.selection).run().results.dist_matrix
        print(self.dist_matrix)
        np.save(f"{self.chosen_struct}-{self.struct_type}-rmsd_diffmat.npy", self.dist_matrix)
        
        return self.dist_matrix
    
    def sns_plot(self, mat1: np.ndarray):
        """
        Generate and save an RMSD matrix heatmap.
        """
        if mat1.shape[0] > 500:
            logging.warning("Matrix is large; consider downsampling for better visualization.")
        
        logging.info("Generating heatmap.")
        sns.heatmap(mat1, cmap='inferno', cbar_kws={'label': 'RMSD (Angstrom)'})
        plt.xlabel("Frames")
        plt.ylabel("Frames")
        plt.title(f"RMSD diffusion matrix for {self.struct_type} {self.chosen_struct}")
        plt.savefig(self.save_dir / f"{self.struct_type}-{self.chosen_struct}-RMSD_matrix.png")
        plt.close()
        logging.info("Heatmap saved.")

    def use_rmsd_matrix(self):
        """
        Check whether the trajectory is equilibrated by looking at the RMSD matrix.
        """
        logging.info("Starting RMSD matrix calculation.")
        dist_matrix = self.calc_rmsd_matrix()
        
        logging.info("Calculation completed.")
        self.sns_plot(dist_matrix)



class BlockAnalysisEquilibration(LoadFiles):
    def __init__(self, data: np.ndarray, chunk_length: int):
        super().__init__(save_dir=".", analysis_mode="block")  # Optional
        self.data = np.asarray(data)
        self.chunk_length = chunk_length 
        self.block_bse = None
        self.block_sizes = None


    def block_analysis(self) -> float:
        """
        Compute the block standard error (BSE) for a specific chunk length.
        """
        n = len(self.data)
        if self.chunk_length < 2 or self.chunk_length > n // 2:
            return np.nan

        considered_length = (n // self.chunk_length) * self.chunk_length
        n_blocks = considered_length // self.chunk_length

        if n_blocks < 2:
            return np.nan

        blocks = self.data[:considered_length].reshape((n_blocks, self.chunk_length))
        block_means = blocks.mean(axis=1)
        bse = block_means.std(ddof=1) / np.sqrt(n_blocks)

        return bse

    def assess_convergence(self, min_chunk=2, max_chunk=None):
        """
        Perform block analysis across a range of chunk lengths.
        """
        n = len(self.data)
        max_chunk = max_chunk or n // 2
        self.block_sizes = np.arange(min_chunk, max_chunk)
        self.block_bse = np.array([self.block_analysis(cl) for cl in self.block_sizes])

        return self.block_bse, self.block_sizes

    def estimate_autocorrelation(self, total_steps=None):
        """
        Estimate autocorrelation time and number of independent samples.
        """
        if self.block_bse is None or self.block_sizes is None:
            raise ValueError("Run assess_convergence() before this.")

        valid_bse = self.block_bse[~np.isnan(self.block_bse)]
        if len(valid_bse) < 2:
            return np.nan, np.nan

        BSE_final = valid_bse[-1]
        variance = np.var(valid_bse)

        if variance == 0:
            return np.nan, np.nan

        total_steps = total_steps or len(self.data)
        tau = (BSE_final ** 2) / variance * total_steps
        N_eff = variance / (BSE_final ** 2)

        return tau, N_eff

    def plot_bse(self):
        """
        Plot the block standard error (BSE) vs block size.
        """
        if self.block_bse is None or self.block_sizes is None:
            raise ValueError("Run assess_convergence() before plotting.")

        plt.figure(figsize=(8, 5))
        plt.plot(self.block_sizes, self.block_bse, 'o-', label='Block Std Error')
        plt.xlabel("Block size")
        plt.ylabel("BSE")
        plt.title("Block Analysis of Time Series")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    def use_block_analysis(self):
        max_chunk = max(2, len(self.data) // 20)  # Ensure at least 20 blocks
        block_bse, block_sizes = self.assess_convergence(max_chunk=max_chunk)
        tau, N_eff = self.estimate_autocorrelation(total_steps=len(self.data))

        logging.info(f"Estimated autocorrelation time: {tau:.2f}")
        logging.info(f"Estimated number of independent samples: {N_eff:.2f}")

        self.plot_bse()

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Assess equilibration of MD trajectories.")
    
    parser.add_argument("-s", "--save_dir", required=True, type=str, help="Directory to save outputs.")
    parser.add_argument("-m", "--membrane", required=False, type=str, help="Membrane label.")
    parser.add_argument("-t", "--struct_type", required=True, type=str, help="Type of molecule name.")
    parser.add_argument("-a", "--analysis_mode", required=True, choices=['rmsd_mat', 'block'],
                        help="Choose 'rmsd_mat' for RMSD matrix, 'block' for block analysis.")

    # RMSD mode args
    parser.add_argument("-g", "--gro", type=str, help="Path to the GRO file.")
    parser.add_argument("-x", "--xtc", type=str, help="Path to the XTC trajectory file.")
    parser.add_argument("-sel", "--selection", type=str, help="Atom selection string.")

    # Block mode args
    parser.add_argument("-bd", "--block_data", type=str, help="Path to .npy file for block analysis.")
    parser.add_argument("-bl", "--chunk_length", type=int, help="Chunk length for block analysis.")

    args = parser.parse_args()

    if args.analysis_mode == "rmsd_mat":
        if not all([args.gro, args.xtc, args.selection]):
            raise ValueError("GRO, XTC, and selection must be provided for RMSD matrix mode.")

        rmsd_mat = RMSDMatrixEquilibration(
            gro_file=args.gro,
            traj_file=args.xtc,
            selection=args.selection
        )
        rmsd_mat.chosen_struct = args.membrane  # Needed for plotting title/path
        rmsd_mat.struct_type = args.struct_type
        rmsd_mat.save_dir = Path(args.save_dir)
        rmsd_mat.use_rmsd_matrix()

    elif args.analysis_mode == "block":
        if not all([args.block_data, args.chunk_length]):
            raise ValueError("Block data file and chunk length are required for block mode.")
        block_path = Path(args.block_data)
        if not block_path.exists():
            raise FileNotFoundError(f"Block data file not found: {block_path}")

        data = np.load(block_path)
        block_analysis = BlockAnalysisEquilibration(data=data, chunk_length=args.chunk_length)
        block_analysis.use_block_analysis()
