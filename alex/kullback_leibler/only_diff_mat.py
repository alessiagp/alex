import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align
import argparse


# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class KLProbabilities:
    def __init__(self, gro_file: str, traj_file: str, selection: str, struct: str, save_dir: str, struct_type: str):
        self.gro_file = Path(gro_file)
        self.traj_file = Path(traj_file)
        self.selection = selection.strip()
        self.struct = struct
        self.struct_type = struct_type
        
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

        # Initialize MDAnalysis Universe
        self.t = mda.Universe(str(self.gro_file), str(self.traj_file))
        self.dist_matrix: np.ndarray | None = None

        # Validate selection
        try:
            selected_atoms = self.t.select_atoms(self.selection)
            if len(selected_atoms) == 0:
                raise ValueError(f"Selection '{self.selection}' did not match any atoms.")
        except Exception as e:
            raise ValueError(f"Selection error: {e}")

    def rmsd_matrix(self) -> np.ndarray:
        """
        Align the trajectory with respect to the average structure and compute the pairwise distance matrix.
        """
        logging.info(f"Starting alignment for {self.struct}")
        avg = align.AverageStructure(self.t, select=self.selection).run()
        align.AlignTraj(self.t, avg.results.universe, select=self.selection, in_memory=True).run()
        logging.info("Alignment completed.")

        logging.info("Calculating pairwise distance matrix.")
        self.dist_matrix = diffusionmap.DistanceMatrix(self.t, select=self.selection).run().results.dist_matrix

        out_file = self.save_dir / f"{self.struct}-{self.struct_type}-RMSD_matrix.npy"
        np.save(out_file, self.dist_matrix)
        logging.info(f"Distance matrix saved to {out_file}")

        return self.dist_matrix

    def sns_plot(self, mat1: np.ndarray):
        """
        Generate and save an RMSD matrix heatmap.
        """
        if mat1.shape[0] > 500:
            logging.warning("Matrix is large; consider downsampling for better visualization.")
        
        logging.info("Generating heatmap.")
        plt.figure(figsize=(10, 8))
        sns.heatmap(mat1, cmap="inferno", cbar_kws={"label": "RMSD (Å)"})
        plt.xlabel("Frames")
        plt.ylabel("Frames")
        plt.title(f"RMSD diffusion matrix for {self.struct_type} {self.struct}")

        out_png = self.save_dir / f"{self.struct}-{self.struct_type}-RMSD_matrix.png"
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()
        logging.info(f"Heatmap saved to {out_png}")

    def processing(self):
        """
        Execute the KL probabilities calculation pipeline.
        """
        logging.info("Starting RMSD matrix calculation.")
        dist_matrix = self.rmsd_matrix()
        self.sns_plot(dist_matrix)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate the RMSD matrix of a molecular dynamics trajectory.")
    
    parser.add_argument("-g", "--gro", required=True, type=str, help="Path to the GRO file.")
    parser.add_argument("-x", "--xtc", required=True, type=str, help="Path to the XTC trajectory file.")
    parser.add_argument("-s", "--save_dir", required=True, type=str, help="Directory to save outputs.")
    parser.add_argument("-c", "--struct", required=True, type=str, help="Molecule name.")
    parser.add_argument("-st", "--struct_type", required=True, type=str, help="Type of the structure (e.g. apo or name of ligand for holo).")
    parser.add_argument("-sel", "--selection", required=True, type=str, help="Selection string for atoms.")
    
    args = parser.parse_args()

    gro_file = Path(args.gro)
    traj_file = Path(args.xtc)
    save_dir = Path(args.save_dir)

    if not gro_file.exists():
        raise FileNotFoundError(f"The GRO file '{gro_file}' does not exist.")
    if not traj_file.exists():
        raise FileNotFoundError(f"The XTC file '{traj_file}' does not exist.")

    save_dir.mkdir(parents=True, exist_ok=True)

    klp = KLProbabilities(
        gro_file=str(gro_file),
        traj_file=str(traj_file),
        save_dir=str(save_dir),
        struct_type=args.struct_type,
        struct=args.struct,
        selection=args.selection,
    )
    
    klp.processing()
