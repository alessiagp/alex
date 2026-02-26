import numpy as np
import logging
from pathlib import Path
import argparse

import MDAnalysis as mda
from MDAnalysis.analysis import diffusionmap, align
from scipy.spatial.distance import squareform

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


class KLProbabilities:
    def __init__(
        self,
        gro_file: str,
        traj_file: str,
        selection: str,
        struct: str,
        save_dir: str,
        struct_type: str,
        stride: int = 1,
        align_to_avg: bool = False
    ):
        if stride < 1:
            raise ValueError("stride must be >= 1")

        self.gro_file = Path(gro_file)
        self.traj_file = Path(traj_file)
        self.selection = selection.strip()
        self.struct = struct
        self.struct_type = struct_type
        self.stride = stride
        self.align_to_avg = align_to_avg

        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

        # Initialize MDAnalysis Universe
        self.t = mda.Universe(str(self.gro_file), str(self.traj_file))
        self.ref = mda.Universe(str(self.gro_file))

        # Validate selection
        try:
            selected_atoms = self.t.select_atoms(self.selection)
            if len(selected_atoms) == 0:
                raise ValueError(f"Selection '{self.selection}' did not match any atoms.")
        except Exception as e:
            raise ValueError(f"Selection error: {e}")

    def rmsd_matrix(self) -> np.ndarray:
        """
        Align trajectory and compute condensed pairwise RMSD distance matrix.
        """
        logging.info(f"Starting alignment for {self.struct} (stride={self.stride}, align_to_avg={self.align_to_avg})")
        logging.info(
        "Alignment mode: %s", "average structure" if align_to_avg else "GRO reference")

        if self.align_to_avg:
            logging.info("Aligning trajectory to average structure.")
            avg = align.AverageStructure(
                self.t,
                select=self.selection,
                step=self.stride
            ).run()

            align.AlignTraj(
                self.t,
                avg.results.universe,
                select=self.selection,
                in_memory=True,
                step=self.stride
            ).run()
        else:
            logging.info("Aligning trajectory to reference GRO structure.")
            align.AlignTraj(
                self.t,
                self.ref,
                select=self.selection,
                in_memory=True,
                step=self.stride
            ).run()

        logging.info("Alignment completed.")
        logging.info("Calculating pairwise distance matrix.")

        full_matrix = diffusionmap.DistanceMatrix(
            self.t,
            select=self.selection,
            step=self.stride
        ).run().results.dist_matrix

        condensed = squareform(
            full_matrix,
            force="tovector",
            checks=False
        )

        out_file = (
            self.save_dir
            / f"{self.struct_type}-{self.struct}-RMSD_matrix_stride{self.stride}_condensed.npy"
        )
        np.save(out_file, condensed)

        logging.info(
            f"Condensed distance matrix saved to {out_file} "
            f"(length={condensed.shape[0]})"
        )

        return condensed

    def processing(self):
        logging.info("Starting RMSD matrix calculation.")
        self.rmsd_matrix()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the RMSD matrix of a molecular dynamics trajectory."
    )

    parser.add_argument("-g", "--gro", required=True, help="Path to the GRO file.")
    parser.add_argument("-x", "--xtc", required=True, help="Path to the XTC trajectory file.")
    parser.add_argument("-s", "--save_dir", required=True, help="Directory to save outputs.")
    parser.add_argument("-c", "--struct", required=True, help="Molecule name.")
    parser.add_argument("-st", "--struct_type", required=True, help="Structure type (e.g. holo/apo).")
    parser.add_argument("-sel", "--selection", required=True, help="Atom selection string.")
    parser.add_argument("--stride", type=int, default=1, help="Stride for trajectory frames.")
    
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--align-to-avg", action="store_true", help="Align trajectory to the average structure.")

    group.add_argument("--align-to-ref", action="store_true", help="Align trajectory to the reference GRO structure (default).")

    args = parser.parse_args()

    align_to_avg = False

    if args.align_to_avg:
        align_to_avg = True
    elif args.align_to_ref:
        align_to_avg = False

    gro_file = Path(args.gro)
    traj_file = Path(args.xtc)
    save_dir = Path(args.save_dir)

    if not gro_file.exists():
        raise FileNotFoundError(f"The GRO file '{gro_file}' does not exist.")
    if not traj_file.exists():
        raise FileNotFoundError(f"The XTC file '{traj_file}' does not exist.")

    klp = KLProbabilities(
        gro_file=str(gro_file),
        traj_file=str(traj_file),
        save_dir=str(save_dir),
        struct_type=args.struct_type,
        struct=args.struct,
        selection=args.selection,
        stride=args.stride,
        align_to_avg=align_to_avg
    )

    klp.processing()