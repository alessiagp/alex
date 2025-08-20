import mdtraj as md
import argparse
import os

class GromacsConverter:
    def __init__(self, xtc_path, gro_path, xyz_filename, hydrogens=False):
        self.xtc_path = xtc_path
        self.gro_path = gro_path
        self.xyz_filename = xyz_filename
        self.hydrogens = hydrogens
        self._validate_files()

    def _validate_files(self):
        if not os.path.exists(self.xtc_path):
            raise FileNotFoundError(f"XTC file not found: {self.xtc_path}")
        if not os.path.exists(self.gro_path):
            raise FileNotFoundError(f"GRO file not found: {self.gro_path}")

    def convert(self):
        print("Loading trajectory...")
        full_traj = md.load_xtc(self.xtc_path, top=self.gro_path)

        if not self.hydrogens:
            print("Selecting non-hydrogen atoms...")
            no_h = full_traj.topology.select('element != H')
            heavy_traj = full_traj.atom_slice(no_h)

            print(f"Full trajectory has {full_traj.n_atoms} atoms.")
            print(f"Filtered trajectory has {len(no_h)} heavy atoms.")
            print(f"Number of frames to convert: {heavy_traj.n_frames}")

            print(f"Saving XYZ file: {self.xyz_filename}")
            with md.formats.XYZTrajectoryFile(self.xyz_filename, 'w') as f:
                f.write(heavy_traj.xyz * 10)  # Convert to Angstroms
            print("Conversion completed successfully.")
        else:
            print(f"Full trajectory has {full_traj.n_atoms} atoms.")
            print(f"Number of frames to convert: {full_traj.n_frames}")

            print(f"Saving XYZ file: {self.xyz_filename}")
            with md.formats.XYZTrajectoryFile(self.xyz_filename, 'w') as f:
                f.write(full_traj.xyz * 10)  # Convert to Angstroms
            print("Conversion completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GROMACS XTC to XYZ, optionally removing hydrogen atoms.")
    parser.add_argument("--xtc", required=True, help="Path to the XTC file")
    parser.add_argument("--gro", required=True, help="Path to the GRO file")
    parser.add_argument("--xyz", required=True, help="Path to output XYZ file")
    parser.add_argument("--hydrogens", action="store_true", help="Include hydrogen atoms in the output")

    args = parser.parse_args()

    converter = GromacsConverter(args.xtc, args.gro, args.xyz, args.hydrogens)
    converter.convert()
