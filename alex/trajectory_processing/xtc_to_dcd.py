import mdtraj as md
import argparse
import os

class GromacsConverter:
    def __init__(self, xtc_path, gro_path, dcd_filename):
        self.xtc_path = xtc_path
        self.gro_path = gro_path
        self.dcd_filename = dcd_filename
        self._validate_files()

    def _validate_files(self):
        if not os.path.exists(self.xtc_path):
            raise FileNotFoundError(f"XTC file not found: {self.xtc_path}")
        if not os.path.exists(self.gro_path):
            raise FileNotFoundError(f"GRO file not found: {self.gro_path}")

    def convert(self):
        print("Loading trajectory...")
        t = md.load_xtc(self.xtc_path, top=self.gro_path)

        print(f"Trajectory has {t.n_atoms} atoms.")
        print(f"Number of frames to convert: {t.n_frames}")

        print(f"Saving DCD file: {self.dcd_filename}")
        t.save_dcd(f"{self.dcd_filename}.dcd")
        print("Conversion completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GROMACS XTC to DCD")
    parser.add_argument("--xtc", required=True, help="Path to the XTC file")
    parser.add_argument("--gro", required=True, help="Path to the GRO file")
    parser.add_argument("--dcd", required=True, help="Path to output DCD file (no extension)")

    args = parser.parse_args()

    converter = GromacsConverter(args.xtc, args.gro, args.dcd)
    converter.convert()
