import os
import re
import sys
import numpy as np


class RandomSMapProcessor:
    def __init__(self, filepath, nsmap=48):
        self.filepath = filepath
        self.nsmap_expected = nsmap
        self.smaps = []

        # Extract protein name from:
        # protein_name_random_N102.dat
        filename = os.path.basename(filepath)

        match = re.match(r"(.+)_random_N\d+\.dat$", filename)

        if not match:
            print(
                f"Error: Input file '{filename}' does not match the expected "
                f"format 'protein_name_random_N102.dat'."
            )
            sys.exit(1)

        self.protein_name = match.group(1)

        # Output in the same directory as the input file
        output_dir = os.path.dirname(os.path.abspath(filepath))

        self.smap_filepath = os.path.join(
            output_dir,
            f"{self.protein_name}-smaps.txt"
        )

        self._validate_inputs()

    def _validate_inputs(self):
        """Check that the input file exists."""
        if not os.path.isfile(self.filepath):
            print(f"Error: File '{self.filepath}' does not exist.")
            sys.exit(1)

    def _extract_smaps(self):
        """Extract every random_smap value from the file."""

        extracted_smaps = []

        with open(self.filepath, "r") as f:
            for line in f:

                match = re.search(
                    r"random_smap\s*=\s*"
                    r"([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
                    line
                )

                if match:
                    extracted_smaps.append(
                        float(match.group(1))
                    )

        return extracted_smaps

    def process_file(self):
        """Extract random_smap values and calculate statistics."""

        print(f"\nReading file: {self.filepath}")

        self.smaps = self._extract_smaps()

        nsmap = len(self.smaps)

        print(f"Found {nsmap} random_smap values.")

        if nsmap != self.nsmap_expected:
            print(
                f"Warning: Expected {self.nsmap_expected} values, "
                f"but found {nsmap}."
            )

        if nsmap == 0:
            print("No random_smap values found. Exiting.")
            return

        self._write_results()
        self._print_statistics()

    def _write_results(self):
        """Write individual random_smap values."""

        with open(self.smap_filepath, "w") as f:
            for i, smap in enumerate(self.smaps, start=1):
                f.write(f"{i}\t{smap:.6f}\n")

        print(f"\nValues written to:")
        print(f"  {self.smap_filepath}")

    def _print_statistics(self):
        """Calculate and print average and sample standard deviation."""

        smap_array = np.array(self.smaps)

        mean_smap = np.mean(smap_array)
        std_smap = np.std(smap_array, ddof=1)

        print("\nStatistics:")
        print(f"  Average random_smap:      {mean_smap:.6f}")
        print(f"  Standard deviation:       {std_smap:.6f}")


# ==============================
# Main execution
# ==============================

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print(
            "Usage: python3 random_smap_analysis.py "
            "<random_output_file>"
        )
        sys.exit(1)

    filepath = sys.argv[1]

    processor = RandomSMapProcessor(
        filepath,
        nsmap=48
    )

    processor.process_file()