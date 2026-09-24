import os
import re
import sys
import numpy as np


class SMapProcessor:
    def __init__(self, opt_name, nsmap=48):
        self.workdir = os.getcwd()
        self.directory = os.path.join(self.workdir, "optimize-results")

        self.opt_name = opt_name
        self.nsmap = nsmap

        self.smaps = []
        self.smap_filepath = ""
        self.stats_filepath = ""

        self._validate_inputs()
        self._set_output_files()

    def _validate_inputs(self):
        """Validate input arguments and ensure directory exists."""
        if not os.path.exists(self.directory):
            print(f"Error: Directory '{self.directory}' does not exist.")
            sys.exit(1)

    def _set_output_files(self):
        """Set output file paths."""
        self.smap_filepath = os.path.join(
            self.directory,
            f"last_smaps_{self.opt_name}.txt"
        )

        self.stats_filepath = os.path.join(
            self.directory,
            f"last_smaps_stats_{self.opt_name}.txt"
        )

    def _extract_last_smap(self, filepath):
        """Extract the last_smap value from an optimization file."""

        with open(filepath, "r") as f:
            lines = f.readlines()

        # Search backwards, since last_smap is expected near the end
        for line in reversed(lines):
            match = re.search(
                r"last_smap\s+([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
                line
            )

            if match:
                return float(match.group(1))

        print(f"Warning: No last_smap found in {os.path.basename(filepath)}")
        return None

    def process_files(self):
        """Iterate through optimization files and extract last_smap values."""

        print("\nSearching for optimization files...")

        for file in sorted(os.listdir(self.directory)):

            if not file.startswith(self.opt_name):
                continue

            filepath = os.path.join(self.directory, file)

            # Skip directories
            if not os.path.isfile(filepath):
                continue

            smap = self._extract_last_smap(filepath)

            if smap is not None:
                self.smaps.append(smap)
                print(f"{file}: last_smap = {smap}")

        print(f"\nFound {len(self.smaps)} last_smap values.")

        if len(self.smaps) != self.nsmap:
            print(
                f"Warning: Expected {self.nsmap} values, "
                f"but found {len(self.smaps)}."
            )

        self._write_results()

    def _write_results(self):
        """Write individual smaps, average, and standard deviation."""

        if not self.smaps:
            print("Warning: No smap values found.")
            return

        # Convert to numpy array
        smap_array = np.array(self.smaps)

        # Average
        mean_smap = np.mean(smap_array)

        # Standard deviation of the sample
        std_smap = np.std(smap_array, ddof=1)

        # ----------------------------------
        # Write individual last_smap values
        # ----------------------------------

        with open(self.smap_filepath, "w") as f:
            for smap in self.smaps:
                f.write(f"{smap:.6f}\n")

        # ----------------------------------
        # Write statistics
        # ----------------------------------

        with open(self.stats_filepath, "w") as f:
            f.write(f"Average last_smap: {mean_smap:.6f}\n")
            f.write(f"Standard deviation: {std_smap:.6f}\n")

        print(f"\nResults written successfully:")
        print(f"  {self.smap_filepath}")
        print(f"  {self.stats_filepath}")

        print(f"\nAverage last_smap = {mean_smap:.6f}")
        print(f"Standard deviation = {std_smap:.6f}")


# ==============================
# Main execution
# ==============================

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print(
            "Error: Missing argument.\n"
            "Usage: python3 smap_analysis.py <optimization_name>"
        )
        sys.exit(1)

    opt_name = sys.argv[1]

    processor = SMapProcessor(opt_name)
    processor.process_files()