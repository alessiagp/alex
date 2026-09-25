import os
import re
import sys
import numpy as np


class SMapProcessor:
    def __init__(self, opt_name, nsmap=48):
        self.workdir = os.getcwd()
        self.directory = os.path.join(self.workdir, "optimize-results")

        self.opt_name = opt_name
        self.nsmap_expected = nsmap

        self.smaps = []

        self.smap_filepath = os.path.join(
            self.directory,
            f"last_smaps_{self.opt_name}.txt"
        )

        self._validate_inputs()

    def _validate_inputs(self):
        """Validate input arguments and ensure directory exists."""
        if not os.path.exists(self.directory):
            print(f"Error: Directory '{self.directory}' does not exist.")
            sys.exit(1)

    def _extract_last_smap(self, filepath):
        """Extract the last_smap value from an optimization file."""

        with open(filepath, "r") as f:
            lines = f.readlines()

        # Search backwards because last_smap is near the end
        for line in reversed(lines):

            match = re.search(
                r"last_smap\s+"
                r"([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
                line
            )

            if match:
                return float(match.group(1))

        print(
            f"Warning: No last_smap found in "
            f"{os.path.basename(filepath)}"
        )

        return None

    def process_files(self):
        """Iterate through optimization files and extract last_smap values."""

        print("\nSearching for optimization files...")

        for file in sorted(os.listdir(self.directory)):

            filepath = os.path.join(self.directory, file)

            if not os.path.isfile(filepath):
                continue

            if not file.startswith(self.opt_name):
                continue

            smap = self._extract_last_smap(filepath)

            if smap is not None:
                self.smaps.append(smap)

                print(
                    f"{file}: last_smap = {smap:.6f}"
                )

        nsmap = len(self.smaps)

        print(f"\nFound {nsmap} last_smap values.")

        if nsmap != self.nsmap_expected:
            print(
                f"Warning: Expected {self.nsmap_expected} values, "
                f"but found {nsmap}."
            )

        if nsmap == 0:
            print("No last_smap values found.")
            return

        self._write_results()
        self._print_statistics()

    def _write_results(self):
        """Write individual last_smap values."""

        with open(self.smap_filepath, "w") as f:

            for i, smap in enumerate(self.smaps, start=1):
                f.write(f"{smap:.6f}\n")

        print("\nValues written to:")
        print(f"  {self.smap_filepath}")

    def _print_statistics(self):
        """Calculate and print average and sample standard deviation."""

        smap_array = np.array(self.smaps)

        mean_smap = np.mean(smap_array)

        # Sample standard deviation
        std_smap = np.std(smap_array, ddof=1)

        print("\nStatistics:")
        print(f"  Average last_smap:      {mean_smap:.6f}")
        print(f"  Standard deviation:    {std_smap:.6f}")


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

    processor = SMapProcessor(
        opt_name,
        nsmap=48
    )

    processor.process_files()