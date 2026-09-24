import os
import re
import sys
import numpy as np


class RandomMappingProcessor:
    def __init__(self, filepath, natoms, nmaps=48):
        self.filepath = filepath
        self.natoms = natoms
        self.nmaps_expected = nmaps

        self.mappings = []

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

        # Write outputs in the same directory as the input file
        output_dir = os.path.dirname(os.path.abspath(filepath))

        self.map_filepath = os.path.join(
            output_dir,
            f"{self.protein_name}-MAPPINGS.txt"
        )

        self.probs_filepath = os.path.join(
            output_dir,
            f"{self.protein_name}-probabilities.txt"
        )

        self.error_filepath = os.path.join(
            output_dir,
            f"{self.protein_name}-error.txt"
        )

        self._validate_inputs()

    def _validate_inputs(self):
        """Check that the input file exists."""
        if not os.path.isfile(self.filepath):
            print(f"Error: File '{self.filepath}' does not exist.")
            sys.exit(1)

    def _extract_mappings(self):
        """
        Extract every mapping following 'conv mapping'.

        Example:

            conv mapping
            6 17 23 39 41 ...

        Only the row immediately following 'conv mapping' is extracted.
        """

        extracted_mappings = []

        with open(self.filepath, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):

            if line.strip() != "conv mapping":
                continue

            if i + 1 >= len(lines):
                print("Warning: 'conv mapping' found at end of file.")
                continue

            mapping_line = lines[i + 1].strip()

            if re.fullmatch(r"\d+(?:\s+\d+)*", mapping_line):
                mapping = list(map(int, mapping_line.split()))
                extracted_mappings.append(mapping)

            else:
                print(
                    f"Warning: Invalid mapping after 'conv mapping' "
                    f"at line {i + 1}."
                )

        return extracted_mappings

    def make_counts(self, mapping_matrix, nmaps):
        """Calculate the frequency with which each atom appears."""

        if nmaps == 0:
            print("Warning: No mappings found.")
            return [0] * self.natoms

        counts = [0] * self.natoms

        for row in mapping_matrix:
            for value in row:

                if 0 <= value < self.natoms:
                    counts[value] += 1
                else:
                    print(
                        f"Warning: Atom index {value} is outside "
                        f"the expected range 0-{self.natoms - 1}."
                    )

        return [count / nmaps for count in counts]

    def process_file(self):
        """Extract mappings and calculate probabilities."""

        print(f"\nReading file: {self.filepath}")

        self.mappings = self._extract_mappings()

        nmaps = len(self.mappings)

        print(f"Found {nmaps} mappings.")

        if nmaps != self.nmaps_expected:
            print(
                f"Warning: Expected {self.nmaps_expected} mappings, "
                f"but found {nmaps}."
            )

        if nmaps == 0:
            print("No mappings found. Exiting.")
            return

        self._write_results()

    def _write_results(self):
        """Write mappings, probabilities, and 95% confidence intervals."""

        nmaps = len(self.mappings)

        # ----------------------------------
        # Write mappings
        # ----------------------------------

        print("\nWriting mappings...")

        with open(self.map_filepath, "w") as f:
            for mapping in self.mappings:
                f.write(" ".join(map(str, mapping)) + "\n")

        # ----------------------------------
        # Calculate probabilities
        # ----------------------------------

        print("Calculating probabilities...")

        probabilities = self.make_counts(
            self.mappings,
            nmaps
        )

        with open(self.probs_filepath, "w") as f:
            f.write("\n".join(map(str, probabilities)))

        # ----------------------------------
        # Calculate 95% confidence interval
        # ----------------------------------

        print("Calculating standard error and 95% confidence interval...")

        prob_array = np.array(probabilities)

        variance = prob_array * (1.0 - prob_array)
        standard_error = np.sqrt(variance / nmaps)
        errors_95ci = 1.96 * standard_error

        with open(self.error_filepath, "w") as f:
            f.write("\n".join(map(str, errors_95ci)))

        print("\nResults written successfully:")
        print(f"  {self.map_filepath}")
        print(f"  {self.probs_filepath}")
        print(f"  {self.error_filepath}")


# ==============================
# Main execution
# ==============================

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print(
            "Usage: python3 random_maps.py "
            "<random_output_file> <num_atoms>"
        )
        sys.exit(1)

    filepath = sys.argv[1]
    natoms = int(sys.argv[2])

    processor = RandomMappingProcessor(
        filepath,
        natoms,
        nmaps=48
    )

    processor.process_file()