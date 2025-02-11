import os
import re
import sys

class MappingProcessor:
    def __init__(self, opt_name, natoms):
        self.workdir = os.getcwd()
        self.directory = os.path.join(self.workdir, "optimize-results")
        self.opt_name = opt_name
        self.natoms = natoms

        self.mappings = []
        self.map_filepath = ""
        self.probs_filepath = ""

        self._validate_inputs()
        self._ensure_required_files()  # Ensure required files exist

    def _validate_inputs(self):
        """Validate input arguments and ensure directory exists."""
        if not os.path.exists(self.directory):
            print(f"Error: Directory '{self.directory}' does not exist.")
            sys.exit(1)

    def _ensure_required_files(self):
        """Check if '48-MAPPINGS' and 'probabilities' files exist, create them if missing."""
        mapping_found = False
        probabilities_found = False

        for file in os.listdir(self.directory):
            if file.startswith("48-MAPPINGS"):
                self.map_filepath = os.path.join(self.directory, file)
                mapping_found = True
            elif file.startswith("probabilities"):
                self.probs_filepath = os.path.join(self.directory, file)
                probabilities_found = True

        # If files are missing, create empty ones
        if not mapping_found:
            self.map_filepath = os.path.join(self.directory, f"48-MAPPINGS_{self.opt_name}.txt")
            open(self.map_filepath, 'w').close()
            print(f"Created missing file: {self.map_filepath}")

        if not probabilities_found:
            self.probs_filepath = os.path.join(self.directory, f"probabilities_{self.opt_name}.txt")
            open(self.probs_filepath, 'w').close()
            print(f"Created missing file: {self.probs_filepath}")

    def make_counts(self, mapping_matrix, nmaps):
        """Memoization algorithm to count occurrences of each atom in each mapping."""
        memo_list = [0] * self.natoms
        for row in mapping_matrix:
            for value in row:
                memo_list[value] += 1
        return [x / nmaps for x in memo_list] 

    def process_files(self):
        """Iterate through directory files and process mappings."""
        for file in os.listdir(self.directory):
            print(file)
            filepath = os.path.join(self.directory, file)

            if file.startswith(self.opt_name):
                self._process_mapping_file(filepath)

            # Files were already checked and created if missing
            if len(self.mappings) == 48:
                self._write_results()
                break  # No need to continue processing more files

    def _process_mapping_file(self, filepath):
        """Extract the lowest mapping from an optimization file."""
        print("\nSA optimization file found. Writing lowest mapping into file...")

        with open(filepath, 'r') as f:
            lines = f.read()

        match = re.search(r'conv mapping\n([\d\s]+)last_smap', lines, re.DOTALL)
        if match:
            numbers = list(map(int, match.group(1).split()))
            self.mappings.append(numbers)
        else:
            print(f"Warning: No mapping found in {os.path.basename(filepath)}")

        print("Length of mapping file so far:", len(self.mappings))

    def _write_results(self):
        """Write mappings and probabilities to respective files."""
        print("\nWriting mapping into 48-MAPPINGS file...")
        with open(self.map_filepath, 'w') as mf:
            for row in self.mappings:
                mf.write(" ".join(map(str, row)) + "\n")  # Cleaner formatting

        print("\nCalculating probabilities...")
        probabilities = self.make_counts(self.mappings, len(self.mappings))

        with open(self.probs_filepath, 'w') as pf:
            pf.write(" ".join(map(str, probabilities)))

        print(f"\nResults written successfully to {self.map_filepath} and {self.probs_filepath}")

# ==============================
# Main execution
# ==============================
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Error: Missing arguments. Usage: python script.py <optimization_name> <num_atoms>")        
        sys.exit(1)

    opt_name = sys.argv[1]
    natoms = int(sys.argv[2])

    processor = MappingProcessor(opt_name, natoms)
    processor.process_files()
