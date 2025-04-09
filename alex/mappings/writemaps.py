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
        self._ensure_required_files()

    def _validate_inputs(self):
        """Validate input arguments and ensure directory exists."""
        if not os.path.exists(self.directory):
            print(f"Error: Directory '{self.directory}' does not exist.")
            sys.exit(1)

    def _ensure_required_files(self):
        """Check if '48-MAPPINGS' and 'probabilities' files exist, create them if missing."""
        mapping_files = [f for f in os.listdir(self.directory) if f.startswith("48-MAPPINGS")]
        probability_files = [f for f in os.listdir(self.directory) if f.startswith("probabilities")]

        self.map_filepath = os.path.join(self.directory, mapping_files[0]) if mapping_files else os.path.join(self.directory, f"48-MAPPINGS_{self.opt_name}.txt")
        self.probs_filepath = os.path.join(self.directory, probability_files[0]) if probability_files else os.path.join(self.directory, f"probabilities_{self.opt_name}.txt")

        # Create empty files if missing
        if not mapping_files:
            open(self.map_filepath, 'w').close()
            print(f"Created missing file: {self.map_filepath}")
        if not probability_files:
            open(self.probs_filepath, 'w').close()
            print(f"Created missing file: {self.probs_filepath}")

    def make_counts(self, mapping_matrix, nmaps):
        """Memoization algorithm to count occurrences of each atom in each mapping."""
        if nmaps == 0:
            print("Warning: No mappings found. Returning zeroed probabilities.")
            return [0] * self.natoms

        memo_list = [0] * self.natoms
        for row in mapping_matrix:
            for value in row:
                memo_list[value] += 1
        return [x / nmaps for x in memo_list]

    def process_files(self):
        """Iterate through directory files and process mappings."""
        for file in os.listdir(self.directory):
            filepath = os.path.join(self.directory, file)
            if file.startswith(self.opt_name):
                self._process_mapping_file(filepath)
            
            if len(self.mappings) >= 48:
                self._write_results()

    def _process_mapping_file(self, filepath):
        """Extract the lowest mapping from an optimization file."""
        print("\nSA optimization file found. Writing lowest mapping into file...")
        
        with open(filepath, 'r') as f:
            lines = f.read()

        match = re.search(r'conv mapping\s*\n([\d\s]+?)\s*last_smap', lines, re.DOTALL)
        if match:
            numbers = list(map(int, match.group(1).split()))
            print(f"Extracted mapping: {numbers}")  # Debug print
            self.mappings.append(numbers)
        else:
            print(f"Warning: No mapping found in {os.path.basename(filepath)}")

        print("Length of mapping file so far:", len(self.mappings))

    def _write_results(self):
        """Write mappings and probabilities to respective files."""
        print("\nWriting mapping into 48-MAPPINGS file...")
        if not self.mappings:
            print("Warning: No mappings found. Writing empty file.")

        with open(self.map_filepath, 'w') as mf:  # Overwrite (instead of append)
            for row in self.mappings:
                mf.write(" ".join(map(str, row)) + "\n")

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
        print("Error: Missing arguments. Usage: python3 writemaps.py <optimization_name> <num_atoms>")        
        sys.exit(1)

    opt_name = sys.argv[1]
    natoms = int(sys.argv[2])

    processor = MappingProcessor(opt_name, natoms)
    processor.process_files()
