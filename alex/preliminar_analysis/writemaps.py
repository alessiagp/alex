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

    def _validate_inputs(self):
        """
        Validate input arguments and ensure directory exists.
        """
        if not os.path.exists(self.directory):
            print(f"Error: Directory '{self.directory}' does not exist.")
            sys.exit(1)
    
    def make_counts(self, mapping_matrix, nmaps):
        """
        Memoization algorithm to count occurrences of each atom in each mapping.
        """
        memo_list = [0] * self.natoms
        for row in mapping_matrix:
            for value in row:
                memo_list[value] += 1
        return [x / nmaps for x in memo_list] 

    def process_files(self):
        """
        Iterate through directory files and process mappings.
        """
        for file in os.listdir(self.directory):
            print(file)
            filepath = os.path.join(self.directory, file)

            if file.startswith(self.opt_name):
                self._process_mapping_file(filepath)

            elif file.startswith("probabilities"):
                self.probs_filepath = filepath
                print("\nFound probabilities file:", self.probs_filepath)

            elif file.startswith("48-MAPPINGS"):
                self.map_filepath = filepath

            # When mappings reach 48, process and write results
            if len(self.mappings) == 48:
                self._write_results()
                break  # No need to continue processing more files

    def _process_mapping_file(self, filepath):
        """
        Extract the lowest mapping from an optimization file.
        """
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
        """
        Write mappings and probabilities to respective files.
        """
        if not self.map_filepath:
            print("Error: 48-MAPPINGS file not found.")
            sys.exit(1)

        print("\nWriting mapping into 48-MAPPINGS file...")
        with open(self.map_filepath, 'w') as mf:
            for row in self.mappings:
                mf.write(" ".join(map(str, row)) + "\n")  # Cleaner formatting

        print("\nCalculating probabilities...")
        probabilities = self.make_counts(self.mappings, self.natoms, len(self.mappings))

        if self.probs_filepath:
            with open(self.probs_filepath, 'w') as pf:
                pf.write(" ".join(map(str, probabilities)))
        else:
            print("Warning: No probabilities file found.")

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
