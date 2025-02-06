import os
import re

class SmapExtractor:
    def __init__(self, directory):
        """
        Initialize the extractor with the target directory.
        """
        self.directory = os.path.join(os.getcwd(), directory)
        self.pattern = r"SMAP=(\d+\.\d+)"
        self.all_smap_values = []
        self.entropy_file = None  # Will be assigned later

        self._validate_directory()

    def _validate_directory(self):
        """
        Check if the directory exists.
        """
        if not os.path.exists(self.directory):
            print(f"Error: Directory '{self.directory}' does not exist.")
            exit(1)

    def extract_smap_values(self):
        """
        Extract SMAP values from .dat files in the directory.
        """
        for filename in os.listdir(self.directory):
            file_path = os.path.join(self.directory, filename)

            # Find the entropy file for later writing
            if filename.startswith("ENTROPY"):
                self.entropy_file = file_path

            # Extract SMAP values from .dat files
            if filename.endswith(".dat"):
                with open(file_path, "r") as file:
                    text = file.read()
                    smap_values = re.findall(self.pattern, text)

                if smap_values:
                    self.all_smap_values.append(smap_values)

        print("Extracted SMAP Values:", self.all_smap_values)
        print("Number of entropies retrieved:", len(self.all_smap_values))

    def write_to_entropy_file(self):
        """
        Write extracted SMAP values to the entropy file.
        """
        if not self.entropy_file:
            print("Error: No ENTROPY file found in the directory.")
            return

        with open(self.entropy_file, "w") as ef:
            for series in self.all_smap_values:
                ef.write(" ".join(series) + "\n")

        print(f"SMAP values successfully written to '{self.entropy_file}'")

# ==============================
# Main execution
# ==============================
if __name__ == "__main__":
    extractor = SmapExtractor("optimize-results")
    extractor.extract_smap_values()
    extractor.write_to_entropy_file()
