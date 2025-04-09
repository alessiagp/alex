from pathlib import Path
import re
import argparse
from typing import Optional

class SmapExtractor:
    def __init__(self, directory: str, entropy_suffix: Optional[str] = None):
        """
        Initialize the extractor with the target directory and optional entropy output suffix.
        """
        self.directory = Path.cwd() / directory
        self.pattern = re.compile(r"SMAP=(\d+\.\d+)")
        self.all_smap_values = []
        self.entropy_file: Path | None = None
        self.entropy_suffix = entropy_suffix

        self._validate_directory()

    def _validate_directory(self):
        """
        Check if the directory exists.
        """
        if not self.directory.exists():
            raise FileNotFoundError(f"Directory '{self.directory}' does not exist.")

    def extract_smap_values(self):
        """
        Extract SMAP values from .dat files in the directory.
        """
        for file_path in self.directory.iterdir():
            if file_path.is_file():
                if file_path.name.startswith("ENTROPY") and self.entropy_file is None:
                    self.entropy_file = file_path

                if file_path.suffix == ".dat":
                    try:
                        text = file_path.read_text()
                        smap_values = self.pattern.findall(text)
                        if smap_values:
                            self.all_smap_values.append(smap_values)
                    except Exception as e:
                        print(f"Error reading file {file_path.name}: {e}")

        print("Extracted SMAP Values:", self.all_smap_values)
        print("Number of entropies retrieved:", len(self.all_smap_values))

    def write_to_entropy_file(self):
        """
        Write extracted SMAP values to the entropy file.
        Creates a new entropy file with optional custom name if none exists.
        """
        if not self.entropy_file:
            filename = (
                f"ENTROPY_{self.entropy_suffix}.txt"
                if self.entropy_suffix
                else "ENTROPY_OUTPUT.txt"
            )
            self.entropy_file = self.directory / filename
            print(f"No ENTROPY file found. Creating new file: {self.entropy_file.name}")

        try:
            with self.entropy_file.open("w") as ef:
                for series in self.all_smap_values:
                    ef.write(" ".join(series) + "\n")
            print(f"SMAP values successfully written to '{self.entropy_file}'")
        except Exception as e:
            print(f"Error writing to file {self.entropy_file}: {e}")


    def main():
        parser = argparse.ArgumentParser(description="Extract SMAP values and write to an ENTROPY file.")
        parser.add_argument("filename", help="Suffix to create ENTROPY_<filename>.txt")
        parser.add_argument(
            "-d", "--directory", default="optimize-results",
            help="Directory to scan for .dat and ENTROPY files (default: optimize-results)"
        )
        args = parser.parse_args()

        extractor = SmapExtractor(directory=args.directory, entropy_suffix=args.filename)
        extractor.extract_smap_values()
        extractor.write_to_entropy_file()

    if __name__ == "__main__":
        main()
