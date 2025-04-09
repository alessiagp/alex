import argparse
import re
import logging
from pathlib import Path
from typing import Optional

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

class SmapExtractor:
    def __init__(self, directory: str, entropy_suffix: Optional[str] = None):
        self.directory = Path.cwd() / directory
        self.pattern = re.compile(r"SMAP=(\d+\.\d+)")
        self.all_smap_values = []
        self.entropy_file: Path | None = None
        self.entropy_suffix = entropy_suffix

        self._validate_directory()

    def _validate_directory(self):
        if not self.directory.exists():
            raise FileNotFoundError(f"Directory '{self.directory}' does not exist.")

    def extract_smap_values(self):
        for file_path in self.directory.iterdir():
            if file_path.is_file():
                if file_path.name.startswith("ENTROPY") and self.entropy_file is None:
                    self.entropy_file = file_path

                if file_path.suffix == ".dat" and file_path.name.startswith(self.entropy_suffix or ""):
                    try:
                        text = file_path.read_text()
                        smap_values = self.pattern.findall(text)
                        if smap_values:
                            self.all_smap_values.append(smap_values)
                    except Exception as e:
                        logging.error(f"Error reading file {file_path.name}: {e}")

        logging.info(f"Extracted SMAP Values: {self.all_smap_values}")
        logging.info(f"Number of entropy series retrieved: {len(self.all_smap_values)}")

    def write_to_entropy_file(self):
        if not self.entropy_file:
            filename = (
                f"ENTROPY_{self.entropy_suffix}.txt"
                if self.entropy_suffix
                else "ENTROPY_OUTPUT.txt"
            )
            self.entropy_file = self.directory / filename
            logging.info(f"No ENTROPY file found. Creating new file: {self.entropy_file.name}")

        try:
            with self.entropy_file.open("w") as ef:
                for series in self.all_smap_values:
                    ef.write(" ".join(series) + "\n")
            logging.info(f"SMAP values successfully written to '{self.entropy_file}'")
        except Exception as e:
            logging.error(f"Error writing to file {self.entropy_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Extract SMAP values and write to an ENTROPY file.")
    parser.add_argument("filename", help="Suffix to create ENTROPY_<filename>.txt and filter .dat files")
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
