import os
import sys

class FileFormatter:
    def __init__(self, opt_name):
        """
        Initialize the file formatter with input and output file paths.
        """
        self.workdir = os.path.join(os.getcwd(), "optimize-results")
        self.opt_name = opt_name
        self.input_file = os.path.join(self.workdir, f"probabilities_{self.opt_name}.txt")
        self.output_file = os.path.join(self.workdir, f"PROBS_{self.opt_name}.txt")

        self._validate_input_file()
        self._ensure_output_file()  # Ensure output file exists

    def _validate_input_file(self):
        """
        Check if the input file exists.
        """
        if not os.path.exists(self.input_file):
            print(f"Error: Input file '{self.input_file}' not found.")
            sys.exit(1)

    def _ensure_output_file(self):
        """
        Create the output file if it does not exist.
        """
        if not os.path.exists(self.output_file):
            open(self.output_file, 'w').close()
            print(f"Created missing output file: {self.output_file}")

    def format_file(self):
        """
        Reads input file, replaces commas with newlines, and writes to output file.
        """
        with open(self.input_file, 'r') as f:
            content = f.read()

        # Replace commas with newlines
        formatted_content = content.replace(',', '\n')

        with open(self.output_file, 'w') as f:
            f.write(formatted_content)

        print(f"File successfully formatted and saved as '{self.output_file}'")

# ==============================
# Main Execution
# ==============================
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: Missing argument. Usage: python script.py <optimization_name>")
        sys.exit(1)

    opt_name = sys.argv[1]  # Now correctly uses only opt_name
    formatter = FileFormatter(opt_name)
    formatter.format_file()
