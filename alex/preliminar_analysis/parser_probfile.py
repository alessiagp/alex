import os
import sys

class FileFormatter:
    def __init__(self, input_filename, output_filename):
        """
        Initialize the file formatter with input and output file paths.
        """
        self.workdir = os.path.join(os.getcwd(), "optimize-results")
        self.input_file = os.path.join(self.workdir, input_filename)
        self.output_file = os.path.join(self.workdir, output_filename)

        self._validate_input_file()

    def _validate_input_file(self):
        """
        Check if the input file exists.
        """
        if not os.path.exists(self.input_file):
            print(f"Error: Input file '{self.input_file}' not found.")
            sys.exit(1)

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
# Main execution
# ==============================
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Error: Missing arguments. Usage: python script.py <input_filename> <output_filename>")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    formatter = FileFormatter(input_filename, output_filename)
    formatter.format_file()
