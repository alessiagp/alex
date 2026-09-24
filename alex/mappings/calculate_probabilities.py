import sys
import numpy as np


def load_mapping_file(filepath):
    """
    Load mapping file where each row contains atom indices separated by spaces.
    """
    mappings = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line:
                mappings.append(list(map(int, line.split())))

    return np.array(mappings)


def compute_probabilities(mapping_matrix, natoms):
    """
    Compute probabilities exactly as in the original script,
    but using numpy for speed.
    Compute the 95% confidence interval error for atomic probabilities 
    derived from binary Bernoulli trials (optimizations).
    """

    nmaps = mapping_matrix.shape[0]

    if nmaps == 0:
        raise ValueError("No mappings found.")

    # Flatten all atom indices
    flat_atoms = mapping_matrix.flatten()

    # Count occurrences of each atom
    counts = np.bincount(flat_atoms, minlength=natoms)

    # Normalize by number of mappings
    probabilities = counts / nmaps
    
    # Calculate the standard error and 95% confidence interval
    variance = probabilities * (1.0 - probabilities)
    standard_error = np.sqrt(variance / nmaps)
    errors_95ci = 1.96 * standard_error

    return probabilities, errors_95ci

def main():

    if len(sys.argv) < 3:
        print("Usage: python compute_probabilities.py <mapping_file> <num_atoms> <output_prefix>")
        sys.exit(1)

    mapping_file = sys.argv[1]
    natoms = int(sys.argv[2])
    output_prefix = sys.argv[3]

    mapping_matrix = load_mapping_file(mapping_file)

    probabilities, errors = compute_probabilities(mapping_matrix, natoms)

    output_file_prob=f"{output_prefix}_probabilities.txt"
    output_file_err=f"{output_prefix}_errors.txt"
    
    with open(output_file_prob, "w") as f:
        f.write("\n".join(map(str, probabilities)))

    with open(output_file_err, "w") as f:
        f.write("\n".join(map(str, errors)))
        
    print(f"Probabilities written to {output_file_prob}")
    print(f"95% confidence errors written to {output_file_err}")

if __name__ == "__main__":
    main()