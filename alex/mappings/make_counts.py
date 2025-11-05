import numpy as np
from collections import OrderedDict
    
#------------------------------------------------------------------------------------------------#

def atom_average_counts(mapping_matrix, nmaps, atom_to_aa_dict):
    """
    Memoisation algorithm to retrieve the occurrences of each amino acid in each mapping by calculating the average occurrence of atoms in each aminoacid.
    Takes as input the mapping matrix obtained by processing data with the `MAPPING` flag, the total number of mappings, and the protein dictionary assigning aminoacids to the delimiting atomnums.
    """
    
    # Create an ordered dictionary to store the counts for each amino acid
    amino_acid_counts = OrderedDict()
 
    # Flatten the NumPy array and iterate through the values
    for atom in mapping_matrix.flat:
        for amino_acid, atom_interval in atom_to_aa_dict.items():
            if atom_interval[0] <= atom < atom_interval[1]:
                amino_acid_counts[amino_acid] = amino_acid_counts.get(amino_acid, 0) + 1
    
    # Manually sort the probabilities based on the order of amino acids in AA_dict
    amino_acid_probabilities = OrderedDict()
    for amino_acid, interval in atom_to_aa_dict.items():
        count = amino_acid_counts.get(amino_acid, 0)
        probability = count / len(mapping_matrix.flat) #counts divided by tot number of mappings
        amino_acid_probabilities[amino_acid] = probability
    
    return amino_acid_probabilities
    
#------------------------------------------------------------------------------------------------#

def prob_maxatom(atomistic_probs, natoms, atom_to_aa_dict):
    """
    For each amino acid, return the maximum atomic probability inside its atom interval.
    
    Assumes:
    - atom_to_aa_dict ranges are 1-based and inclusive
    - ranges cover all atoms exactly once (checked automatically)
    """

    if len(atomistic_probs) != natoms:
        raise ValueError(
            f"The number of probabilities ({len(atomistic_probs)}) does not match natoms ({natoms})."
        )

    covered_atoms = []
    for start, end in atom_to_aa_dict.values():
        covered_atoms.extend(range(start, end + 1))

    covered_atoms = sorted(covered_atoms)

    if covered_atoms[0] != 1 or covered_atoms[-1] != natoms:
        raise ValueError(
            f"Dictionary does not cover all atoms (expected 1–{natoms}, "
            f"found {covered_atoms[0]}–{covered_atoms[-1]})."
        )

    for i in range(1, len(covered_atoms)):
        if covered_atoms[i] != covered_atoms[i - 1] + 1:
            raise ValueError(
                f"Gap or overlap detected between atoms {covered_atoms[i - 1]} and {covered_atoms[i]}."
            )

    aa_counts = []
    for aa, (start, end) in atom_to_aa_dict.items():
        temp = atomistic_probs[start - 1:end]  # convert 1-based to 0-based
        aa_counts.append(max(temp))

    return aa_counts
