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
    Simple algorithm to retrieve the occurrences of each amino acid in each mapping 
    according to the most retained atom in the amino acid.

    Takes as input the processed probability file of length equal to number of atoms, the number of atoms in the structure, and the protein dictionary assigning aminoacids to the delimiting atomnums. 
    
    1. use files from atomistic occurrence probabilities, iterate through them and select the aa according to AA_dict
    2. keep only the highest atom probability among all amino acids --> they are going to be ordered already
    """

    aa_counts = []
 
    if len(atomistic_probs) != natoms:
        raise ValueError("The number of probabilities in the file is different than the number of atoms")
    else:
        # select aa ranges
        for amino_acid, atom_interval in atom_to_aa_dict.items():
        # select data corresponding to ranges in prob file --> slice list
            temp_aa = atomistic_probs[atom_interval[0]:atom_interval[1]]
        # find max among these and append to counts
            aa_counts.append(max(temp_aa))

    return aa_counts
