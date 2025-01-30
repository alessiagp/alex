import os
import re
import sys

mappings=[]
entropies = []
workdir = os.getcwd()
directory=f"{workdir}/optimize-results/"
map_filepath=""
probs_filepath=""

opt_name=str(sys.argv[1])

def make_counts(mapping_matrix, natoms, nmaps): #memoisation algorithm
	"""Memoisation algorithm to retrieve the occurrences of each atom in each mapping"""	
	memo_list = [0] * natoms 
	for row in mapping_matrix: #for each of the 48 mappings
		for value in row: #for each atom in the mapping
			memo_list[value]+=1
	memo_list[:] = [x / nmaps for x in memo_list]
	return memo_list

for FILE in os.listdir(directory):
    print(FILE)

    """Retrieve LOWEST MAPPING from each of the 48 files"""
    if FILE.startswith(opt_name):
        SA_filepath = os.path.join(directory, FILE)
        print("\nSA optimisation file found. Writing lowest mapping into file...")
        with open(SA_filepath, 'r') as f:
            lines = f.read() #.replace('\n',"")
            numbers_string = re.search(r'conv mapping\n([\d\s]+)last_smap', lines, re.DOTALL).group(1)
            numbers = list(map(int, numbers_string.split()))
            mappings.append(numbers)

    print("Length of mapping file so far: ", len(mappings))

    """Calculate probabilities for all atoms in the protein and copy results into file"""
    if FILE.startswith('probabilities'):
        probs_filepath = os.path.join(directory, FILE)
        print("\nfound probabilities file!")

    """Recording mapping file to global variable"""
    if FILE.startswith('48-MAPPINGS'):
       map_filepath = os.path.join(directory, FILE)

    """ Write mapping of current trajectory into file"""
    if len(mappings) == 48:
        print("\nwriting mapping into 48-MAPPINGS file...")
        with open(map_filepath, 'w') as MF:
            for row in mappings:
                MF.write((str(row).strip("[]"))+"\n")
        print("\ncalculating probabilities...")
        probabilities = make_counts(mappings, 93, len(mappings))
        with open(probs_filepath, 'w') as PF:
            PF.write(str(probabilities).strip("[]"))


