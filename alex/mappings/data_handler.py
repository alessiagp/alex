import os
import numpy as np
from collections import OrderedDict
    
#------------------------------------------------------------------------------------------------#

def load_data(workdir, filename):
    """
    Function that takes the working directory and the filename as input.
    The two accepted filenames indicate the formatting that will be produced:
     - `filename=PROBS` to have data formatted to calculate atomistic occurrence probability (plain data plot) and to calculate aminoacid occurrence probability by max atom probability;
     - `filename=48-MAPPINGS` to have data formatted to calculate aminoacid occurrence probability by average probability of atoms in the aminoacid.
    """
# generalise from 48-MAPPINGS to just MAPPINGS
    data_list = []
    for FILE in os.listdir(workdir):
        if filename == 'PROBS':
            if FILE.startswith(filename):
                filepath = os.path.join(workdir, FILE)
                data = np.loadtxt(filepath)
                data_list.append(data)
        elif filename == '48-MAPPINGS':
            if FILE.startswith(filename):
                filepath = os.path.join(workdir, FILE)
                data = np.loadtxt(filepath, delimiter=',')
                data_list.append(data.astype(int))
        else:
            raise ValueError(f"Filename {filename} is not in correct format")
    return data_list
    
#------------------------------------------------------------------------------------------------#

