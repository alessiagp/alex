import pandas as pd
import numpy as np
    
#------------------------------------------------------------------------------------------------#

def mem_parser(filename):
    """
    Function to extract and format memory data consumed by EXCOGITO in each process. 
    Takes as input the filename and returns a list with formatted values.
    """
    mem_list=[]
    with open(filename, 'r') as f:
        lines=f.readlines()
        for line in lines:
            line=str(line).strip('\n').rstrip('g')
            mem_list.append(float(line))
    return mem_list
    
#------------------------------------------------------------------------------------------------#

def extract_entropies(traj_entropy_file):
    """
    Function taking as input the entropy file, opens it, and ultimately formats it into a Pandas dataframe with columns `mean`, `std`, and `MC steps`.
    """
    E=[]

    with open(traj_entropy_file, 'r') as f:
        text = f.readlines()
        for row in text:
            row_proc=row.split(" ")
            row_proc.remove("\n")
            E.append(row_proc)
    E=np.float_(E)
    
    E=pd.DataFrame(E)
    E.astype(float)
    E=E.T #transposed
    
    E['mean'] = E.mean(axis=1)
    E['std'] = E.std(axis=1)
    E['MC steps']=E.index+1
    
    return E
    
#------------------------------------------------------------------------------------------------#

def entropy_data_for_plotting(extracted_entropy, confidence_interval):
    """
    Function preparing entropy data for plotting, takes as input the file processed by the `extract_entropy` function and the desired confidence interval.
    
    Returns a tuple with (xaxis, lower_bound, upper_bound)
    """
    xaxis = range(0, len(extracted_entropy['mean']))
    
    lower_bound = np.asarray(extracted_entropy['mean']-(np.asarray(extracted_entropy['std']*confidence_interval)))
    
    upper_bound = np.asarray(extracted_entropy['mean']+(np.asarray(extracted_entropy['std']*confidence_interval)))
    
    return xaxis, lower_bound, upper_bound

#------------------------------------------------------------------------------------------------#
    
def min_position(lst):
    """
    Function to retrieve the index of the minimum of a list.
    """
    arr = np.array(lst)
    min_pos = np.argmin(arr) + 1
    return min_pos
