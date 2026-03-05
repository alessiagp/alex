from collections import OrderedDict
import numpy as np

def parse_gro_to_dict(gro_file, keep_resname=True):
    """
    Parse a GROMACS .gro file and return an OrderedDict mapping
    residue (ResName+ResID) -> (start_atom_index, end_atom_index),
    where indices are global atom indices in the .gro file.
    """
    AA_dict = OrderedDict()
    
    with open(gro_file, "r") as f:
        lines = f.readlines()
    
    # First line = title, second line = atom count, last line = box
    atom_lines = lines[2:-1]
    
    prev_resid = None
    prev_resname = None
    start_idx = None
    
    for i, line in enumerate(atom_lines, start=1):  # start=1 for global atom index
        resid = int(line[0:5])
        resname = line[5:10].strip()
        
        if resid != prev_resid:
            # close previous residue block
            if prev_resid is not None and keep_resname==True:
                AA_dict[f"{prev_resname}{prev_resid}"] = (start_idx, i-1)
            elif prev_resid is not None and keep_resname==False:
                AA_dict[f"{prev_resid}"] = (start_idx, i-1)
            # start new block
            start_idx = i
            prev_resid = resid
            prev_resname = resname
    
    # close last residue
    if prev_resid is not None and keep_resname==True:
        AA_dict[f"{prev_resname}{prev_resid}"] = (start_idx, len(atom_lines))
    elif prev_resid is not None and keep_resname==False:
        AA_dict[f"{prev_resid}"] = (start_idx, i-1)

    return AA_dict