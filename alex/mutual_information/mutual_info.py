import os
import sys
import logging
import pickle as pkl
from mdigest.core.parsetrajectory import MDS  
from mdigest.core.correlation import DynCorr
from mdigest.core.dcorrelation import DihDynCorr
import mdigest.core.toolkit as toolkit

import configparser

from parse_mi_params import parse_param_file


def to_pickle(dataframe, output):
    """Save a dataframe to a pickle file."""
    pkl.dump(dataframe, open(output, 'wb'))

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class MutualInformation:
    def __init__(self, structure_file: str, trajectory_list: list, n_replicas: int, selection: str, stride: int, save_dir: str,
                 param_list: list = None):
        self.structure_file = structure_file
        self.trajectory_list = trajectory_list
        self.selection = selection
        self.n_replicas = n_replicas
        self.stride = stride
        self.save_dir = save_dir
        self.param_list = param_list or []   # <--- store parameter sets
        
        self.mds = MDS()
        
        self.filename = None
        self.correlation_type = None
        self.cache_dir = os.path.join(self.save_dir, "CACHE")
        self.results_dir = os.path.join(self.save_dir, "RESULTS", "correlations")

    def initialise_system(self, selection_1='protein and name CA', selection_2='protein', initial=0, final=-1): 
        try: 
            logging.info("Initializing the molecular system...") 
            self.mds.set_num_replicas(self.n_replicas) 
            self.mds.load_system(self.structure_file, self.trajectory_list) 
            self.mds.align_traj(inmem=True, selection=self.selection) 
            self.mds.set_selection(selection_1, selection_2) 
            self.mds.stride_trajectory(initial=initial, final=final, step=self.stride) 
            logging.info("System initialized successfully.") 
        except Exception as e: 
            logging.error(f"Failed to initialize the system: {e}") 
            raise

    def correlations(self):
        """
        Run correlations using all parameter sets provided in param_list.
        """
        if not self.param_list:
            raise ValueError("No parameter list provided to MutualInformation.")
        
        results = {}
        for i, params in enumerate(self.param_list, start=1):
            correlation_type = params.get("correlation_type", "CA").lower()
            self.correlation_type = correlation_type
            logging.info(f"Computing correlations set {i} with parameters: {params}")
            
            try:
                if correlation_type == 'ca':
                    dyncorr = DynCorr(self.mds)
                    dyncorr.parse_dynamics(**params)
                    results[i] = dyncorr
                elif correlation_type == 'dih':
                    dihdyncorr = DihDynCorr(self.mds)
                    dihdyncorr.parse_dih_dynamics(**params)
                    results[i] = dihdyncorr
                else:
                    raise ValueError("Invalid correlation_type: Use 'CA' or 'dih'.")
            except Exception as e:
                logging.error(f"Failed to compute correlations set {i}: {e}")
                raise
        return results


    def save_data(self, corr, filename=''):
        self.filename = filename
        toolkit.folder_exists(self.cache_dir)
        toolkit.folder_exists(self.results_dir)
        
        file_path = os.path.join(self.cache_dir, self.filename)
        try:
            corr.save_class(file_name_root=file_path)
            logging.info(f"Correlation data saved to {file_path}")
        except Exception as e:
            logging.error(f"Failed to save correlation data: {e}")
            raise

    def processing(self):
        logging.info("Initializing the system...")
        self.initialise_system()
        
        logging.info("Starting correlation calculations from parameter list...")
        results = self.correlations()
        
        # Save each result
        for i, corr in results.items():
            filename = f"corr_set_{i}_{self.correlation_type}.pkl"
            self.save_data(corr, filename=filename)
        
        logging.info("All correlations completed and saved.")
        return results


def parse_trajectory_list(input_data):
    if isinstance(input_data, str):
        trajectory_list = input_data.split(',')
    elif isinstance(input_data, list):
        trajectory_list = input_data
    else:
        raise TypeError("Invalid input type. Expected a string or a list.")
    
    trajectory_list = [traj.strip() for traj in trajectory_list]
    for traj in trajectory_list:
        if not os.path.isfile(traj):
            raise ValueError(f"Trajectory file '{traj}' does not exist.")
    return trajectory_list

if __name__ == '__main__':
    if len(sys.argv) < 7:
        raise ValueError("Usage: mutual_info.py <structure_file> <trajectory_list> <selection> <n_replicas> <stride> <save_dir> <param_file>")
    
    structure_file = sys.argv[1].strip()
    trajectory_list = parse_trajectory_list(sys.argv[2])
    selection = sys.argv[3].strip()
    n_replicas = int(sys.argv[4])
    stride = int(sys.argv[5])
    save_dir = sys.argv[6].strip()
    param_file = sys.argv[7].strip()

    # Parse param file
    param_list = parse_param_file(param_file)

    mi_exec = MutualInformation(
        structure_file=structure_file,
        trajectory_list=trajectory_list,
        selection=selection,
        n_replicas=n_replicas,
        stride=stride,
        save_dir=save_dir,
        param_list=param_list
    )

    mi_exec.processing()