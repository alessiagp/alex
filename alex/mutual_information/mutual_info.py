import os
import sys
import logging
import pickle as pkl
from mdigest.core.parsetrajectory import MDS  
from mdigest.core.correlation import DynCorr
from mdigest.core.dcorrelation import DihDynCorr
import mdigest.core.toolkit as toolkit
import configparser

import configparser

import configparser
import json
import yaml   # requires: pip install pyyaml

def parse_param_file(file_path):
    """
    Parse a parameter file (.ini, .yaml/.yml, .json).
    
    Returns:
        list of dict
    """
    if file_path.endswith(".ini"):
        return _parse_ini(file_path)
    elif file_path.endswith((".yaml", ".yml")):
        return _parse_yaml(file_path)
    elif file_path.endswith(".json"):
        return _parse_json(file_path)
    else:
        raise ValueError("Unsupported file format. Use .ini, .yaml/.yml, or .json")


def _convert_value(value):
    """Convert string to bool, int, float, or keep as string."""
    if isinstance(value, str):
        if value.lower() in ["true", "false"]:
            return value.lower() == "true"
        try:
            if "." in value:
                return float(value)
            return int(value)
        except ValueError:
            return value
    return value

def _parse_ini(file_path):
    config = configparser.ConfigParser()
    config.optionxform = str  # preserve case
    config.read(file_path)

    param_list = []
    for section in config.sections():
        params = {}
        for key, value in config[section].items():
            params[key] = _convert_value(value)
        param_list.append(params)
    return param_list

def _parse_yaml(file_path):
    with open(file_path, "r") as f:
        data = yaml.safe_load(f)
    # YAML can be a list of dicts or a single dict
    if isinstance(data, list):
        return data
    elif isinstance(data, dict):
        return [data]
    else:
        raise ValueError("YAML must define a list or dict of parameter sets.")


def _parse_json(file_path):
    with open(file_path, "r") as f:
        data = json.load(f)
    # JSON can be a list of dicts or a single dict
    if isinstance(data, list):
        return data
    elif isinstance(data, dict):
        return [data]
    else:
        raise ValueError("JSON must define a list or dict of parameter sets.")


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