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
    elif file_path.endswith(".txt"):
        return _parse_txt(file_path)
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
