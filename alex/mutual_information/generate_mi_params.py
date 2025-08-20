import argparse
import json
import yaml
import re

# Allowed options
VALID_CORR_TYPES = ["CA", "dih"]
VALID_LMI = ["gaussian"]

def validate_mi(value: str) -> bool:
    """Return True if value is a valid MI method string"""
    if value == "gaussian":
        return True
    # Matches knn_i_j where i and j are integers
    return bool(re.fullmatch(r"knn_\d+_\d+", value))

def ask_param(prompt, default=None, type_=str, choices=None, validator=None):
    """Helper: ask user with default, type conversion, and validation"""
    while True:
        if default is not None:
            q = f"{prompt} [{default}]: "
        else:
            q = f"{prompt}: "
        ans = input(q).strip()

        # Default
        if ans == "" and default is not None:
            return default

        # Type conversion
        try:
            if type_ == bool:
                val = ans.lower() in ["true", "yes", "1", "y"]
            else:
                val = type_(ans)
        except Exception:
            print(f"Invalid type. Expected {type_.__name__}")
            continue

        # Choice validation
        if choices and val not in choices:
            print(f"Invalid choice. Allowed: {choices}")
            continue

        # Custom validator
        if validator and not validator(val):
            print(f"Invalid value: {val}")
            continue

        return val


def interactive_param_set():
    """Ask user for one validated parameter set"""
    print("\n--- Define a new parameter set ---")
    params = {}
    params["correlation_type"] = ask_param(
        "Correlation type", default="CA", choices=VALID_CORR_TYPES
    )
    params["scale"] = ask_param("Scale displacements? (True/False)", True, bool)
    params["normalize"] = ask_param("Normalize data? (True/False)", True, bool)
    params["mean_center"] = ask_param("Mean-center (for dihedrals)?", True, bool)
    params["LMI"] = ask_param("LMI method", "gaussian", choices=VALID_LMI)
    params["MI"] = ask_param("MI method (gaussian or knn_i_j)", default="knn_5_1", validator=validate_mi)
    params["DCC"] = ask_param("Compute DCC? (True/False)", True, bool)
    params["PCC"] = ask_param("Compute PCC? (True/False)", True, bool)
    params["COV_DISP"] = ask_param("Compute covariance/dispersion? (True/False)", True, bool)
    return params

def save_params(param_list, fmt, filename):
    """Save parameter sets to chosen format"""
    if fmt == "ini":
        with open(filename, "w") as f:
            for i, params in enumerate(param_list, start=1):
                f.write(f"[set{i}]\n")
                for k, v in params.items():
                    f.write(f"{k} = {v}\n")
                f.write("\n")
    elif fmt == "txt":
        with open(filename, "w") as f:
            for i, params in enumerate(param_list, start=1):
                for k, v in params.items():
                    f.write(f"{k} = {v}\n")
                if i < len(param_list):
                    f.write("---\n")
    elif fmt == "json":
        with open(filename, "w") as f:
            json.dump(param_list, f, indent=2)
    elif fmt in ("yaml", "yml"):
        with open(filename, "w") as f:
            yaml.safe_dump(param_list, f)
    else:
        raise ValueError("Unsupported format")
    print(f"\nSaved {len(param_list)} parameter sets to {filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interactively generate parameter files for MutualInformation.")
    parser.add_argument("--format", choices=["ini", "txt", "json", "yaml"], default="ini", help="Output file format")
    parser.add_argument("--output", default=None, help="Output filename")
    args = parser.parse_args()

    param_list = []
    while True:
        param_list.append(interactive_param_set())
        more = input("Add another parameter set? (y/n): ").strip().lower()
        if more not in ["y", "yes"]:
            break

    filename = args.output or f"params.{args.format}"
    save_params(param_list, args.format, filename)
