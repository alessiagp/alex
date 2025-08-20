import argparse
import json
import yaml

def ask_param(prompt, default=None, type_=str):
    """Helper to ask user for input with default + type conversion"""
    if default is not None:
        prompt = f"{prompt} [{default}]: "
    else:
        prompt = f"{prompt}: "
    ans = input(prompt).strip()
    if ans == "" and default is not None:
        return default
    if type_ == bool:
        return ans.lower() in ["true", "yes", "1", "y"]
    return type_(ans)

def interactive_param_set():
    """Ask user for one parameter set interactively"""
    print("\n--- Define a new parameter set ---")
    params = {}
    params["correlation_type"] = ask_param("Correlation type (CA/dih)", "CA")
    params["scale"]            = ask_param("Scale displacements? (True/False)", True, bool)
    params["normalize"]        = ask_param("Normalize data? (True/False)", True, bool)
    params["mean_center"]      = ask_param("Mean-center (for dihedrals)?", True, bool)
    params["LMI"]              = ask_param("LMI method", "gaussian")
    params["MI"]               = ask_param("MI method", "knn_5_2")
    params["DCC"]              = ask_param("Compute DCC? (True/False)", True, bool)
    params["PCC"]              = ask_param("Compute PCC? (True/False)", True, bool)
    params["COV_DISP"]         = ask_param("Compute covariance/dispersion? (True/False)", True, bool)
    return params

def save_params(param_list, fmt, filename):
    """Save parameter sets to given format"""
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
    print(f"\n✅ Saved {len(param_list)} parameter sets to {filename}")

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
