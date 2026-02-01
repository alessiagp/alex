import pathlib
import numpy as np
import pandas as pd

def extract_entropies(entropy_filename : pathlib.PosixPath):
    """
    Extract entropy values from the file and format them into a DataFrame.
    """
    E = []
    with entropy_filename.open() as f:
        for row in f:
            row_proc = row.strip().split(" ")
            E.append(row_proc)

    E = np.float64(E)
    df = pd.DataFrame(E).astype(float).T

    # Add statistical columns
    df['mean'] = df.mean(axis=1)
    df['std'] = df.std(axis=1)
    df['MC steps'] = df.index + 1  # Monte Carlo steps (starting at 1)

    return df

def plot_confidence_intervals(entropy_df, confidence_z=1.96):
    """
    Prepare entropy data for plotting with confidence intervals.
    """

    mean = entropy_df['mean']
    std = entropy_df['std']
    xaxis = range(len(mean))

    margin_error = confidence_z * (std / np.sqrt(len(mean)))
    lower_bound = mean - margin_error
    upper_bound = mean + margin_error

    return xaxis, lower_bound.to_numpy(), upper_bound.to_numpy()