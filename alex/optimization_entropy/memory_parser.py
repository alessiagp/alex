import os
import numpy as np
import pandas as pd

class EntropyAnalyzer:
    def __init__(self, entropy_filename, confidence_z=1.96):
        """
        Initialize the analyzer with file paths and confidence interval Z-score.
        """
        self.entropy_filename = entropy_filename
        self.confidence_z = confidence_z
        self.entropy_data = None

        self._validate_entropy_file()

    def _validate_entropy_file(self):
        """
        Check if entropy input file exist.
        """
        if not os.path.exists(self.entropy_filename):
            raise FileNotFoundError(f"Error: File '{self.entropy_filename}' not found.")

    def extract_entropies(self):
        """
        Extract entropy values from the file and format them into a DataFrame.
        """
        E = []
        with open(self.entropy_filename, 'r') as f:
            for row in f:
                row_proc = row.strip().split(" ")
                E.append(row_proc)

        E = np.float64(E)  # Convert to NumPy float array
        df = pd.DataFrame(E).astype(float).T  # Transpose & ensure float type

        # Add statistical columns
        df['mean'] = df.mean(axis=1)
        df['std'] = df.std(axis=1)
        df['MC steps'] = df.index + 1  # Monte Carlo steps (starting at 1)

        self.entropy_data = df  # Store in class for further use
        return df

    def prepare_plot_data(self):
        """
        Prepare entropy data for plotting with confidence intervals.
        """
        if self.entropy_data is None:
            raise ValueError("Entropy data has not been extracted yet. Call `extract_entropies()` first.")

        mean = self.entropy_data['mean']
        std = self.entropy_data['std']
        xaxis = range(len(mean))

        margin_error = self.confidence_z * (std / np.sqrt(len(mean)))
        lower_bound = mean - margin_error
        upper_bound = mean + margin_error

        self.xaxis = xaxis

        return xaxis, lower_bound.to_numpy(), upper_bound.to_numpy()
    
    def entropy_assessment(self):
        """
        Assess the goodness of the optimization by checking the ratio of the differences in absolute value between the value of mapping entropy
        at half the optimization and the end, and the value of mapping entropy at half and the beginning.
        D = |Smap_1/2 - Smap_f|/|Smap_1/2 - Smap_i|
        """
        Smap_half = self.entropy_data['mean'][(len(self.xaxis)//2)-1]
        Smap_i = self.entropy_data['mean'][0]
        Smap_f = self.entropy_data['mean'][len(self.xaxis)-1]

        ratio = abs(Smap_half - Smap_f)/abs(Smap_half - Smap_i)
        print(f"\nRatio of differences: {ratio}")


    @staticmethod
    def min_position(lst):
        """
        Retrieve the index of the minimum value in a list (1-based index).
        """
        return np.argmin(np.array(lst)) + 1

    def run(self):
        """
        Execute the full analysis for entropy.
        Returns the entropy dataframe, the xaxis, the lower, and the upper bound.
        """

        entropy_df = self.extract_entropies()
        print("\nExtracted Entropy Data:\n", entropy_df.head())

        xaxis, lower, upper = self.prepare_plot_data()
        print("\nPrepared Plot Data:", xaxis, lower, upper)

        min_pos = self.min_position(entropy_df['mean'])
        print("\nMinimum entropy position:", min_pos)

        self.entropy_assessment()

        return entropy_df, xaxis, lower, upper

class MemoryAnalyzer:
    def __init__(self, mem_filename):
        """
        Initialize the analyzer with file paths and confidence interval Z-score.
        """
        self.mem_filename = mem_filename
        self._validate_memory_files()

    def _validate_memory_file(self):
        """
        Check if memory input file exist.
        """
        if not os.path.exists(self.mem_filename):
            raise FileNotFoundError(f"Error: File '{self.mem_filename}' not found.")

    def parse_memory(self):
        """
        Extract and format memory data consumed by EXCOGITO.
        """
        mem_list = []
        with open(self.mem_filename, 'r') as f:
            for line in f:
                formatted_value = line.strip().rstrip('g')  # Remove newline & 'g' suffix
                mem_list.append(float(formatted_value))
        return mem_list
    
    def run(self):
        """
        Execute the full analysis for the consumed memory.
        Returns the memory data.
        """
        memory_data = self.parse_memory()
        print("Memory Data:", memory_data)
        return memory_data
    
