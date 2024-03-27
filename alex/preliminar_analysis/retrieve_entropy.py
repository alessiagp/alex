import re
import os

# Directory where your files are located
curr_dir=os.getcwd()
directory = curr_dir+'/optimize-results'

# Define the regular expression pattern
pattern = r'SMAP=(\d+\.\d+)'

# List to store all SMAP values
all_smap_values = []

# Loop through the files in the directory
for filename in os.listdir(directory):
    Smap=[]
    if filename.endswith(".dat"):  # Change the file extension as needed
        file_path = os.path.join(directory, filename)
        with open(file_path, 'r') as file:
            text = file.read()
            smap_values = re.findall(pattern, text)
            Smap.extend(smap_values)
    if len(Smap) != 0:
        all_smap_values.append(Smap)

# Check for values retrieved
print(all_smap_values)
print("Number of entropies retrieved: ", len(all_smap_values))

# we write the matrix into the ENTROPY-1-1-100 file
for filename in os.listdir(directory):
    if filename.startswith("ENTROPY"):
        entropy_path = os.path.join(directory, filename)
        with open(entropy_path, 'w') as EF:
            for series in all_smap_values:
                for value in series:
                    EF.write(str(value).strip("''")+" ")
                EF.write("\n")
