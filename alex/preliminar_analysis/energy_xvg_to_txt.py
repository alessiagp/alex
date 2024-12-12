import re
import sys

filename = str(sys.argv[1]) #name of energy file without extension
input_file = f'{filename}.xvg'
output_file = f'{filename}.txt'

numbers=[]
with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line and not line.startswith('#'):
            columns = line.split()
            if len(columns) >= 2:
                try:
                    number = float(columns[1])
                    numbers.append(number)
                except ValueError:
                    continue

with open(output_file, 'w') as f:
    for number in numbers:
        f.write(str(number) + '\n')

