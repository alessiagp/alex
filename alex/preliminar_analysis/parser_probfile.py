import sys
import os

workdir = f'{os.getcwd()}/optimize-results'

probabilities_name=str(sys.argv[1])
PROBS_name=str(sys.argv[2])

input_file = f"{workdir}/{probabilities_name}"
output_file = f"{workdir}/{PROBS_name}"

# Read the content from the input file
with open(input_file, 'r') as f:
	
	content = f.read()
	# Remove commas from the content
	content_without_commas = content.replace(', ', '\n')

	# Write the formatted content to the output file with each number on a new line
	with open(output_file, 'w') as f:
		f.write(content_without_commas)
		
		print("File successfully formatted!")

