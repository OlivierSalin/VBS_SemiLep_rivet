import itertools
import glob
import os
import numpy as np
import csv


processes = [
    "WpZ_llqq", "WmZ_llqq", "ZZjj_llqq", "WpZ_vvqq", "WmZ_vvqq", "ZZ_vvqq",
    "WpZ_lvqq", "WmZ_lvqq", "WpWm_lvqq", "WmWm_lvqq", "WpWp_lvqq"
]
#processes= ["WpZ_llqq","ZZ_llqq"]
operators = ["FM0", "FS0", "FT0"]

# Output CSV file
output_file = "aQGC_combinations.csv"

# Generate combinations and write to CSV
with open(output_file, mode="w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    
    for process in processes:
        for operator in operators:
            # Extract the operator and process from the combination
            row = f"MGPy8EG_aQGC{operator}_QUAD_1_{process}"
            writer.writerow([row])

print(f"CSV file '{output_file}' created successfully.")