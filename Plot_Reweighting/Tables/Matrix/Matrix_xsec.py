import pandas as pd
import numpy as np
import os
import re
import itertools
import seaborn as sns
import matplotlib.pyplot as plt

# ==================== USER CONFIGURATION ====================

all_ops_cat = [
    "FM0", "FM1", "FM2", "FM3", "FM4", "FM5", "FM7", "FM8", "FM9", 
    "FS0", "FS1", "FS2",
    "FT0", "FT1", "FT2", "FT3", "FT4", "FT5", "FT6"
]
cross_terms = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2) if op1[:2] == op2[:2]]
Type_MC_all = [
    "EFTDec_Madspin", "EFTDec_polarisation",
    "Rwg_pol_50k", "Rwg_pol_100k", 
    "Rwg_InvSqrtXsec_50k", "Rwg_InvXsec_50k"
]
groups = ["FM", "FS", "FT"]

path_base = "/exp/atlas/salin/ATLAS/VBS_mc/plotting/Plot_Reweighting/Tables/Tables/Uncertainty/Xsec"

# ==================== FUNCTIONS ====================

def parse_file(filepath):
    data = {}
    with open(filepath, "r") as f:
        for line in f:
            match = re.match(r"(.*?):\s*(-?\d*\.?\d+(?:[eE][-+]?\d+)?)", line.strip())
            if match:
                key, value = match.group(1), float(match.group(2))
                data[key] = value
    return data

def build_matrix(data, op_list, process_name="WpZ_llqq"):
    n = len(op_list)
    mat = pd.DataFrame(index=op_list, columns=op_list, dtype=float)
    for i, op1 in enumerate(op_list):
        for j, op2 in enumerate(op_list):
            if i == j:
                key = f"{process_name}_{op1}_QUAD"
            elif i < j:
                key = f"{process_name}_{op1}vs{op2}_CROSS"
            else:
                key = f"{process_name}_{op2}vs{op1}_CROSS"
            mat.loc[op1, op2] = data.get(key, np.nan)
    return mat

def plot_matrix(matrix, outname, title, cmap="YlGnBu", vmin=None, vmax=None):
    plt.figure(figsize=(8, 6))
    sns.heatmap(matrix, annot=True, fmt=".2g", cmap=cmap, vmin=vmin, vmax=vmax,
                linewidths=0.5, square=True, cbar_kws={"label": "σ [fb]"})
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outname + ".png")
    plt.close()

# ==================== MAIN EXECUTION ====================

for type_mc in Type_MC_all:
    lower_type_mc = type_mc.lower()
    file_xsec = f"{path_base}/VBS_cross_section_run2_{lower_type_mc}.txt"
    file_unc = f"{path_base}/VBS_Uncertainty_cross_section_run2_{lower_type_mc}.txt"

    if not os.path.exists(file_xsec):
        print(f"WARNING: Missing file {file_xsec}")
        continue

    data_xsec = parse_file(file_xsec)

    for group in groups:
        ops = [op for op in all_ops_cat if op.startswith(group)]
        matrix = build_matrix(data_xsec, ops)
        outname = f"matrix_{group}_{type_mc}"
        matrix.to_csv(outname + ".csv")
        plot_matrix(matrix, outname, f"{group} Cross Section Matrix ({type_mc})")

print("Done. CSV and PNG files generated.")