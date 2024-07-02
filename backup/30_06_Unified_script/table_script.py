import glob
import os
from array import array

import shutil
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
import pandas as pd   
from pandas.plotting import table

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Ana", default = "WpZ_llqq")
parser.add_option("--DOCUT", default = "YES")

opts, _ = parser.parse_args()

all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    
base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables"

output_table = f"{base_dir}/{opts.Ana}/Plots/Table/test01/DOCUT_{opts.DOCUT}/"

# Check if the directory exists, if not, create it
if not os.path.exists(output_table):
    os.makedirs(output_table)
    

    
Cutflow_paths = {
    "SM": f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_FM0_SM_EXT0/DOCUT_{opts.DOCUT}/",
    **{op: f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_{op}_QUAD_EXT0/DOCUT_{opts.DOCUT}/" for op in all_ops}
}
print(Cutflow_paths["SM"])
all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

Cutflow_path = {}

for op in all_ops_SM:
    cutflow_merged_path = Cutflow_paths[op] + "cutflow_merged.txt"
    cutflow_resolved_path = Cutflow_paths[op] + "cutflow_resolved.txt"
    cutflow_frac_res_path = Cutflow_paths[op] +"frac_after_cuts_error_bar_resolved.txt"
    cutflow_frac_merged_path = Cutflow_paths[op] +"frac_after_cuts_error_bar_merged.txt"

    if os.path.exists(cutflow_merged_path) and os.path.exists(cutflow_resolved_path):
        if os.path.exists(cutflow_frac_res_path) and os.path.exists(cutflow_frac_merged_path):
            Cutflow_path[op]= Cutflow_paths[op]
            Cutflow_path[op + "_merged"]= cutflow_merged_path
            Cutflow_path[op + "_resolved"]= cutflow_resolved_path
            Cutflow_path[op + "_frac_merged"]= cutflow_frac_merged_path 
            Cutflow_path[op + "_frac_res"]= cutflow_frac_res_path


def extract_counts(Cutflow_path, op):
    Info_table = {}
    Info_err_frac={}
    cutflow_merged_path = Cutflow_path[op + "_merged"]
    cutflow_resolved_path = Cutflow_path[op + "_resolved"]
    cutflow_frac_res_path = Cutflow_path[op + "_frac_res"]
    cutflow_frac_merged_path = Cutflow_path[op + "_frac_merged"]
    with open(cutflow_merged_path, 'r') as f:
        Lines = [line for line in f.readlines()]
        Total_line = Lines[2].split()
        Total_nb_events_merged = float(Total_line[0])
        Merged_line = Lines[-2].split()
        Merged_nb_events = float(Merged_line[2])

    with open(cutflow_resolved_path, 'r') as f:
        Lines = [line for line in f.readlines()]
        Total_line = Lines[2].split()
        Total_nb_events_resolved = float(Total_line[0])
        Resolved_line = Lines[-2].split()
        Resolved_nb_events = float(Resolved_line[2])
        
    with open(cutflow_frac_res_path, 'r') as f:
        Lines = [line for line in f.readlines()]
        err_frac_res = float(Lines[0].split()[0])
           
    with open(cutflow_frac_merged_path, 'r') as f:  
        Lines = [line for line in f.readlines()]
        err_frac_merged = float(Lines[0].split()[0])

    # Check if total number of events are the same in both files
    if Total_nb_events_merged == Total_nb_events_resolved:
        Info_table[op + "_total"] = Total_nb_events_merged
        Info_table[op + "_merged"] = Merged_nb_events
        Info_table[op + "_resolved"] = Resolved_nb_events
        Info_err_frac[op + "_err_frac_resolved"] = err_frac_res
        Info_err_frac[op + "_err_frac_merged"] = err_frac_merged
    else:
        print("Total number of events are not the same in both files.")
    return Info_table, Info_err_frac


def get_all_counts(Cutflow_path):
    all_counts = {}
    all_error={}
    for op in Cutflow_path.keys():
        if op.endswith("_merged") or op.endswith("_resolved") or op.endswith("_frac_res") or op.endswith("_frac_merged"):
            continue
        Info_table, Info_err_frac  = extract_counts(Cutflow_path, op)
        all_counts.update(Info_table)
        all_error.update(Info_err_frac)
    return all_counts, all_error

# Usage
all_counts,all_error = get_all_counts(Cutflow_path)



# Convert the dictionary to a DataFrame


def format_value(row):
    if row['Signal region'] == 'merged' or row['Signal region'] == 'resolved':
        percentage = round(row[opts.Ana] / common_total_value * 100, 2)
        
        error_key = f"{row['Operator']}_err_frac_{row['Signal region']}"
        error = round(all_error[error_key] * 100, 2)
        
        return f"({percentage} " + r"$\pm$ " + f"{error})" r"\%"
    else:
        return row[opts.Ana]
    

def join_values(group):
    if group["Operator"].values[0] == "Any":
        return {'Signal Region': 'Total', opts.Ana: '\n' + '\n'.join(map(str, group[opts.Ana])) + '\n'}
    else:   
        return {'Signal Region': 'merged\nresolved', opts.Ana: '\n' + '\n'.join(map(str, group[opts.Ana])) + '\n'}


data = []
for key, value in all_counts.items():
  operator, region = key.split('_')
  data.append({'Operator': operator, 'Signal region': region, opts.Ana: value})

# Create DataFrame
df = pd.DataFrame(data)

# Check if all 'total' values are the same for each operator
total_values = df[df['Signal region'] == 'total'][opts.Ana]
if total_values.nunique() > 1:
  raise ValueError("Not all 'total' values are the same for each operator")

# Get the common 'total' value
common_total_value = total_values.iloc[0]

# Add a new row at the beginning with 'Operator' = 'Any', 'Signal region' = 'total', and 'value' as the common value
df = pd.concat([pd.DataFrame([{'Operator': 'Any', 'Signal region': 'total', opts.Ana: r"100\%"+f" ({common_total_value})"}]), df], ignore_index=True)

# Filter out rows where 'Signal region' is 'total' and 'Operator' is not 'Any'
df = df[~((df['Signal region'] == 'total') & (df['Operator'] != 'Any'))]


# Apply the format_value function to the 'value' column
df[opts.Ana] = df.apply(format_value, axis=1)
df = df.sort_values(by='Operator')

# Group by 'Operator' and apply the join_values function to the 'value' column

df_grouped = df.groupby('Operator').apply(join_values).apply(pd.Series).reset_index()

# Rename the 0 column to 'value'
df_grouped.rename(columns={0: opts.Ana}, inplace=True)

# Sort by Operator with SM operator being on top
df_grouped['sort'] = df_grouped['Operator'].apply(lambda x: 0 if x=='Any' else (1 if x=='SM' else 2))
df_grouped.sort_values(by=['sort', 'Operator'], inplace=True)
df_grouped.drop(columns='sort', inplace=True)

# Create a new figure with a single subplot
fig, ax = plt.subplots(1, 1)
ax.axis('off')

# Remove the index
df_grouped.index = [''] * len(df_grouped)

# Create a table from the DataFrame and add it to the subplot
tab = table(ax, df_grouped, loc='center', cellLoc='center', rowLabels=[])

# Adjust cell padding to fit sublines
cellDict = tab.get_celld()
for i in range(0,len(df_grouped.columns)):
    for j in range(0,len(df_grouped)+1):
        cellDict[(j,i)].set_height(0.1)
        cellDict[(j,i)].PAD = 0.05
        

# Save the figure as a PDF
plt.savefig('./Tables/table_Works04_.pdf', format='pdf')









    