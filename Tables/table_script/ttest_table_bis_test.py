import glob
import os
from array import array

import shutil
import matplotlib.pyplot as plt
import numpy as np
import math

plt.rcParams['text.usetex'] = True
import pandas as pd   
from pandas.plotting import table
import utils_func as uf

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Ana", default = "WpZ_llqq")
parser.add_option("--DOCUT", default = "YES")

opts, _ = parser.parse_args()
base_dir_bis = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables"

output_table = f"{base_dir_bis}/{opts.Ana}/Plots/Table/test01/DOCUT_{opts.DOCUT}/"

all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

processes=["WmZ","WpZ","ZZ"]
decays=["llqq"]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"


# Check if the directory exists, if not, create it
if not os.path.exists(output_table):
    os.makedirs(output_table)
    
lumi = 139
    
Cutflow_paths = {
    "SM": f"{base_dir_bis}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_FM0_SM_EXT0/DOCUT_{opts.DOCUT}/",
}
for process in ["WpZ","WmZ","ZZ"]:
    for decay in ["llqq"]:
        for op in all_ops:
            path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables_file/")
            matches = glob.glob(path)
            #print(matches)
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            Cutflow_paths[op] = matches[0]
print(Cutflow_paths["SM"])

Cutflow_path = {}

for process in ["WpZ","WmZ","ZZ"]:
    for decay in ["llqq"]:
        for op in all_ops_SM:
            if op in Cutflow_paths and Cutflow_paths[op]:
                #print(f"Processing operator {op}")
                base_path = Cutflow_paths[op]
                paths = {
                    "": base_path,
                    "_merged": base_path + "cutflow_merged.txt",
                    "_resolved": base_path + "cutflow_resolved.txt",
                    "_frac_res": base_path + "frac_after_cuts_error_bar_resolved.txt",
                    "_frac_merged": base_path + "frac_after_cuts_error_bar_merged.txt"
                }
                
                if all(os.path.exists(path) for path in paths.values()):
                    Cutflow_path.update({f"{process}_{decay}_{op}{suffix}": path for suffix, path in paths.items()})
                    
#print(Cutflow_path)
def extract_counts(Cutflow_path, op):
    Info_table = {}
    Info_err_frac={}
    cutflow_merged_path = Cutflow_path[f"{process}_{decay}_{op}_merged"]
    cutflow_resolved_path = Cutflow_path[f"{process}_{decay}_{op}_resolved"]
    cutflow_frac_res_path = Cutflow_path[f"{process}_{decay}_{op}_frac_res"]
    cutflow_frac_merged_path = Cutflow_path[f"{process}_{decay}_{op}_frac_merged"]
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
        print("Key: ",op)
        if op.endswith("_merged") or op.endswith("_resolved") or op.endswith("_frac_res") or op.endswith("_frac_merged"):
            continue
        Info_table, Info_err_frac  = extract_counts(Cutflow_path, op)
        all_counts.update(Info_table)
        all_error.update(Info_err_frac)
    return all_counts, all_error

# Usage
all_counts,all_error = get_all_counts(Cutflow_path)
print(all_counts)
print(all_error)

process=opts.Ana.split("_")[0]
decay=opts.Ana.split("_")[1]
xsection_SM=uf.cross_section_fb(EFT_op='FM0', EFT_type='SM', proces=process, dec=decay)
lumi_scalling_SM= round(lumi*xsection_SM, 2)
# Convert the dictionary to a DataFrame


""" 
xsection_dict = {}
def format_and_join_values(group):
    formatted_values = {'total': '', 'merged': '', 'resolved': ''}
    for _, row in group.iterrows():
        if row['Signal region'] in ['total','merged', 'resolved']:
            if row['Operator'] == 'SM':
                lumiscaling = lumi_scalling_SM
            else:
                key=f"{process}_{decay}_{row['Operator']}" 
                if key in xsection_dict:
                    xsection = xsection_dict[key]
                else:
                    xsection = uf.cross_section_fb(EFT_op=row['Operator'], EFT_type='QUAD',  proces=process, dec=decay)
                    xsection_dict[key] = xsection  # Store the calculated cross section in the dictionary
                    #print(f"\nCross section of Operator: {row['Operator']}, EFT_type: QUAD, process: {process}, decay: {decay} is {xsection} fb")
                    
                lumiscaling = lumi*xsection
                #print(f"\nCross section of Operator: {row['Operator']}, EFT_type: QUAD, process: {process}, decay: {decay} is {xsection} fb")
                #print(f"Luminosity scaling factor: {lumiscaling}")
                
            
            total_value = df[(df['Signal region'] == 'total') & (df['Operator'] == row['Operator'])][opts.Ana].values[0]
            yield_= row[opts.Ana]
            yield_normalised = row[opts.Ana]/total_value
            if ((yield_normalised* lumiscaling) > 0.09) :
                yield_scaled = round(yield_normalised * lumiscaling, 2)
            else:
                yield_scaled = round(yield_normalised * lumiscaling, 4)
            percentage = round(yield_normalised * 100, 2)
        
            if row['Signal region'] == 'total':
                formatted_values['total'] = f"{yield_scaled}  ({percentage})" r"\%"
            
            elif row['Signal region'] in ['merged', 'resolved']:
                error_key = f"{row['Operator']}_err_frac_{row['Signal region']}"
                error = round(all_error[error_key] * 100, 2)
                
                formatted_values[row['Signal region']] = f"{yield_scaled}  ({percentage} " + r"$\pm$ " + f"{error})" r"\%"
            
        else:
            formatted_values[row['Signal region']] = row[opts.Ana]
    
    return {'Signal Region': 'total\nmerged\nresolved', opts.Ana: '\n'.join(formatted_values.values())}

data = []
for key, value in all_counts.items():
  operator, region = key.split('_')
  data.append({'Operator': operator, 'Signal region': region, opts.Ana: value})
  #print(data)

# Create DataFrame
df = pd.DataFrame(data)

# Apply the format_value function to the 'value' column
#df[opts.Ana] = df.apply(format_value, axis=1)
df = df.sort_values(by='Operator')

# Group by 'Operator' and apply the join_values function to the 'value' column
#df_grouped = df.groupby('Operator').apply(join_values).apply(pd.Series).reset_index()
df_grouped = df.groupby('Operator').apply(format_and_join_values).apply(pd.Series).reset_index()
# Rename the 0 column to 'value'
df_grouped.rename(columns={0: opts.Ana}, inplace=True)

# Sort by Operator with SM operator being on top
df_grouped['sort'] = df_grouped['Operator'].apply(lambda x: 0 if x=='SM' else 1)
df_grouped.sort_values(by=['sort', 'Operator'], inplace=True)
df_grouped.drop(columns='sort', inplace=True)

fig, ax = plt.subplots(1, 1)



# Remove the index
df_grouped.index = [''] * len(df_grouped)

# Create a table from the DataFrame and add it to the subplot
 
# Adjust cell padding to fit sublines
cellDict = tab.get_celld()
for i in range(0,len(df_grouped.columns)):
    for j in range(0,len(df_grouped)+1):
        cellDict[(j,i)].set_height(0.06)
        cellDict[(j,i)].PAD = 0.03
        
fig.set_size_inches(12, 12)

# Hide the axes
ax.axis('off')

# Save the figure as a PDF
plt.savefig('./Tables/table_Lumi04_.pdf', format='pdf')

# Convert the DataFrame to a LaTeX table
latex_table = df_grouped.style.to_latex()

# Write the LaTeX table to a file
with open('./Tables/table_Lumi04_.tex', 'w') as f:
    f.write(latex_table)









     """