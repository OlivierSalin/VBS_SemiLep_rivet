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
    
Cutflow_paths = {}

for decay in ["llqq"]:
    for process in ["WpZ","WmZ","ZZ"]:
        for op in all_ops_SM:
            if op=="SM":
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_first/")
            else:
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_first/")
            matches = glob.glob(path)
            #print(matches)
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            Cutflow_paths[f"{process}_{decay}_{op}"] = matches[0]
#print(Cutflow_paths)
all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

Cutflow_paths_dict = {}
all_counts_dict = {}
all_error_dict = {}

for decay in ["llqq"]:
    for process in ["WpZ","WmZ","ZZ"]:
        phys_process = f"{process}_{decay}"
        Cutflow_path = {}
        for op in all_ops_SM:
            if f"{process}_{decay}_{op}" in Cutflow_paths and Cutflow_paths[f"{process}_{decay}_{op}"]:
                #print(f"Processing operator {op}")
                cutflow_merged_path = Cutflow_paths[f"{process}_{decay}_{op}"] + "cutflow_merged.txt"
                cutflow_resolved_path = Cutflow_paths[f"{process}_{decay}_{op}"] + "cutflow_resolved.txt"
                cutflow_frac_res_path = Cutflow_paths[f"{process}_{decay}_{op}"] +"frac_after_cuts_error_bar_resolved.txt"
                cutflow_frac_merged_path = Cutflow_paths[f"{process}_{decay}_{op}"] +"frac_after_cuts_error_bar_merged.txt"

                if os.path.exists(cutflow_merged_path) and os.path.exists(cutflow_resolved_path):
                    if os.path.exists(cutflow_frac_res_path) and os.path.exists(cutflow_frac_merged_path):
                        Cutflow_path[op]= Cutflow_paths[f"{process}_{decay}_{op}"]
                        Cutflow_path[op + "_merged"]= cutflow_merged_path
                        Cutflow_path[op + "_resolved"]= cutflow_resolved_path
                        Cutflow_path[op + "_frac_merged"]= cutflow_frac_merged_path 
                        Cutflow_path[op + "_frac_res"]= cutflow_frac_res_path

        #print(Cutflow_path["SM"])
        #print(Cutflow_path)
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
        Cutflow_paths_dict[phys_process] = Cutflow_path
        all_counts_dict[phys_process] = all_counts
        all_error_dict[phys_process] = all_error
#print(all_counts_dict)
#print(all_error_dict)




xsection_dict = {}
def format_and_join_values(group, phys_process):
    process=phys_process.split("_")[0]
    decay=phys_process.split("_")[1]
    #print(f"Processing {process}_{decay}")
    formatted_values = {'total': '', 'merged': '', 'resolved': ''}
    for _, row in group.iterrows():
        if phys_process in row:
            if row['Signal region'] in ['total','merged', 'resolved']:
                if row['Operator'] == 'SM':
                    xsection_SM=uf.cross_section_fb(EFT_op='FM0', EFT_type='SM', proces=process, dec=decay)
                    lumi_scalling_SM= round(lumi*xsection_SM, 2)
                    lumiscaling = lumi_scalling_SM
                else:
                    key=f"{process}_{decay}_{row['Operator']}" 
                    if key in xsection_dict:
                        xsection = xsection_dict[key]
                    else:
                        xsection = uf.cross_section_fb(EFT_op=row['Operator'], EFT_type='QUAD',  proces=process, dec=decay)
                        xsection_dict[key] = xsection  # Store the calculated cross section in the dictionary
                    
                    lumiscaling = lumi*xsection
                
            total_value = df[(df['Signal region'] == 'total') & (df['Operator'] == row['Operator'])][phys_process].values
            if total_value.size > 0:
                total_value = total_value[0]
                #print(f"Phys process {phys_process}, total value: {total_value}")
                yield_ = row[phys_process]
                #print("Phys process %d operator %d, Yield %d",phys_process, row['Operator'],yield_)
                if pd.notnull(yield_) and total_value != 0:
                    yield_normalised = yield_/total_value
                    if ((yield_normalised* lumiscaling) > 0.09) :
                        yield_scaled = round(yield_normalised * lumiscaling, 2)
                    else:
                        yield_scaled = round(yield_normalised * lumiscaling, 4)
                    percentage = round(yield_normalised * 100, 2)
            
                    if row['Signal region'] == 'total':
                        formatted_values['total'] = f"{yield_scaled}  ({percentage})" r"\%"
                    
                    elif row['Signal region'] in ['merged', 'resolved']:
                        all_error_=all_error_dict[phys_process]
                        error_key = f"{row['Operator']}_err_frac_{row['Signal region']}"
                        error = round(all_error_[error_key] * 100, 2)
                        
                        formatted_values[row['Signal region']] = f"({percentage} " + r"$\pm$ " + f"{error})" r"\%"
                
            else:
                formatted_values[row['Signal region']] = row[phys_process]
    
    return pd.Series({'Signal Region': 'total\nmerged\nresolved', phys_process: '\n'.join(formatted_values.values())})

data = []
phys_processes = ["WmZ_llqq","WpZ_llqq","ZZ_llqq"]
df_grouped_list = []            
for phys_process, all_counts in all_counts_dict.items():
    data=[]
    for key, value in all_counts.items():
        operator, region = key.split('_')
        data.append({'Operator': operator, 'Signal region': region, phys_process: value})

        # Create DataFrame
    df = pd.DataFrame(data)

    df = df.sort_values(by='Operator')

    # Apply the format_and_join_values function to each physics process column       
    df_grouped = df.groupby('Operator').apply(lambda group: format_and_join_values(group, phys_process)).reset_index()
    df_grouped_list.append(df_grouped)

# Merge the DataFrames
df_final = pd.concat(df_grouped_list, axis=1)
df_final = df_final.loc[:,~df_final.columns.duplicated()]

df_final = df_final.sort_values(by='Operator')
df_final['sort'] = df_final['Operator'] == 'SM'
df_final = df_final.sort_values(by=['sort', 'Operator'], ascending=[False, True]).drop('sort', axis=1)

# Create a figure and a subplot
fig, ax = plt.subplots(1, 1)

# Remove the index
df_final.index = [''] * len(df_final)

# Create a table from the DataFrame and add it to the subplot
tab = table(ax, df_final, loc='center', cellLoc='center', rowLabels=[])

# Adjust cell padding to fit sublines
cellDict = tab.get_celld()
for i in range(0,len(df_final.columns)):
    for j in range(0,len(df_final)+1):
        cellDict[(j,i)].set_height(0.07)
        cellDict[(j,i)].PAD = 0.02
        
fig.set_size_inches(12, 12)

# Hide the axes
ax.axis('off')

# Save the figure as a PDF
plt.savefig('./Tables/table_Proces04_.pdf', format='pdf')









    