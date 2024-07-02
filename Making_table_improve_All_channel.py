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
parser.add_option("--nb_lep", default = 2)
parser.add_option("--Name", default= "")
parser.add_option("--All_channel", default = False)
parser.add_option("--All_ops", default= True)

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
Details_dir= "Tables/All_channel/Stats/"

# Check if the directory exists, if not, create it
if not os.path.exists(output_table):
    os.makedirs(output_table)
    
processes_all = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
decays_all = ["llqq",'lvqq','vvqq']

if opts.All_ops:
    all_ops_SM =  ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
else:
    all_ops_SM = ["SM","FM0","FM2","FS1","FT1","FT5","FT8"]

processes_1lep = ["WmZ","WpZ","WmWm","WpWm","WpWp"]
decays_1lep = ['lvqq'] 

processes_2lep=["WpZ","WmZ","ZZ"]
decays_2lep = ['llqq']   

processes_0lep=["WpZ","WmZ","ZZ"]
decays_0lep = ['vvqq']

Lepton_decay = {
    "0Lepton": [f"{process}_vvqq" for process in processes_0lep],  
    "1Lepton": [f"{process}_lvqq" for process in processes_1lep],
    "2Lepton": [f"{process}_llqq" for process in processes_2lep]
}

lepton_channel = {
    "0Lepton": [processes_0lep,decays_0lep],
    "1Lepton": [processes_1lep,decays_1lep],
    "2Lepton": [processes_2lep,decays_2lep]
}

nb_lepton= int(opts.nb_lep)
if opts.All_channel:
    processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
    decays = ['lvqq','llqq','vvqq']
else:
    if nb_lepton == 1:
        print("1 lepton")
        processes = ["WmZ","WpZ","WmWm","WpWm","WpWp"]
        decays = ['lvqq']
    elif nb_lepton == 2:
        processes=["WpZ","WmZ","ZZ"]
        decays = ['llqq']
    elif nb_lepton == 0:
        processes=["WpZ","WmZ","ZZ"]
        decays = ['vvqq']
    
lumi = 139
    
Cutflow_paths = {}

for decay in decays_all:
    for process in processes_all:
        if uf.possible_process(process,decay):
            for op in all_ops_SM:
                if op=="SM":
                    path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/{Details_dir}/combined/")
                else:
                    path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/{Details_dir}/")
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

all_ops_SM = ["SM","FM0","FM2","FS1","FT1","FT5","FT8"]

Cutflow_paths_dict = {}
all_counts_dict_channel = {}
all_error_dict_channel = {}
all_counts_dict_norm_channel={}
Sum_process_channel={}
Sum_error_channel={}





for channel, processes_decays in lepton_channel.items():
    processes = processes_decays[0]
    decays = processes_decays[1]
    all_counts_dict = {}
    all_error_dict = {}
    all_counts_dict_norm={}
    for decay in decays:
        for process in processes:
            phys_process = f"{process}_{decay}"
            Cutflow_path = {}
            for op in all_ops_SM:
                if f"{process}_{decay}_{op}" in Cutflow_paths and Cutflow_paths[f"{process}_{decay}_{op}"]:
                    #print(f"Processing operator {op}")
                    if op == "SM":
                        cutflow_merged_path = Cutflow_paths[f"{process}_{decay}_{op}"] + "cutflow_merged_total.txt"
                        cutflow_resolved_path = Cutflow_paths[f"{process}_{decay}_{op}"] + "cutflow_resolved_total.txt"
                        cutflow_frac_res_path = Cutflow_paths[f"{process}_{decay}_{op}"] +"frac_err_resolved.txt"
                        cutflow_frac_merged_path = Cutflow_paths[f"{process}_{decay}_{op}"] +"frac_err_merged.txt"
                    else:
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
                    if op == "SM":
                        Total_nb_events_merged = float(Total_line[1])
                    else:
                        Total_nb_events_merged = float(Total_line[0])
                    Merged_line = Lines[-2].split()
                    Merged_nb_events = float(Merged_line[2])

                with open(cutflow_resolved_path, 'r') as f:
                    Lines = [line for line in f.readlines()]
                    Total_line = Lines[2].split()
                    if op == "SM":
                        Total_nb_events_resolved = float(Total_line[1])
                    else:
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



    for decay in decays:
        for process in processes:
            phys_process = f"{process}_{decay}"
            dict_norm={}
            for op in all_ops_SM:
                if op=="SM":
                    xsection_fb= uf.take_xsec_fb2("FM0", "SM", process, decay)
                    if xsection_fb is None:
                        #xsection_fb = uf.cross_section_fb("FM0", "SM", process, decay)
                        if xsection_fb is None:
                            continue
                else:
                    xsection_fb= uf.take_xsec_fb2(op, "QUAD", process, decay)
                    if xsection_fb is None:
                        #xsection_fb = uf.cross_section_fb(op, "QUAD", process, decay)
                        if xsection_fb is None:
                            continue
                lumi_scale = lumi * xsection_fb
                total_events = all_counts_dict.get(phys_process, {}).get(f"{op}_total", 0)
                merged_events = all_counts_dict.get(phys_process, {}).get(f"{op}_merged", 0)
                resolved_events = all_counts_dict.get(phys_process, {}).get(f"{op}_resolved", 0)
                if total_events > 0:  # Check to avoid division by zero
                    dict_norm[op + "_total"] = lumi_scale * total_events / total_events
                    dict_norm[op + "_merged"] = lumi_scale * merged_events / total_events
                    dict_norm[op + "_resolved"] = lumi_scale * resolved_events / total_events
                else:
                    continue   
            all_counts_dict_norm[phys_process] = dict_norm
    all_counts_dict_norm_channel[channel] = all_counts_dict_norm
    all_error_dict_channel[channel] = all_error_dict

    Sum_process={}
    for phys_process, dict_norm in all_counts_dict_norm.items():
        if phys_process in Lepton_decay[channel]:
            #print(phys_process)
            for op_key, value in dict_norm.items():
                # Skip adding WmZ and WpZ contributions to SM to avoid double counting
                if 'WmZ' in phys_process and 'SM' in op_key:
                    print(f"Skipping {phys_process} {op_key}")
                    continue
                if 'WmWm' in phys_process and 'SM' in op_key:
                    print(f"Skipping {phys_process} {op_key}")
                    continue
                if 'WpWp' in phys_process and 'SM' in op_key:
                    print(f"Skipping {phys_process} {op_key}")
                    continue
                
                if op_key in Sum_process:
                    Sum_process[op_key] += value
                else:
                    Sum_process[op_key] = value
    Sum_process_channel[channel] = Sum_process
    
            
    Sum_error = {}
    for phys_process, dict_norm in all_error_dict.items():
        if phys_process in Lepton_decay[channel]:
            for op_key, value in dict_norm.items():
                # Skip adding WmZ and WpZ contributions to SM to avoid double counting
                if 'WmZ' in phys_process and 'SM' in op_key:
                    print(f"Skipping {phys_process} {op_key}")
                    continue
                if 'WmWm' in phys_process and 'SM' in op_key:
                    print(f"Skipping {phys_process} {op_key}")
                    continue
                if 'WpWp' in phys_process and 'SM' in op_key:
                    print(f"Skipping {phys_process} {op_key}")
                    continue
                if op_key in Sum_error:
                    # Sum the inverse of the square of the errors
                    Sum_error[op_key] += 1 / (value ** 2)
                else:
                    Sum_error[op_key] = 1 / (value ** 2)
    for op_key in Sum_error:
        Sum_error[op_key] = (Sum_error[op_key] ** -0.5)
    Sum_error_channel[channel] = Sum_error
                

# Convert the sum back by taking the square root of the inverse

    
    all_counts_dict_norm[channel] = Sum_process
    #print(channel)
    #print(all_counts_dict_norm[channel])
    all_error_dict[channel] = Sum_error
    all_error_dict_channel[channel][channel] = Sum_error_channel[channel]

    all_counts_dict_norm_All= all_counts_dict_norm
    for phys_process, dict_norm in all_counts_dict_norm.items():
        updates = {}  # Collect updates here
        for op in all_ops_SM:
            # Check if the operator-related keys are missing and prepare updates
            if not any(op in key for key in dict_norm):
                updates[op + "_total"] = 0
                updates[op + "_merged"] = 0
                updates[op + "_resolved"] = 0
        
        # Apply updates outside the inner loop
        all_counts_dict_norm_All[phys_process].update(updates)
        
    all_counts_dict_norm_channel[channel] = all_counts_dict_norm_All


#print(all_counts_dict_norm_channel["0Lepton"])
xsection_dict = {}
def format_and_join_values(group, phys_process,channel):

    #print(f"Processing {process}_{decay}")
    formatted_values = {'total': '', 'merged': '', 'resolved': ''}
    for _, row in group.iterrows():
        if phys_process in row:
            total_value = df[(df['Signal region'] == 'total') & (df['Operator'] == row['Operator'])][phys_process].values
            if total_value.size > 0:
                total_value = total_value[0]
                yield_ = row[phys_process]
                if pd.notnull(yield_) and total_value != 0:
                    if row['Signal region'] == 'total':
                        yield_normalised = yield_
                        percentage = round(yield_/total_value * 100, 2)
                    #print("Phys process %d operator %d, Yield %d",phys_process, row['Operator'],yield_)
                    else:
                        yield_normalised = yield_/total_value
                        percentage = round(yield_normalised * 100, 2)
                    if ((yield_normalised) > 0.09) :
                        yield_scaled = round(yield_normalised, 2)
                    else:
                        yield_scaled = round(yield_normalised, 4)
                    
            
                    if row['Signal region'] == 'total':
                        formatted_values['total'] = f"{yield_scaled}  ({percentage})" r"\%"
                    
                    elif row['Signal region'] in ['merged', 'resolved']:
                        all_error_=all_error_dict_channel[channel][phys_process]
                        error_key = f"{row['Operator']}_err_frac_{row['Signal region']}"
                        error = round(all_error_[error_key] * 100, 2)
                        
                        formatted_values[row['Signal region']] = f"({percentage} " + r"$\pm$ " + f"{error})" r"\%"
                
            else:
                formatted_values[row['Signal region']] = row[phys_process]
    
    return pd.Series({'Signal Region': 'total\nmerged\nresolved', phys_process: '\n'.join(formatted_values.values())})


for channel, all_counts_dict_norm_All in all_counts_dict_norm_channel.items():
    data = []
    phys_processes = [f"{process}_{decay}" for process in processes for decay in decays]
    df_grouped_list = []            
    for phys_process, all_counts in all_counts_dict_norm_All.items():
        #print(phys_process)
        data=[]
        for key, value in all_counts.items():
            operator, region = key.split('_')
            data.append({'Operator': operator, 'Signal region': region, phys_process: value})

            # Create DataFrame
        df = pd.DataFrame(data)

        df = df.sort_values(by='Operator')

        # Apply the format_and_join_values function to each physics process column       
        df_grouped = df.groupby('Operator').apply(lambda group: format_and_join_values(group, phys_process,channel)).reset_index()
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
            cellDict[(j,i)].set_height(0.06)
            cellDict[(j,i)].PAD = 0.0255
            
    fig.set_size_inches(12, 12)

    # Hide the axes
    ax.axis('off')


    output_dir_='./Tables/Test_less_op/'
    os.makedirs(output_dir_, exist_ok=True)
    # Save the figure as a PDF
    plt.savefig(output_dir_+f'table_{channel}_{opts.Name}.pdf', format='pdf')


df_grouped_list = []  
for channel, all_counts_dict_norm_All in all_counts_dict_norm_channel.items():
    phys_processes = [f"{process}_{decay}" for process in processes for decay in decays]
              
    for phys_process, all_counts in all_counts_dict_norm_All.items():
        if phys_process == channel:
            data=[]
            for key, value in all_counts.items():
                operator, region = key.split('_')
                data.append({'Operator': operator, 'Signal region': region, phys_process: value})

                # Create DataFrame
            df = pd.DataFrame(data)

            df = df.sort_values(by='Operator')

            # Apply the format_and_join_values function to each physics process column       
            df_grouped = df.groupby('Operator').apply(lambda group: format_and_join_values(group, phys_process,channel)).reset_index()
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
        cellDict[(j,i)].set_height(0.06)
        cellDict[(j,i)].PAD = 0.02
        
fig.set_size_inches(12, 12)

# Hide the axes
ax.axis('off')



# Save the figure as a PDF
plt.savefig(output_dir_+f'table_All_channel_sum_{opts.Name}.pdf', format='pdf')







    