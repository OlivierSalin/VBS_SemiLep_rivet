import glob
import os
import subprocess
from array import array
import ROOT
import itertools
import time
import shutil
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import re

processes = ["WpZ","ZZ","WmZ"]
decays = ["llqq"]

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]

detail_path=f"/Tables/VBS_jet_algo_plot_bis/Highest_mass_/2Lepton/"
detail_path=f"/Tables/VBS_jet_algo_eta_mCut200_plot/Highest_eta_mCut200/2Lepton/"
#output_dir="/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/BDT_sample/Polar/"
output_dir=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path}/Combined/"





def merge_cutflows(paths, out_file,process,decay):
    # Initialize a dictionary to hold the total counts
    
    Cutflow_type = ["merged", "resolved"]
    # Iterate over each path
    for cutflow_type in Cutflow_type:
        total_counts = {}
        for path_ in paths:
            path = path_ + f"cutflow_{cutflow_type}.txt"
            # Open the cutflow file
            with open(path, 'r') as f:
                # Skip the first line (header)
                next(f)
                next(f)
                # Read the remaining lines
                first_line = True
                # Read the remaining lines
                for line in f:
                    # Split the line into fields
                    fields = line.split()
                    if not first_line:
                        
                        if len(fields) < 2 or not fields[1]:
                            break
                        # Get the cut name
                        cut_name = fields[1]
                        #print(cut_name)
                        if fields[2]=="veto":
                            count = float(fields[3])
                        else:
                            count = float(fields[2])
                        #print(count)

                        # Add the count to the total for this cut
                        if cut_name in total_counts:
                            total_counts[cut_name] += count
                        else:
                            total_counts[cut_name] = count
                    else:
                        # If this is the first line, just add the count to the total
                        count = float(fields[0])
                        if "Total" in total_counts:
                            total_counts["Total"] += count
                        else:
                            total_counts["Total"]=count


                    # Update the flag for the first line
                    first_line = False

        # Calculate the total count across all cuts
        total_count = total_counts["Total"]

        # Initialize the list of lines for the output file
        lines = [f"{process}_{decay}_selections cut-flow:\n"]
        if cutflow_type == "merged":
            lines.append("                            Count       A_cumu         A_incr\n")
        else:
            lines.append("                                      Count    A_cumu    A_incr\n")

        # Initialize the previous count for calculating incremental percentages
        prev_count = total_count

        # For each cut name and total count
        for cut_name, count in total_counts.items():
            # Skip the calculation if prev_count is zero to avoid division by zero
            if prev_count == 0:
                continue

            # Calculate the cumulative and incremental percentages
            cumu_percent = 100 * count / total_count
            incr_percent = 100 * count / prev_count
            
            cumul_percent = f"{cumu_percent:.1f}%"
            incre_percent = f"{incr_percent:.1f}%"

            # Add a line to the output list
            if cut_name == "Total":
                cut_name = ""
                incre_percent="-"
            if cutflow_type == "merged":
                lines.append(f"Pass {cut_name:<15}\t{count:>9}\t{cumul_percent:>8}\t{incre_percent:>6}\n")
            else:
                lines.append(f"Pass {cut_name:<25}\t{count:>9}\t{cumul_percent:>8}\t{incre_percent:>6}\n")

            # Update the previous count
            prev_count = count
        out_file_ = out_file + f"{process}_{decay}_cutflow_{cutflow_type}_aQGC_total.txt"
        # Write the lines to the output file
        with open(out_file_, 'w') as f:
            f.writelines(lines)

def path_ntuple_jetAlgo(dir_details,Processes,Decays,Operators):
    Path= {}
    base_dir= f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{dir_details}/"
    for process in Processes:
        for decay in Decays:
            for op in Operators:
                proc_dec = f"{process}_{decay}"
                run_dir= base_dir + f"{proc_dec}/"
                if op=="SM":
                    path_ntuple =run_dir + f"/FM0_SM/"
                else:
                    path_ntuple =run_dir + f"/{op}_QUAD/"
                     
                if os.path.exists(path_ntuple):
                    print(f"For process {proc_dec}, the path exists.: {path_ntuple}")
                    Path[f"{process}_{decay}_{op}"] = path_ntuple
                else:
                    print(f"For process {proc_dec}_{op}, the path does not exist.")
    return Path
            
def combine_ntuples_jetAlgo(dir_details, Processes, Decays, Operators,output_dir,truth=False,aQGC=False):
    # Get the paths of the ntuple files
    paths = path_ntuple_jetAlgo(dir_details, Processes, Decays, Operators)
    ntuple_files_aQGC = []
    ntuple_files_SM = []

    # Define the directory where the files will be combined
    combined_dir = output_dir
    individual_dir = output_dir + "individual/"

    for process in Processes:
        for decay in Decays:
            
            proc_decay= f"{process}_{decay}"
            combined_dir_bis = combined_dir + f"{proc_decay}/"
            # Create the combined directory if it doesn't exist
            if not os.path.exists(combined_dir):
                os.makedirs(combined_dir)

            if not os.path.exists(combined_dir_bis):
                os.makedirs(combined_dir_bis)
                
            # Combine the ntuple files using hadd
            if aQGC:
                
                ntuple_files = [value + "/ntuple_rivet.root" for key, value in paths.items() if proc_decay in key and "SM" not in key]
                subprocess.run(["hadd", "-f", combined_dir_bis + f"/aQGC_ntuple_rivet.root"] + ntuple_files)
                ntuple_files_aQGC.extend(ntuple_files)
                
                path_cutflow = [value  for key, value in paths.items() if proc_decay in key and "SM" not in key]
                merge_cutflows(path_cutflow, combined_dir_bis,process,decay)
                
                ntuple_SM = [value + "/ntuple_rivet.root" for key, value in paths.items() if proc_decay in key and "SM" in key]
                ntuple_files_SM.extend(ntuple_SM)
                subprocess.run(["hadd", "-f", combined_dir_bis + f"/SM_ntuple_rivet.root"] + ntuple_SM)
                
                ntuple_files = [value + "/ntuple_rivet.root" for key, value in paths.items() if proc_decay in key]
                subprocess.run(["hadd", "-f", combined_dir_bis + f"/ntuple_rivet.root"] + ntuple_files)                
            else:
                ntuple_files = [value + "/ntuple_rivet.root" for key, value in paths.items() if proc_decay in key]
                subprocess.run(["hadd", "-f", combined_dir_bis + f"/ntuple_rivet.root"] + ntuple_files)
                


            if not os.path.exists(individual_dir):
                os.makedirs(individual_dir)
                
            for op in Operators:
                for key, value in paths.items():
                    if f"{proc_decay}_{op}" in key:
                        src_file = value + "/ntuple_rivet.root"
                        dest_file = individual_dir + f"/{proc_decay}_{op}_ntuple_rivet.root"
                        shutil.copy(src_file, dest_file)


            print(f"All ntuple files have been combined in the directory: {combined_dir}")
    
    if ntuple_files_aQGC:
        subprocess.run(["hadd", "-f", combined_dir + "/All_aQGC_ntuple_rivet.root"] + ntuple_files_aQGC)

    # Combine all SM ntuple files
    if ntuple_files_SM:
        subprocess.run(["hadd", "-f", combined_dir + "/All_SM_ntuple_rivet.root"] + ntuple_files_SM)

    ntuple_files_all = [values + "/ntuple_rivet.root" for values in paths.values()]
    subprocess.run(["hadd", "-f", combined_dir + f"/All_ntuple_rivet.root"] + ntuple_files_all)
    
    with open(combined_dir + f"/All_paths.txt", "w") as file:
        for path in ntuple_files_all:
            file.write(path + "\n")
            
            
#combine_ntuples_jetAlgo(detail_path,processes, decays, all_ops_cat,output_dir,truth=False,aQGC=True)