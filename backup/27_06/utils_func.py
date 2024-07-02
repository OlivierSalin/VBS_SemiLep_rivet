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
import pyAMI.client
import pyAMI.atlas.api as AtlasAPI
client = pyAMI.client.Client('atlas')
AtlasAPI.init()

def extract_EFT_op_proces_dec(conf):
    """
    This function takes conf as a parameter,
    and returns the EFT operator, process, and decay.
    """
    # Split the conf name with the underscore 
    parts = conf.split('_')
    EFT_op,EFT_type, proc, decay = parts[3], parts[4], parts[1], parts[2]
    
    return EFT_op, EFT_type, proc, decay

def extract_EFT_op_proces_dec_bis(conf):
    """
    This function takes conf as a parameter,
    and returns the EFT operator, process, and decay.
    """
    pattern = r"MGPy8EG_aQGC(.*)_(.*)_1_(.*)_(.*)"
    match = re.search(pattern, conf)
    if match:
        EFT_op = match.group(1)
        EFT_type = match.group(2)
        proc = match.group(3)
        decay = match.group(4)
        
    elif conf.startswith("user."):
        # If the conf matches the Sample_name_ pattern
        parts = conf.split('_')
        EFT_op,EFT_type, proc, decay = parts[3], parts[4], parts[1], parts[2]
    else:
        raise ValueError("Invalid conf format")

    return EFT_op, EFT_type, proc, decay


def cross_section_fb(EFT_op, EFT_type,  proces, dec):
    """
    This function takes an EFT operator, a process, and a decay as parameters,
    and returns the cross section for the corresponding dataset.
    """
    
    all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
               "FS02","FS1",
               "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    all_ops2 = ["FM","FS","FT"]
    all_ops_ = [op + "_QUAD" for op in all_ops]
    Processes = [proces]
    Decay = [dec]

    organized_info = {}

    if EFT_type == 'QUAD':
        for process in Processes:
            for decay in Decay:  
                pattern = f"%%MGPy8EG_aQGCF%%_{process}_{decay}%%"
                List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                fields=['ldn', 'cross_section', 'dataset_number'], 
                                                limit=[1, 100], type='EVNT', project='mc16%')

                if List_ami:
                    for dataset in List_ami:
                        operator = next((op for op in all_ops_ if op in dataset['ldn']), None)
                        if operator:
                            key = f"{operator}_{process}_{decay}"
                            if key not in organized_info:
                                organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                            organized_info[key]['name'].append(dataset['ldn'])
                            organized_info[key]['cross_section'].append(dataset['cross_section']) # Cross section in nb
                            organized_info[key]['dataset_number'].append(dataset['dataset_number'])
                else:
                    print(f'No datasets found for {process}_{decay}\n')
                        
    elif (EFT_type == 'SM'):
        for process in Processes:
            for decay in Decay:
                if process == 'WmZ' or process == 'WpZ':
                    process_ = 'WZ'
                if process == "ZZ":
                    process_= "ZZ"
                pattern = f'%{process_}jj_llqq_EW6%'
                List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                fields=['ldn', 'cross_section', 'dataset_number'], 
                                                limit=[1, 100], type='EVNT', project='mc16%')
                #print(List_ami)
                if List_ami:
                    for dataset in List_ami:
                        key = f"{EFT_op}_{EFT_type}_{process}_{decay}"
                        if key not in organized_info:
                            organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                        organized_info[key]['name'].append(dataset['ldn'])
                        organized_info[key]['cross_section'].append(dataset['cross_section'])

    # Construct the key from the EFT operator, process, and decay
    key_ = f"{EFT_op}_{EFT_type}_{proces}_{dec}"

    # Check if the key exists in the organized_info dictionary
    if key_ in organized_info:
        # If the key exists, return the first cross section
        xsection_nb = organized_info[key_]['cross_section'][0] # Return the cross section in nb
        xsection_fb = float(xsection_nb) * 1e6 # Convert the cross section to fb
        return xsection_fb
    else:
        # If the key does not exist, print an error message and return None
        #print(organized_info[key_]['cross_section'][0])
        print(f'Error: No datasets found for {key_}')
        return None


def take_xsec_fb(EFT_op, EFT_type, proces, dec):
    # Create the key
    key = f"{proces}_{dec}_{EFT_op}_{EFT_type}"
    
    # Open the text file and read the cross sections
    with open('VBS_xsection_test.txt', 'r') as f:
        for line in f:
            # Split the line into key and value
            key_file, xsection_fb = line.strip().split(': ')
            # If the key matches the input, return the cross section
            if key_file == key:
                return float(xsection_fb)
    
    # If the key was not found in the file, return None
    return None

def take_xsec_fb2(EFT_op, EFT_type, proces, dec):
    # Create the key
    key = f"{proces}_{dec}_{EFT_op}_{EFT_type}"
    
    # Open the text file and read the cross sections
    with open('VBS_xsection.txt', 'r') as f:
        for line in f:
            try:
                # Split the line into key and value
                key_file, xsection_fb_str = line.strip().split(': ')
                # If the key matches the input, attempt to convert the cross section to float
                if key_file == key:
                    xsection_fb = float(xsection_fb_str)
                    return xsection_fb
            except ValueError:
                # If conversion to float fails or the line does not conform to expected format
                print(f"No cross section found for {key_file}.")
                continue
    
    # If the key was not found in the file or the value is not a float, return None
    print(f"Key '{key}' not found or value is not a float.")
    


def plot_histograms(output_plot,desired_num_bins, file_path, label):
    # Open the first ROOT file to retrieve the list of parameters (branches)
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)
    first_root_file = ROOT.TFile(file_path, "READ")
    keys = [key.GetName() for key in first_root_file.GetListOfKeys()]
    print(keys)
    first_root_file.Close()

    for parameter_to_plot in keys:    
        # Create a new TCanvas
        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 800, 600)

        # Create a THStack
        hs = ROOT.THStack("hs", "Distribution of " + parameter_to_plot)
        
        
        # Loop over the files to retrieve and stack the histograms

        root_file = ROOT.TFile(file_path, "READ")
    
        # Retrieve the histogram
        histogram = root_file.Get(parameter_to_plot)
        histogram.SetDirectory(0)  # Detach from the file
        current_num_bins = histogram.GetNbinsX()
        
        rebin_factor = int(current_num_bins / desired_num_bins)
        if rebin_factor > 0:
            histogram.Rebin(rebin_factor)
            
         
        # Add the histogram to the stack
        hs.Add(histogram)
        


        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # x1, y1, x2, y2

        # Set legend style
        legend.SetFillColor(0) 
        legend.SetFillStyle(0)  
        legend.SetBorderSize(0)  
        legend.SetTextSize(0.03)  

        # Add entry to the legend
        legend.AddEntry(histogram, label, "l")
    
        # Close the ROOT file
        root_file.Close()
        
        # Draw the histograms
        #hs.Draw("nostack")
        hs.Draw("hist")
        hs.GetXaxis().SetTitle(parameter_to_plot )
        hs.GetYaxis().SetTitle("Events / bin")
        
        # Set log scale for y-axis
        canvas.SetLogy()
        
        # Draw the legend
        legend.Draw()
        canvas.Draw()
        
        # Update the canvas
        canvas.Update()
        
        # Save the canvas to a file
        canvas.SaveAs(output_plot+parameter_to_plot + "_hist.png")
        
def possible_process(proc, decay):
    
    valid_combinations = {
        "WmZ": ["llqq", "lvqq", "vvqq"],
        "WpZ": ["llqq", "lvqq", "vvqq"],
        "ZZ": ["llqq", "vvqq"],  # Only ll and vv decays are valid for ZZ
        "WmWm": ["lvqq"],
        "WpWm": ["lvqq"],
        "WpWp": ["lvqq"],
    } 

    return decay in valid_combinations.get(proc, []) 
            
def extract_counts_cutflow(dir_path):
    Efficiency_SR = {}
    Info_err_frac={}

    cutflow_merged_path = dir_path + "cutflow_merged.txt"
    cutflow_resolved_path = dir_path + "cutflow_resolved.txt"
    cutflow_frac_res_path = dir_path + "Info/frac_after_cuts_error_bar_resolved.txt"
    cutflow_frac_merged_path = dir_path + "Info/frac_after_cuts_error_bar_merged.txt"

    with open(cutflow_merged_path, 'r') as f:
        lines = f.readlines()
        Efficiency_SR["total"] = float(lines[2].split()[0])
        Efficiency_SR["merged"] = float(lines[-2].split()[2])

    with open(cutflow_resolved_path, 'r') as f:
        lines = f.readlines()
        Efficiency_SR["resolved"] = float(lines[-2].split()[2])

    with open(cutflow_frac_res_path, 'r') as f:
        Info_err_frac["err_frac_resolved"] = float(f.readline().split()[0])

    with open(cutflow_frac_merged_path, 'r') as f:
        Info_err_frac["err_frac_merged"] = float(f.readline().split()[0])

    return Efficiency_SR, Info_err_frac

def merge_cutflows(paths, out_file,process,decay):
    # Initialize a dictionary to hold the total counts
    total_counts = {}

    # Iterate over each path
    for path_ in paths:
        path = path_ + "cutflow_merged.txt"
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
    lines.append("                            Count    A_cumu    A_incr\n")

    # Initialize the previous count for calculating incremental percentages
    prev_count = total_count

    # For each cut name and total count
    for cut_name, count in total_counts.items():
        # Calculate the cumulative and incremental percentages
        cumu_percent = 100 * count / total_count
        incr_percent = 100 * count / prev_count
        
        cumul_percent = f"{cumu_percent:.1f}%"
        incre_percent = f"{incr_percent:.1f}%"

        # Add a line to the output list
        if cut_name == "Total":
            cut_name = ""
            incre_percent="-"
        lines.append(f"Pass {cut_name:<15}\t{count:>9}\t{cumul_percent:>8}\t{incre_percent:>6}\n")

        # Update the previous count
        prev_count = count
    out_file_ = out_file + "cutflow_merged_total.txt"
    # Write the lines to the output file
    with open(out_file_, 'w') as f:
        f.writelines(lines)
        
        
        
def combine_parts(Path,Parts, process, decay,Plot=False):
    path_parts= {}
    total_events = 0
    merged_events, resolved_events = 0, 0
    frac_err_merged_inv2, frac_err_resolved_inv2 = 0, 0
    
    working_parts = []
    working_path = {}
    for part in Parts:
        path = Path[f"{process}_{decay}"] + f"/{part}/DOCUT_YES/"
        path_parts[f"{process}_{decay}_{part}"] = path
        print(f"Path :{path}")
        if not (os.path.isfile(path + "/hists.root") and os.path.isfile(path + "cutflow_merged.txt") and os.path.isfile(path + "/ntuple_rivet.root") and 
                os.path.isfile(path + "Info/frac_after_cuts_error_bar_merged.txt")):
            continue
        working_parts.append(part)
        working_path[f"{process}_{decay}_{part}"]=path
        
        SR_eff, frac_err = extract_counts_cutflow(path)
        print(f"merged_events in {part}: {SR_eff['merged']}")
        total_events += SR_eff["total"]
        merged_events += SR_eff["merged"]
        resolved_events += SR_eff["resolved"]
        frac_err_merged_inv2 += frac_err["err_frac_merged"]**-2
        frac_err_resolved_inv2 += frac_err["err_frac_resolved"]**-2
    print(f"working parts: {working_parts}")
    if frac_err_merged_inv2 > 0:
        frac_err_merged = (frac_err_merged_inv2)**-0.5
    else:
        frac_err_merged = 0
        
    if frac_err_resolved_inv2 > 0:
        frac_err_resolved = (frac_err_resolved_inv2)**-0.5
    else:
        frac_err_resolved = 0

    # Create the combined directory if it doesn't exist
    combined_dir = Path[f"{process}_{decay}"] + "/combined/DOCUT_YES/"
    os.makedirs(combined_dir, exist_ok=True)
    merge_cutflows(list(working_path.values()), combined_dir,process,decay)

    with open(combined_dir + "/cross_section.txt", 'w') as f:
        f.write(str(take_xsec_fb("FM0","SM", process, decay)))

    with open(combined_dir + "/total_events.txt", 'w') as f:
        f.write(str(total_events))

    # Write each variable to a separate text file
    with open(combined_dir + "/merged_events.txt", 'w') as f:
        f.write(str(merged_events))

    with open(combined_dir + "/resolved_events.txt", 'w') as f:
        f.write(str(resolved_events))

    with open(combined_dir + "/frac_err_merged.txt", 'w') as f:
        f.write(str(frac_err_merged))

    with open(combined_dir + "/frac_err_resolved.txt", 'w') as f:
        f.write(str(frac_err_resolved))

    # Combine the hist.root files
    hist_files = [path + "/hists.root" for path in working_path.values()]
    subprocess.run(["hadd", "-f", combined_dir + "/hist.root"] + hist_files)

    if Plot:
        label_plot= f"{process}_{decay}_SM"
        plot_histograms(output_plot=combined_dir+ "/plots/" , desired_num_bins=200, file_path=combined_dir + "/hist.root", label=label_plot)

    # Combine the ntuple_rivet.root files
    ntuple_files = [path + "/ntuple_rivet.root" for path in working_path.values()]
    subprocess.run(["hadd", "-f", combined_dir + "/ntuple_rivet.root"] + ntuple_files)
    
    ntuple_truth_files = [path + "/ntuple_truth.root" for path in working_path.values() if os.path.isfile(path + "/ntuple_truth.root")]
    if ntuple_truth_files:
        subprocess.run(["hadd", "-f", combined_dir + "/ntuple_truth.root"] + ntuple_truth_files)
    else:
        print("No ntuple_truth.root files found.")
    # Create a backup directory with a unique name
    timestamp = time.strftime("%H%M")
    backup_dir = Path[f"{process}_{decay}"] + f"/Backup/combined_{total_events}_{timestamp}/"
    shutil.copytree(combined_dir, backup_dir)
 
 
 
    
def find_conf_dir(base_dir,conf):  
    if conf.startswith("user."):
        search_com= base_dir + f"/*{conf}*EXT0"
        print("searching for dir with pattern", search_com)
        conf_dir_arr = glob.glob(search_com)
        print("found possibilities for dir", conf_dir_arr)
        conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
        if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)

    else:
        pattern = r"MGPy8EG_aQGC(.*)_(.*)_1_(.*)_(.*)"
        match = re.search(pattern, conf)
        if(match):
            search_com = base_dir + f"/*{conf}*"
            print("searching for dir with pattern", search_com)
            conf_dir_arr = glob.glob(search_com)
            print("found possibilities for dir", conf_dir_arr)
            conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
            if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)
        

    return conf_dir 
    
    
    
def path_ntuple(dir_details,Processes,Decays,Operators):
    Path= {}
    for process in Processes:
        for decay in Decays:
            for op in Operators:
                proc_dec = f"{process}_{decay}"
                base_dir= f"/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/{proc_dec}/"
                if op=="SM":
                    conf_= f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}"
                else:
                    conf_= f"MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}"   
                conf_dir= find_conf_dir(base_dir,conf_)

                path_ = conf_dir + f"/DOCUT_YES/{dir_details}/"
                if op=="SM":
                    path_ntuple = path_ + "/combined"
                    #path_ntuple = base_dir + f"mc16_13TeV.503172.MGPy8EG_aQGCFM0_SM_1_ZZ_llqq.merge.EVNT.e8263_e7400/DOCUT_YES/Tables/BDT_ntuple_01/combined/ntuple_rivet.root"
                else:  
                    path_ntuple = path_ 
                     
                if os.path.exists(path_ntuple):
                    Path[f"{process}_{decay}_{op}"] = path_ntuple
                else:
                    print(f"For process {proc_dec}, the path does not exist.")
    return Path

def write_xsex_merged(dir_details,Processes,Decays,Operators,output_dir):
    Paths = path_ntuple(dir_details,Processes,Decays,Operators)
    with open(output_dir + "VBS_xsection_merged_fb.txt", "w") as f:
        for process in Processes:
            for decay in Decays:
                for op in Operators:
                    key = f"{process}_{decay}_{op}"
                    if key in Paths:
                        path= Paths[key]
                        if op == "SM":
                            total_events_path = path + "/total_events.txt"
                            merge_events_path = path + "/merged_events.txt"
                            xsection_path = path + "/cross_section.txt"
                            if os.path.exists(total_events_path):
                                with open(total_events_path, 'r') as file:
                                    total_events = float(file.read().strip())
                            else:
                                total_events = None  # or handle the missing file appropriately

                            if os.path.exists(merge_events_path):
                                with open(merge_events_path, 'r') as file:
                                    merged_events = float(file.read().strip())
                            else:
                                merged_events = None  # or handle the missing file appropriately

                            if os.path.exists(xsection_path):
                                with open(xsection_path, 'r') as file:
                                    xsection = float(file.read().strip())
                                if total_events is not None and merged_events is not None:
                                    xsection_merged = xsection * (merged_events / total_events)
                                else:
                                    xsection_merged = None  # Handle division by zero or missing values appropriately
                        else:
                            xsection_merged_path = path+ "/Info/xsec_times_frac_merged_fb.txt"
                            #print(f"Path for xsection: {xsection_merged_path}")
                            if os.path.exists(xsection_merged_path):
                                with open(xsection_merged_path, 'r') as file:
                                    xsection_merged = float(file.read().strip())
                                    #print(f"Xsection merged from path: {xsection_merged}")
                            else:
                                xsection_merged = None
                                print(f"No path found for {key}")
                                
                        if xsection_merged is not None:
                            f.write(f"{key}: {xsection_merged}\n")
                        else:
                            f.write(f"{key}: Data not available\n")                            

                    else:
                        print(f"No path found for {key}")

def combine_ntuples(dir_details, Processes, Decays, Operators,output_dir,truth=False):
    # Get the paths of the ntuple files
    paths = path_ntuple(dir_details, Processes, Decays, Operators)

    # Define the directory where the files will be combined
    combined_dir = output_dir
    individual_dir = output_dir + "individual/"


    for process in Processes:
        for decay in Decays:
            proc_decay= f"{process}_{decay}"
            
            # Create the combined directory if it doesn't exist
            if not os.path.exists(combined_dir):
                os.makedirs(combined_dir)

            # Combine the ntuple files using hadd
            ntuple_files = [value + "/ntuple_rivet.root" for key, value in paths.items() if proc_decay in key]
            subprocess.run(["hadd", "-f", combined_dir + f"/{proc_decay}_ntuple_rivet.root"] + ntuple_files)

            if not os.path.exists(individual_dir):
                os.makedirs(individual_dir)
                
            for op in Operators:
                for key, value in paths.items():
                    if f"{proc_decay}_{op}" in key:
                        src_file = value + "/ntuple_rivet.root"
                        dest_file = individual_dir + f"/{proc_decay}_{op}_ntuple_rivet.root"
                        shutil.copy(src_file, dest_file)

            if truth:
                ntuple_truth_files = [value + "/ntuple_truth.root" for key, value in paths.items() if proc_decay in key and os.path.isfile(value + "/ntuple_truth.root")]
                if ntuple_truth_files:
                    subprocess.run(["hadd", "-f", combined_dir + f"/{proc_decay}_ntuple_truth.root"] + ntuple_truth_files)
                else:
                    print(f"No ntuple_truth.root files found for {proc_decay}")

            # Write the paths to a text file
            with open(combined_dir + f"/{proc_decay}_paths.txt", "w") as file:
                for path in ntuple_files:
                    file.write(path + "\n")

            print(f"All ntuple files have been combined in the directory: {combined_dir}")
    
    
    ntuple_files_all = [values + "/ntuple_rivet.root" for values in paths.values()]
    subprocess.run(["hadd", "-f", combined_dir + f"/All_ntuple_rivet.root"] + ntuple_files_all)
    
    with open(combined_dir + f"/All_paths.txt", "w") as file:
        for path in ntuple_files_all:
            file.write(path + "\n")

                    



