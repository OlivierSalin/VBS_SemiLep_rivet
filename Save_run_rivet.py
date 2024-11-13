import os
import shutil
import glob
import re
import utils_func as uf
from optparse import OptionParser
parser = OptionParser()

parser.add_option("--All_channel", default = False)
parser.add_option("--nb_lep", default = 2)
parser.add_option("--Name", default= "test")
opts, _ = parser.parse_args()

print(f'opts.All_channel: {opts.All_channel}, opts.nb_lep: {opts.nb_lep}, opts.Name: {opts.Name}')

all_ops = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_SM = ["FM0"]

Processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
Decay = ["llqq",'lvqq','vvqq']

nb_lepton= int(opts.nb_lep)


all_ops_cat = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02","FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

all_ops_cat=["FM0","FM2","FS1","FT1","FT5"]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/"


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
    

def file_for_table(mydir, op,process,decay):
    key= f"{process}_{decay}"
    bases = base_dir
    run_dir = mydir + f"Tables/Vhad_first1/{opts.Name}/"
    if op == "SM":
        operator= "FM0_SM"
    else: 
        operator = f"{op}_QUAD"
    run_dir2 = bases + f"Tables/Vhad_first1/{opts.Name}/{opts.nb_lep}Lepton/{key}/{operator}/"
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(run_dir2, exist_ok=True)

    files_to_copy = [
        "hists.root",
        "ntuple_rivet.root",
        #"ntuple_truth.root",
        "Info/frac_cuts_merged.txt",
        "Info/frac_after_cuts_error_bar_merged.txt",
        "Info/frac_cuts_resolved.txt",
        "Info/frac_after_cuts_error_bar_resolved.txt",
        "cutflow_merged.txt",
        "cutflow_resolved.txt",
        "cutflow_merged_img.png",
        "cutflow_resolved_img.png",
        "Info/Cross_section_fb.txt",
    ]

    for file in files_to_copy:
        source_file = os.path.join(mydir, file)
        if os.path.exists(source_file):
            if "Info/" in file:
                file = re.sub(r"Info/", "", file)
            destination_file = os.path.join(run_dir, file)
            destination_file2 = os.path.join(run_dir2, file)
            try:
                shutil.copy2(source_file, destination_file)
                shutil.copy2(source_file, destination_file2)
            except FileNotFoundError:
                print(f"Warning: {source_file} not found, skipping.")

    # Copy the repository mydir + "/Info/"
    info_dir = os.path.join(run_dir, "Info")
    info_dir2 = os.path.join(run_dir2, "Info")
    if os.path.exists(info_dir):
        try:
            shutil.rmtree(info_dir)
        except Exception as e:
            print(f"Warning: Unable to delete {info_dir}. Error: {e}")
    try:
        shutil.copytree(os.path.join(mydir, "Info"), info_dir)
    except FileNotFoundError:
        print(f"Warning: {os.path.join(mydir, 'Info')} not found, skipping.")
        
    if os.path.exists(info_dir2):
        try:
            shutil.rmtree(info_dir2)
        except Exception as e:
            print(f"Warning: Unable to delete {info_dir2}. Error: {e}")
    try:
        shutil.copytree(os.path.join(mydir, "Info"), info_dir2)
    except FileNotFoundError:
        print(f"Warning: {os.path.join(mydir, 'Info')} not found, skipping.")

    if op == "SM":
        my_dir_comb = mydir.replace("DOCUT_YES/", "combined/DOCUT_YES/")
        combined_dir = os.path.join(run_dir, "combined/")
        combined_dir2 = os.path.join(run_dir2, "combined/")
        os.makedirs(combined_dir2, exist_ok=True)
        if os.path.exists(my_dir_comb):
            if os.path.exists(combined_dir):
                try:
                    shutil.rmtree(combined_dir)
                except Exception as e:
                    print(f"Warning: Unable to delete {combined_dir}. Error: {e}")
            try:
                shutil.copytree(my_dir_comb, combined_dir)
            except FileNotFoundError:
                print(f"Warning: {my_dir_comb} not found, skipping.")

        if os.path.exists(my_dir_comb):
            if os.path.exists(combined_dir2):
                try:
                    shutil.rmtree(combined_dir2)
                except Exception as e:
                    print(f"Warning: Unable to delete {combined_dir2}. Error: {e}")
            try:
                shutil.copytree(my_dir_comb, combined_dir2)
            except FileNotFoundError:
                print(f"Warning: {my_dir_comb} not found, skipping.")
#base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"

for process in processes:
    for decay in decays:
        if uf.possible_process(process,decay):
            for operator in all_ops_cat:
                # Construct the directory pattern
                if operator == "SM":
                    dir_pattern = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_SM_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/")
                else:
                    dir_pattern = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{operator}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/")
                print(dir_pattern)
                # Use glob to find the directories that match the pattern
                matching_dirs = glob.glob(dir_pattern)
                print(matching_dirs)
                
                
                # Call the function for each matching directory
                for dir_path in matching_dirs:
                    file_for_table(dir_path,operator,process,decay)