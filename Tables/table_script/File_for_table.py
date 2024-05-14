import os
import shutil
import glob
import re




all_ops = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
all_ops_SM = ["FM0"]
processes=["WmZ","WpZ","ZZ"]
decays=["llqq"]



def file_for_table(mydir):
    run_dir = mydir + "Tables/Tables_file_first/"
    os.makedirs(run_dir, exist_ok=True)

    files_to_copy = [
        "hists.root",
        "frac_after_cuts_error_bar_merged.txt",
        "frac_after_cuts_error_bar_resolved.txt",
        "cutflow_merged.txt",
        "cutflow_resolved.txt",
        "cutflow_merged_img.png",
        "cutflow_resolved_img.png",
        "Cross_section_fb.txt",
    ]

    for file in files_to_copy:
        shutil.copy2(mydir + file, run_dir)

    # Copy the repository mydir + "/plots/"
    shutil.copytree(mydir + "/plots/", run_dir + "plots/")
    
base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"

for process in processes:
    for decay in decays:
        for operator in all_ops:
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
                file_for_table(dir_path)