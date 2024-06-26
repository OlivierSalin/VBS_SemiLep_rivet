import os
import shutil
import glob
import re




all_ops = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_SM = ["FM0"]
processes=["WpZ","ZZ"]
processes=["WpZ"]
decays=["llqq"]



def file_for_table(mydir,op):
    run_dir = mydir + "Tables/Spin/Polar/"
    #run_dir = mydir + "Tables/BDT_sample/BDT_sample_40000/"
    os.makedirs(run_dir, exist_ok=True)

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
        shutil.copy2(mydir + file, run_dir)

    # Copy the repository mydir + "/plots/"
    if os.path.exists(run_dir + "/Info"):
        shutil.rmtree(run_dir + "/Info")
    shutil.copytree(mydir + "/Info", run_dir + "/Info")
    if op=="SM":
        my_dir_comb= mydir.replace("DOCUT_YES/", "combined/DOCUT_YES/")
        if os.path.exists(my_dir_comb):
            if os.path.exists(run_dir + "combined/"):
                shutil.rmtree(run_dir + "combined/")
            shutil.copytree(my_dir_comb, run_dir + "combined/")
    
base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/"

for process in processes:
    for decay in decays:
        for operator in all_ops_cat:
            # Construct the directory pattern
            if operator == "SM":
                dir_pattern = os.path.join(base_dir, f"{process}_{decay}/user.osalin.MadGraph_{process}_{decay}_FM0_SM*/DOCUT_YES/")
            else:
                dir_pattern = os.path.join(base_dir, f"{process}_{decay}/user.osalin.MadGraph_{process}_{decay}_{operator}_QUAD*/DOCUT_YES/")
            print(dir_pattern)
            # Use glob to find the directories that match the pattern
            matching_dirs = glob.glob(dir_pattern)
            print(matching_dirs)
            
            # Call the function for each matching directory
            for dir_path in matching_dirs:
                file_for_table(dir_path,operator)