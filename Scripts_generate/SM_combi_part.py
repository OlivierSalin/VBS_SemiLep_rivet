import subprocess
import os
from multiprocessing import Pool
import lib_utils as lu
import utils_func as uf
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--nb_part", default = 0)
opts, _ = parser.parse_args()

all_ops = ["SM"]
#all_ops=["FM0","FM2"]
processes = ["WpZ","ZZ"]
decays = ["llqq"]

nb_part = int(opts.nb_part)

if nb_part == 0:
    Parts=[""]
else:
    Parts = [f"part_{i+1}" for i in range(nb_part)]
    

path_base= "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"
Path_output={}
number_events = 100
total_events = number_events * nb_part

for process in processes:
    for decay in decays:
        proc_dec = f"{process}_{decay}"
        base_dir= f"/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/{proc_dec}/"
        conf_= f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}"
        conf_dir, _, _ = lu.find_evnt_dir_and_file_bis(base_dir,conf_)
        Path_output[proc_dec] = conf_dir
        
path_build_command = "-I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/ `root-config --cflags --libs`"
for process in processes:
    for decay in decays:
        build_command = ["rivet-build", f"Rivet{process}_{decay}.so", f"{process}_{decay}.cc", f"{path_build_command}"]
        #subprocess.run(build_command)
        
def run_command_SM(part):
    for process in processes:
        for decay in decays:
            if part == "":
                command = ["python", "run_rivet.py", "--evtMax", "10000", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES","--part",f"{part}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
            else:
                command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES","--part",f"{part}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
            subprocess.run(command)

""" all_ops_SM=["FM0"]
processes_SM = ["WmZ", "WpZ","ZZ"]        
def run_command_SM(process):
    command = ["python", "run_rivet.py", "--evtMax", "10000", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_llqq", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
    subprocess.run(command) """

if __name__ == "__main__":
    #with Pool() as p:
        #p.map(run_command_SM, Parts)
        
    for process in processes:
        for decay in decays:
            uf.combine_parts(Path_output,Parts,process,decay,Plot=False)
            
        #p.map(run_command_SM, processes_SM)
        
    
        
        
