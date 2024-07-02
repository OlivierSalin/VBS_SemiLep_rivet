import subprocess
import os
from multiprocessing import Pool
import lib_utils as lu
import utils_func as uf
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--nb_part", default = 0)
parser.add_option("--nb_events", default = 20000)
parser.add_option("--aQGC", default= True)
parser.add_option("--Name", default= "All_ops")

opts, _ = parser.parse_args()

all_ops = ["SM"]
#all_ops=["FM0","FM2"]
processes = ["WpZ","ZZ"]
processes = ["WmZ"]
decays = ["llqq"]

processes = ["WmZ","WpZ","ZZ"]
decays = ['vvqq']

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02","FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
nb_part = int(opts.nb_part)

if nb_part == 0:
    Parts=[""]
else:
    Parts = [f"part_{i+1}" for i in range(nb_part)]
    

path_base= "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"
Path_output={}
number_events = int(opts.nb_events)
total_events = number_events * nb_part

for process in processes:
    for decay in decays:
        proc_dec = f"{process}_{decay}"
        base_dir= f"/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/{proc_dec}/"
        conf_= f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}"
        conf_dir, _, _ = lu.find_evnt_dir_and_file_bis(base_dir,conf_)
        Path_output[proc_dec] = conf_dir
        for part in Parts:
            path = Path_output[f"{process}_{decay}"] + f"/{part}/DOCUT_YES/"
            if not os.path.exists(path):
                os.makedirs(path)
        
path_build_command = "-I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/ `root-config --cflags --libs`"
for process in processes:
    for decay in decays:
        build_command = ["rivet-build", f"Rivet{process}_{decay}.so", f"{process}_{decay}.cc", f"{path_build_command}"]
        subprocess.run(build_command)
        
def run_command_SM(part):
    for process in processes:
        for decay in decays:
            if uf.possible_process(process,decay):
                if part == "":
                    command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES","--part",f"{part}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
                else:
                    command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES","--part",f"{part}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
                subprocess.run(command)
            
def run_command_aQGC(op):
    for process in processes:
        for decay in decays:
            if uf.possible_process(process,decay):
                if op == "SM":
                    command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
                else:
                    command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
                subprocess.run(command)

""" all_ops_SM=["FM0"]
processes_SM = ["WmZ", "WpZ","ZZ"]        
def run_command_SM(process):
    command = ["python", "run_rivet.py", "--evtMax", "10000", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_llqq", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
    subprocess.run(command) """

if __name__ == "__main__":

    
    if nb_part != 0:    
        for process in processes:
            for decay in decays:
                uf.combine_parts(Path_output,Parts,process,decay,Plot=True)
            
        #p.map(run_command_SM, processes_SM)
        
    
        
        
