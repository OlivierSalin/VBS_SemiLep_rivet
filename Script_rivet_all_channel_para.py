import subprocess
import os
from multiprocessing import Pool
import lib_utils as lu
import utils_func as uf
import itertools
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--nb_part", default = 0)
parser.add_option("--nb_events", default = 20000)
parser.add_option("--aQGC", default= True)
parser.add_option("--Name", default= "All_ops")
parser.add_option("--All_channel", default = False)
parser.add_option("--nb_lep", default = 1)

opts, _ = parser.parse_args()

nb_lepton= int(opts.nb_lep)
all_ops = ["SM"]
#all_ops=["FM0","FM2"]
processes = ["WpZ","ZZ"]
processes = ["WmZ"]
decays = ["llqq"]

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

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02","FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
#all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
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
        
def run_command_SM(proc_part_tuple):
    process, part = proc_part_tuple
    for decay in decays:
        if uf.possible_process(process,decay):
            if part == "":
                command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES","--part",f"{part}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
            else:
                command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES","--part",f"{part}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
            subprocess.run(command)
            
def run_command_aQGC(proc_op_tuple):
        process, op = proc_op_tuple
        for decay in decays:
            if uf.possible_process(process,decay):
                if op == "SM":
                    command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_{decay}", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
                else:
                    command = ["python", "run_rivet.py", "--evtMax", f"{number_events}", "--conf", f"MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{opts.Name}"]
                subprocess.run(command)


combined_processes_parts = list(itertools.product(processes, Parts))
#print(combined_processes_parts)

# To combine processes with all_ops_cat and apply a hypothetical function
combined_processes_ops = list(itertools.product(processes, all_ops_cat))
#print(combined_processes_ops)

print(len(combined_processes_parts)+len(combined_processes_ops))


if __name__ == "__main__":
    with Pool() as p:
        process_part_tuples = list(itertools.product(processes, Parts))
        p.map(run_command_SM, process_part_tuples)
        if opts.aQGC:
            process_op_tuples = list(itertools.product(processes, all_ops_cat))
            p.map(run_command_aQGC, process_op_tuples)
    
    if nb_part != 0:    
        for process in processes:
            for decay in decays:
                uf.combine_parts(Path_output,Parts,process,decay,Plot=False)
                


            
        #p.map(run_command_SM, processes_SM)
        
    
        
        
